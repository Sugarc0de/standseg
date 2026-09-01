//! GeoPackage output.
//!
//! A GeoPackage is a SQLite database with an agreed set of metadata tables, so
//! this is the one output format the program does not write byte by byte. The
//! SQLite it links is compiled from source into the binary, the same arrangement
//! the GeoTIFF ZSTD support already uses: still one executable, still no system
//! libraries, still nothing to install at the far end.
//!
//! Chosen over a shapefile because a shapefile is five files, a 2 GB ceiling and
//! ten-character field names, and over GeoJSON because RFC 7946 wants WGS 84 and
//! forest imagery is projected. QGIS, ArcGIS, R (`sf`) and Python (`geopandas`,
//! `fiona`) all open a GeoPackage directly.
//!
//! The file is deterministic. `gpkg_contents.last_change` is required by the
//! specification and is the only clock-dependent field in the format, so it is
//! pinned to a constant rather than filled from `now` -- the same command twice
//! gives the same bytes, which is the rule everywhere else in this program.

use std::path::Path;

use rusqlite::{params, Connection};

use crate::geo::{self, Transform};
use crate::io::{IoError, Provenance, Result};
use crate::vector::{self, RegionPolygons};

/// `0x47504B47`, ASCII "GPKG": the `application_id` every reader checks.
const APPLICATION_ID: i32 = 0x4750_4B47;
/// GeoPackage 1.4.0, as `user_version`.
const USER_VERSION: i32 = 10400;

/// See the module note on determinism.
const FIXED_LAST_CHANGE: &str = "1970-01-01T00:00:00.000Z";

/// Where a CRS that is described only by a WKT string, with no EPSG code to
/// call it by, gets filed. Above the reserved range and above anything EPSG
/// issues.
const CUSTOM_SRS_ID: i32 = 100_000;

fn sql(e: rusqlite::Error, path: &Path) -> IoError {
    IoError::new(format!("can't write {}: {e}", path.display()))
}

/// Write a region map as polygons.
///
/// `skip` is the region id that stands for nodata, or `None` to keep every
/// region -- see [`crate::vector::polygonize`].
#[allow(clippy::too_many_arguments)]
pub fn write_region_map(
    path: &Path,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
    skip: Option<u32>,
    transform: Option<Transform>,
    epsg: Option<u32>,
    wkt: Option<&str>,
    layer: &str,
    prov: &Provenance,
) -> Result<usize> {
    let t = transform.unwrap_or(geo::PIXEL_SPACE);
    let polys = vector::polygonize(rband, nlines, nsamps, skip, &t);

    // Pixel area is exact for a rectilinear region map, so a stand's area is
    // its pixel count times this -- no shoelace, no rounding drift, and it
    // agrees with the region map to the pixel.
    let cell = transform.map(|t| geo::pixel_area(&t));

    let (srs_id, org, org_id, definition) = match (epsg, wkt) {
        (Some(code), w) => (
            code as i32,
            "EPSG",
            code as i32,
            geo::wkt_for_epsg(code)
                .or_else(|| w.map(str::to_string))
                .unwrap_or_else(|| "undefined".into()),
        ),
        (None, Some(w)) => (CUSTOM_SRS_ID, "NONE", CUSTOM_SRS_ID, w.to_string()),
        // No georeferencing at all: the coordinates are pixels, and the file
        // says so rather than implying a place.
        (None, None) => (-1, "NONE", -1, "undefined".into()),
    };

    if path.exists() {
        std::fs::remove_file(path)
            .map_err(|e| IoError::new(format!("can't replace {}: {e}", path.display())))?;
    }
    let conn = Connection::open(path).map_err(|e| sql(e, path))?;
    write_all(
        &conn, path, &polys, layer, srs_id, org, org_id, &definition, cell, prov,
    )?;
    conn.close()
        .map_err(|(_, e)| sql(e, path))
        .map(|()| polys.len())
}

#[allow(clippy::too_many_arguments)]
fn write_all(
    conn: &Connection,
    path: &Path,
    polys: &[RegionPolygons],
    layer: &str,
    srs_id: i32,
    org: &str,
    org_id: i32,
    definition: &str,
    cell: Option<f64>,
    prov: &Provenance,
) -> Result<()> {
    let go = |s: &str| conn.execute_batch(s).map_err(|e| sql(e, path));

    go(&format!(
        "PRAGMA application_id = {APPLICATION_ID};
         PRAGMA user_version = {USER_VERSION};
         PRAGMA journal_mode = DELETE;"
    ))?;

    go("CREATE TABLE gpkg_spatial_ref_sys (
            srs_name TEXT NOT NULL,
            srs_id INTEGER NOT NULL PRIMARY KEY,
            organization TEXT NOT NULL,
            organization_coordsys_id INTEGER NOT NULL,
            definition TEXT NOT NULL,
            description TEXT
        );
        CREATE TABLE gpkg_contents (
            table_name TEXT NOT NULL PRIMARY KEY,
            data_type TEXT NOT NULL,
            identifier TEXT UNIQUE,
            description TEXT DEFAULT '',
            last_change DATETIME NOT NULL,
            min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE,
            srs_id INTEGER,
            CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id)
                REFERENCES gpkg_spatial_ref_sys(srs_id)
        );
        CREATE TABLE gpkg_geometry_columns (
            table_name TEXT NOT NULL,
            column_name TEXT NOT NULL,
            geometry_type_name TEXT NOT NULL,
            srs_id INTEGER NOT NULL,
            z TINYINT NOT NULL,
            m TINYINT NOT NULL,
            CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name),
            CONSTRAINT uk_gc_table_name UNIQUE (table_name),
            CONSTRAINT fk_gc_tn FOREIGN KEY (table_name)
                REFERENCES gpkg_contents(table_name),
            CONSTRAINT fk_gc_srs FOREIGN KEY (srs_id)
                REFERENCES gpkg_spatial_ref_sys (srs_id)
        );")?;

    // The three rows the specification requires of every GeoPackage, then ours
    // if it is not already one of them.
    let ins = |id: i32, name: &str, o: &str, oid: i32, def: &str| {
        conn.execute(
            "INSERT OR IGNORE INTO gpkg_spatial_ref_sys
                 (srs_id, srs_name, organization, organization_coordsys_id, definition)
             VALUES (?1, ?2, ?3, ?4, ?5)",
            params![id, name, o, oid, def],
        )
        .map(|_| ())
        .map_err(|e| sql(e, path))
    };
    ins(-1, "Undefined cartesian SRS", "NONE", -1, "undefined")?;
    ins(0, "Undefined geographic SRS", "NONE", 0, "undefined")?;
    ins(
        4326,
        "WGS 84 geodetic",
        "EPSG",
        4326,
        &geo::wkt_for_epsg(4326).unwrap(),
    )?;
    ins(srs_id, &format!("{org}:{org_id}"), org, org_id, definition)?;

    // `area` is left out entirely when there is no transform: a column of pixel
    // counts relabelled as area is worse than no column.
    let area_col = if cell.is_some() { ", area REAL" } else { "" };
    go(&format!(
        "CREATE TABLE \"{layer}\" (
             fid INTEGER PRIMARY KEY AUTOINCREMENT,
             geom MULTIPOLYGON,
             region_id INTEGER NOT NULL,
             n_pixels INTEGER NOT NULL{area_col}
         );"
    ))?;

    let (mut nx, mut ny) = (f64::INFINITY, f64::INFINITY);
    let (mut xx, mut xy) = (f64::NEG_INFINITY, f64::NEG_INFINITY);

    conn.execute_batch("BEGIN").map_err(|e| sql(e, path))?;
    {
        let stmt_sql = if cell.is_some() {
            format!("INSERT INTO \"{layer}\" (geom, region_id, n_pixels, area) VALUES (?1,?2,?3,?4)")
        } else {
            format!("INSERT INTO \"{layer}\" (geom, region_id, n_pixels) VALUES (?1,?2,?3)")
        };
        let mut stmt = conn.prepare(&stmt_sql).map_err(|e| sql(e, path))?;
        for p in polys {
            let (a, b, c, d) = p.envelope();
            nx = nx.min(a);
            ny = ny.min(b);
            xx = xx.max(c);
            xy = xy.max(d);

            let blob = geometry_blob(p, srs_id);
            match cell {
                Some(cell) => stmt.execute(params![
                    blob,
                    p.id as i64,
                    p.npix as i64,
                    p.npix as f64 * cell
                ]),
                None => stmt.execute(params![blob, p.id as i64, p.npix as i64]),
            }
            .map_err(|e| sql(e, path))?;
        }
    }
    conn.execute_batch("COMMIT").map_err(|e| sql(e, path))?;

    if polys.is_empty() {
        nx = 0.0;
        ny = 0.0;
        xx = 0.0;
        xy = 0.0;
    }

    // The command that produced the file, in the place a GeoPackage reader will
    // actually show it -- the layer description. Same role as ENVI `history`
    // and TIFF `ImageDescription`.
    let description = if prov.command.is_empty() {
        String::new()
    } else {
        format!("{} -- {}", prov.software, prov.command)
    };
    conn.execute(
        "INSERT INTO gpkg_contents
             (table_name, data_type, identifier, description, last_change,
              min_x, min_y, max_x, max_y, srs_id)
         VALUES (?1, 'features', ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
        params![layer, description, FIXED_LAST_CHANGE, nx, ny, xx, xy, srs_id],
    )
    .map_err(|e| sql(e, path))?;

    conn.execute(
        "INSERT INTO gpkg_geometry_columns
             (table_name, column_name, geometry_type_name, srs_id, z, m)
         VALUES (?1, 'geom', 'MULTIPOLYGON', ?2, 0, 0)",
        params![layer, srs_id],
    )
    .map_err(|e| sql(e, path))?;

    Ok(())
}

/// A GeoPackage geometry blob: the "GP" header, an envelope, then standard WKB.
fn geometry_blob(p: &RegionPolygons, srs_id: i32) -> Vec<u8> {
    let mut v = Vec::with_capacity(64 + p.parts.iter().flatten().map(|r| r.len() * 16).sum::<usize>());
    v.extend_from_slice(b"GP");
    v.push(0); // version 0, i.e. GeoPackage 1.x binary
    // bit 0 set: little-endian. bits 1-3 = 1: an envelope of [x, y] follows.
    v.push(0b0000_0011);
    v.extend_from_slice(&srs_id.to_le_bytes());
    let (nx, ny, xx, xy) = p.envelope();
    for f in [nx, xx, ny, xy] {
        v.extend_from_slice(&f.to_le_bytes());
    }

    // WKB MultiPolygon, so single- and multi-part stands share one type and
    // readers never have to cope with a mixed-geometry layer.
    v.push(1);
    v.extend_from_slice(&6u32.to_le_bytes());
    v.extend_from_slice(&(p.parts.len() as u32).to_le_bytes());
    for part in &p.parts {
        v.push(1);
        v.extend_from_slice(&3u32.to_le_bytes());
        v.extend_from_slice(&(part.len() as u32).to_le_bytes());
        for ring in part {
            v.extend_from_slice(&(ring.len() as u32).to_le_bytes());
            for &(x, y) in ring {
                v.extend_from_slice(&x.to_le_bytes());
                v.extend_from_slice(&y.to_le_bytes());
            }
        }
    }
    v
}
