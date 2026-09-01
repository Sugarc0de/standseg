//! GeoPackage output, end to end.
//!
//! The claim being tested is not "a file appeared". It is that the polygons are
//! the region map: every pixel covered exactly once, by the polygon of the
//! region that pixel actually belongs to. That is checked by rasterising the
//! polygons back and comparing to the raster the same run wrote, which is the
//! same standard the rest of this repository holds itself to -- the vector
//! output is either identical to the region map or it is wrong.
//!
//! Output goes to `build/out/`, never into `tests/`.

use std::path::{Path, PathBuf};
use std::process::Command;

use rusqlite::Connection;

const BIN: &str = env!("CARGO_BIN_EXE_standseg");

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn outdir(name: &str) -> PathBuf {
    let d = root().join("build/out/gpkg").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn ok(args: &[&str]) -> std::process::Output {
    let o = Command::new(BIN)
        .current_dir(root())
        .args(args)
        .output()
        .expect("run standseg");
    assert!(
        o.status.success(),
        "standseg {args:?} failed:\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
    o
}

/// The 250 x 250 four-band scene the quick start uses. It carries real UTM
/// zone 15N georeferencing, which is what makes it worth testing against.
const DEMO: &str = "tests/golden/misc/temp_byte_bip";
const DEMO_ARGS: &[&str] = &["-t", "10", "-m", ".1", "-n", "15,15,100,2500,2500"];

/// Read an ENVI region map and its header's grid.
fn read_envi(base: &Path) -> (Vec<u32>, usize, usize) {
    // Same rule as `io::envi::header_path`: replace the extension, else append.
    let replaced = base.with_extension("hdr");
    let hdr = std::fs::read_to_string(&replaced)
        .or_else(|_| std::fs::read_to_string(format!("{}.hdr", base.display())))
        .unwrap_or_else(|e| panic!("no header for {}: {e}", base.display()));
    let field = |k: &str| -> usize {
        hdr.lines()
            .find(|l| l.trim_start().starts_with(k))
            .and_then(|l| l.split('=').nth(1))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or_else(|| panic!("no {k} in header"))
    };
    let (ns, nl) = (field("samples"), field("lines"));
    let nbytes = match field("data type") {
        1 => 1,
        12 => 2,
        13 => 4,
        t => panic!("unexpected data type {t}"),
    };
    let raw = std::fs::read(base).expect("region map");
    let band = raw
        .chunks_exact(nbytes)
        .map(|c| {
            let mut b = [0u8; 4];
            b[..nbytes].copy_from_slice(c);
            u32::from_le_bytes(b)
        })
        .collect();
    (band, nl, ns)
}

struct Feature {
    region_id: u32,
    n_pixels: u64,
    area: Option<f64>,
    parts: Vec<Part>,
    envelope: [f64; 4],
    srs_id: i32,
}

/// Rings of one part: the shell, then any holes.
type Part = Vec<Vec<(f64, f64)>>;

/// Decode a GeoPackage geometry blob back to rings. Deliberately written
/// against the specification rather than by reusing the writer, so a mistake
/// shared by both would have to be made twice.
fn decode(blob: &[u8]) -> (i32, [f64; 4], Vec<Part>) {
    assert_eq!(&blob[..2], b"GP", "magic");
    assert_eq!(blob[2], 0, "version");
    assert_eq!(blob[3] & 1, 1, "little endian");
    assert_eq!((blob[3] >> 1) & 0x07, 1, "envelope indicator");
    let srs = i32::from_le_bytes(blob[4..8].try_into().unwrap());

    let f64_at = |o: usize| f64::from_le_bytes(blob[o..o + 8].try_into().unwrap());
    let u32_at = |o: usize| u32::from_le_bytes(blob[o..o + 4].try_into().unwrap());

    let env = [f64_at(8), f64_at(16), f64_at(24), f64_at(32)];
    let mut o = 40;
    assert_eq!(blob[o], 1);
    assert_eq!(u32_at(o + 1), 6, "MultiPolygon");
    let npoly = u32_at(o + 5);
    o += 9;

    let mut parts = Vec::new();
    for _ in 0..npoly {
        assert_eq!(blob[o], 1);
        assert_eq!(u32_at(o + 1), 3, "Polygon");
        let nring = u32_at(o + 5);
        o += 9;
        let mut rings = Vec::new();
        for _ in 0..nring {
            let npt = u32_at(o) as usize;
            o += 4;
            let mut ring = Vec::with_capacity(npt);
            for i in 0..npt {
                ring.push((f64_at(o + i * 16), f64_at(o + i * 16 + 8)));
            }
            o += npt * 16;
            rings.push(ring);
        }
        parts.push(rings);
    }
    assert_eq!(o, blob.len(), "trailing bytes");
    (srs, env, parts)
}

fn features(gpkg: &Path) -> (Vec<Feature>, String, i32) {
    let conn = Connection::open(gpkg).expect("open gpkg");
    let (layer, srs_id): (String, i32) = conn
        .query_row(
            "SELECT table_name, srs_id FROM gpkg_geometry_columns",
            [],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )
        .expect("geometry_columns");
    let has_area = conn
        .prepare(&format!("SELECT * FROM \"{layer}\" LIMIT 0"))
        .unwrap()
        .column_names()
        .contains(&"area");

    let sql = if has_area {
        format!("SELECT region_id, n_pixels, area, geom FROM \"{layer}\" ORDER BY fid")
    } else {
        format!("SELECT region_id, n_pixels, NULL, geom FROM \"{layer}\" ORDER BY fid")
    };
    let mut stmt = conn.prepare(&sql).unwrap();
    let rows: Vec<Feature> = stmt
        .query_map([], |r| {
            let blob: Vec<u8> = r.get(3)?;
            let (srs, envelope, parts) = decode(&blob);
            Ok(Feature {
                region_id: r.get::<_, i64>(0)? as u32,
                n_pixels: r.get::<_, i64>(1)? as u64,
                area: r.get(2)?,
                parts,
                envelope,
                srs_id: srs,
            })
        })
        .unwrap()
        .map(Result::unwrap)
        .collect();
    (rows, layer, srs_id)
}

fn shoelace(ring: &[(f64, f64)]) -> f64 {
    let mut s = 0.0;
    for w in ring.windows(2) {
        s += w[0].0 * w[1].1 - w[1].0 * w[0].1;
    }
    s / 2.0
}

/// Rasterise the polygons back onto the grid by scanline fill at pixel centres.
/// `None` means no polygon claimed that pixel; a repeat claim panics, because
/// overlapping stands would be a silent disaster in a deliverable.
fn rasterize(
    feats: &[Feature],
    nl: usize,
    ns: usize,
    origin: (f64, f64),
    cell: f64,
) -> Vec<Option<u32>> {
    let mut out = vec![None; nl * ns];
    for f in feats {
        let rings: Vec<&Vec<(f64, f64)>> = f.parts.iter().flatten().collect();
        for row in 0..nl {
            let yc = origin.1 - (row as f64 + 0.5) * cell;
            if yc > f.envelope[3] || yc < f.envelope[2] {
                continue;
            }
            let mut xs: Vec<f64> = Vec::new();
            for r in &rings {
                for w in r.windows(2) {
                    let ((ax, ay), (bx, by)) = (w[0], w[1]);
                    if (ay > yc) != (by > yc) {
                        xs.push(ax + (yc - ay) / (by - ay) * (bx - ax));
                    }
                }
            }
            xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for pair in xs.as_chunks::<2>().0 {
                let c0 = ((pair[0] - origin.0) / cell + 0.5).floor().max(0.0) as usize;
                let c1 = (((pair[1] - origin.0) / cell + 0.5).floor() as usize).min(ns);
                for col in c0..c1 {
                    let p = row * ns + col;
                    assert!(
                        out[p].is_none(),
                        "pixel {row},{col} claimed by regions {:?} and {}",
                        out[p],
                        f.region_id
                    );
                    out[p] = Some(f.region_id);
                }
            }
        }
    }
    out
}

/// The whole claim, in one function: the polygons *are* the region map.
fn assert_polygons_match_raster(gpkg: &Path, envi: &Path, origin: (f64, f64), cell: f64) {
    let (feats, _, layer_srs) = features(gpkg);
    let (raster, nl, ns) = read_envi(envi);

    for f in &feats {
        assert_eq!(f.srs_id, layer_srs, "region {} srs", f.region_id);
        let mut area = 0.0;
        for part in &f.parts {
            for (i, ring) in part.iter().enumerate() {
                assert_eq!(ring[0], *ring.last().unwrap(), "ring not closed");
                let a = shoelace(ring);
                if i == 0 {
                    assert!(a > 0.0, "shell {} must be counter-clockwise", f.region_id);
                } else {
                    assert!(a < 0.0, "hole in {} must be clockwise", f.region_id);
                }
                area += a;
            }
        }
        // Exact for a rectilinear coverage; relative because the coordinates
        // can be in the millions.
        let want = f.n_pixels as f64 * cell * cell;
        assert!(
            (area - want).abs() <= 1e-9 * want.max(f.envelope[1].abs()),
            "region {}: geometry area {area} != {want}",
            f.region_id
        );
        if let Some(a) = f.area {
            assert_eq!(a, want, "region {} area column", f.region_id);
        }
    }

    let painted = rasterize(&feats, nl, ns, origin, cell);
    let mut uncovered = 0usize;
    for (p, want) in raster.iter().enumerate() {
        match painted[p] {
            Some(got) => assert_eq!(
                got,
                *want,
                "pixel {} row {} col {}",
                p,
                p / ns,
                p % ns
            ),
            // Only region 0, the nodata region, may go unclaimed.
            None => {
                assert_eq!(*want, 0, "pixel {p} is region {want} but no polygon covers it");
                uncovered += 1;
            }
        }
    }
    assert_eq!(
        uncovered,
        raster.iter().filter(|&&v| v == 0).count(),
        "exactly the nodata pixels may be left uncovered"
    );
}

#[test]
fn phase_one_and_two_maps_become_polygons_that_are_the_raster() {
    let d = outdir("demo");
    let g = d.to_str().unwrap();
    let mut args = vec!["-o", "demo", "--outdir", g, "--format", "gpkg", DEMO];
    args.splice(0..0, DEMO_ARGS.iter().copied());
    ok(&args);

    // The same run again as ENVI, to compare against.
    let mut args = vec!["-o", "ref", "--outdir", g, DEMO];
    args.splice(0..0, DEMO_ARGS.iter().copied());
    ok(&args);

    // The scene's own tie point and pixel size, from its header.
    let origin = (462405.0, 1741815.0);
    for (gp, rf) in [
        ("demo.rmap.51.gpkg", "ref.rmap.51"),
        ("demo.armap.58.gpkg", "ref.armap.58"),
    ] {
        assert_polygons_match_raster(&d.join(gp), &d.join(rf), origin, 30.0);
    }
}

#[test]
fn eight_way_connectivity_makes_multipart_stands_not_crossed_rings() {
    let d = outdir("eight");
    let g = d.to_str().unwrap();
    let mut args = vec!["-8", "-o", "d8", "--outdir", g, "--format", "gpkg", DEMO];
    args.splice(0..0, DEMO_ARGS.iter().copied());
    ok(&args);
    let mut args = vec!["-8", "-o", "ref", "--outdir", g, DEMO];
    args.splice(0..0, DEMO_ARGS.iter().copied());
    ok(&args);

    // -8 lets a region touch itself diagonally, which is the case that turns
    // into a self-intersecting ring if the junction rule is wrong.
    assert_polygons_match_raster(
        &d.join("d8.armap.72.gpkg"),
        &d.join("ref.armap.72"),
        (462405.0, 1741815.0),
        30.0,
    );
    let (feats, _, _) = features(&d.join("d8.armap.72.gpkg"));
    assert!(
        feats.iter().any(|f| f.parts.len() > 1),
        "8-way segmentation should produce at least one multipart stand"
    );
}

#[test]
fn utm_georeferencing_is_carried_from_the_envi_header() {
    let d = outdir("srs");
    let g = d.to_str().unwrap();
    let mut args = vec!["-o", "s", "--outdir", g, "--format", "gpkg", DEMO];
    args.splice(0..0, DEMO_ARGS.iter().copied());
    ok(&args);

    let conn = Connection::open(d.join("s.armap.58.gpkg")).unwrap();
    let app: i32 = conn.query_row("PRAGMA application_id", [], |r| r.get(0)).unwrap();
    assert_eq!(app, 0x4750_4B47, "GPKG application_id");

    // `map info = {UTM, ..., 15, North, WGS-84}` is EPSG:32615.
    let srs: i32 = conn
        .query_row("SELECT srs_id FROM gpkg_geometry_columns", [], |r| r.get(0))
        .unwrap();
    assert_eq!(srs, 32615);
    let (org, def): (String, String) = conn
        .query_row(
            "SELECT organization, definition FROM gpkg_spatial_ref_sys WHERE srs_id = 32615",
            [],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )
        .unwrap();
    assert_eq!(org, "EPSG");
    assert!(def.contains("UTM zone 15N"), "{def}");
    assert!(def.contains("\"central_meridian\",-93"), "{def}");

    // The extent is the scene's own footprint: 250 pixels of 30 m.
    let (minx, maxx, miny, maxy): (f64, f64, f64, f64) = conn
        .query_row(
            "SELECT min_x, max_x, min_y, max_y FROM gpkg_contents",
            [],
            |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?)),
        )
        .unwrap();
    assert_eq!((minx, maxx), (462405.0, 462405.0 + 7500.0));
    assert_eq!((miny, maxy), (1741815.0 - 7500.0, 1741815.0));
}

#[test]
fn a_projection_with_no_epsg_code_keeps_its_wkt() {
    let d = outdir("wkt");
    let g = d.to_str().unwrap();
    ok(&[
        "--rmap",
        "tests/stage2/e2e_masked/input/rmap",
        "--stage2",
        "tests/stage2/e2e_masked/input/layer",
        "--n2",
        "50,8000",
        "-o",
        "m",
        "--outdir",
        g,
        "--format",
        "gpkg",
    ]);

    let gpkg = d.join("m.armap.39.gpkg");
    let conn = Connection::open(&gpkg).unwrap();
    let (srs, org, def): (i32, String, String) = conn
        .query_row(
            "SELECT g.srs_id, s.organization, s.definition
               FROM gpkg_geometry_columns g
               JOIN gpkg_spatial_ref_sys s USING (srs_id)",
            [],
            |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
        )
        .unwrap();
    // Lambert Conformal Conic on NAD83: no code we can name, but the header's
    // WKT is right there and is preserved rather than dropped.
    assert_eq!(org, "NONE");
    assert!(srs >= 100_000);
    assert!(def.contains("Lambert_Conformal_Conic"), "{def}");
    drop(conn);

    // And the geometry still is the region map, nodata excluded.
    assert_polygons_match_raster(
        &gpkg,
        Path::new("tests/stage2/e2e_masked/expected/armap.39"),
        (-1370910.5240, 753648.1064),
        30.0,
    );
}

#[test]
fn masked_pixels_get_no_polygon() {
    let d = outdir("masked");
    let g = d.to_str().unwrap();
    ok(&[
        "--rmap",
        "tests/stage2/e2e_masked/input/rmap",
        "--stage2",
        "tests/stage2/e2e_masked/input/layer",
        "--n2",
        "50,8000",
        "-o",
        "m",
        "--outdir",
        g,
        "--format",
        "gpkg",
    ]);
    let (feats, _, _) = features(&d.join("m.armap.39.gpkg"));
    assert!(
        feats.iter().all(|f| f.region_id != 0),
        "region 0 is nodata and must not become a stand"
    );

    // The fixture is 44.875% nodata; the polygons must account for the rest and
    // nothing more.
    let (raster, nl, ns) = read_envi(Path::new("tests/stage2/e2e_masked/expected/armap.39"));
    let treed = raster.iter().filter(|&&v| v != 0).count();
    let claimed: u64 = feats.iter().map(|f| f.n_pixels).sum();
    assert_eq!(claimed, treed as u64);
    assert_eq!(nl * ns, 40000);
}

#[test]
fn output_is_byte_identical_when_the_command_is() {
    // `gpkg_contents.last_change` is the only clock-dependent field the format
    // has, and it is pinned to a constant so this holds. It has to be the very
    // same command both times, `--outdir` included: the provenance records the
    // real command line, so a different output directory is a different file
    // on purpose.
    let d = outdir("determinism");
    let g = d.to_str().unwrap().to_string();
    let mut args = vec!["-o", "same", "--outdir", &g, "--format", "gpkg", DEMO];
    args.splice(0..0, DEMO_ARGS.iter().copied());

    ok(&args);
    let first = std::fs::read(d.join("same.armap.58.gpkg")).unwrap();
    ok(&args);
    let second = std::fs::read(d.join("same.armap.58.gpkg")).unwrap();

    assert_eq!(first.len(), second.len(), "same command, different size");
    assert!(first == second, "same command twice must give the same bytes");
}

#[test]
fn an_input_with_no_georeferencing_says_so_and_writes_pixel_coordinates() {
    let d = outdir("nogeo");
    let g = d.to_str().unwrap();
    // The IPW container carries no map info at all.
    let o = ok(&[
        "-t",
        "10",
        "-m",
        ".1",
        "-n",
        "15,15,100,2500,2500",
        "-o",
        "p",
        "--outdir",
        g,
        "--format",
        "gpkg",
        "tests/golden/test_3456/input/test_3456.bip.ipw",
    ]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("no georeferencing"),
        "expected a warning on stderr, got:\n{err}"
    );

    let gpkg = std::fs::read_dir(&d)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .find(|p| {
            p.extension().is_some_and(|e| e == "gpkg")
                && p.file_name().unwrap().to_str().unwrap().contains("armap")
        })
        .expect("a gpkg was written");

    let conn = Connection::open(&gpkg).unwrap();
    let srs: i32 = conn
        .query_row("SELECT srs_id FROM gpkg_geometry_columns", [], |r| r.get(0))
        .unwrap();
    assert_eq!(srs, -1, "undefined cartesian, not a guessed CRS");
    // No transform means no area column, rather than a column of pixel counts
    // wearing the word "area".
    let cols = conn
        .prepare("SELECT * FROM p LIMIT 0")
        .unwrap()
        .column_names()
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    assert!(!cols.contains(&"area".to_string()), "{cols:?}");
    assert!(cols.contains(&"n_pixels".to_string()));
}

/// A GeoTIFF's georeferencing lives in tags, not in a `map info` string, and
/// the reader used to drop it -- harmless while every output was a raster that
/// copied the input's header, wrong the moment polygons had to be placed.
#[test]
fn geotiff_tags_are_read_and_reach_the_geopackage() {
    use tiff::encoder::{colortype, TiffEncoder};
    use tiff::tags::Tag;

    let d = outdir("geotiff");
    let src = d.join("scene.tif");

    // A 16 x 16 two-band scene at 30 m in UTM zone 10N, with a deliberately
    // off-origin tie point so a dropped or mis-signed offset would show.
    const N: u32 = 16;
    let mut px = Vec::new();
    for row in 0..N {
        for col in 0..N {
            // Four quadrants, far enough apart that -t 10 keeps them apart.
            px.push(match (row < N / 2, col < N / 2) {
                (true, true) => 10u8,
                (true, false) => 80,
                (false, true) => 160,
                (false, false) => 240,
            });
        }
    }
    {
        let f = std::fs::File::create(&src).unwrap();
        let mut enc = TiffEncoder::new(f).unwrap();
        let mut img = enc.new_image::<colortype::Gray8>(N, N).unwrap();
        img.encoder()
            .write_tag(Tag::Unknown(33550), &[30.0f64, 30.0, 0.0][..])
            .unwrap();
        img.encoder()
            .write_tag(
                Tag::Unknown(33922),
                &[0.0f64, 0.0, 0.0, 500_000.0, 4_650_000.0, 0.0][..],
            )
            .unwrap();
        // GeoKeyDirectory: version 1.1.0, one key -- 3072 ProjectedCSTypeGeoKey
        // = 32610, held inline (tiffTagLocation 0).
        img.encoder()
            .write_tag(
                Tag::Unknown(34735),
                &[1u16, 1, 0, 1, 3072, 0, 1, 32610][..],
            )
            .unwrap();
        img.write_data(&px).unwrap();
    }

    let g = d.to_str().unwrap();
    ok(&[
        "-t", "10", "-m", "1", "-n", "2,2,4", "-o", "t", "--outdir", g, "--format", "gpkg",
        src.to_str().unwrap(),
    ]);

    let gpkg = std::fs::read_dir(&d)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .find(|p| {
            p.extension().is_some_and(|e| e == "gpkg")
                && p.file_name().unwrap().to_str().unwrap().contains("armap")
        })
        .expect("a gpkg was written");

    let conn = Connection::open(&gpkg).unwrap();
    let srs: i32 = conn
        .query_row("SELECT srs_id FROM gpkg_geometry_columns", [], |r| r.get(0))
        .unwrap();
    assert_eq!(srs, 32610, "ProjectedCSTypeGeoKey should give UTM 10N");

    let (minx, maxx, miny, maxy): (f64, f64, f64, f64) = conn
        .query_row(
            "SELECT min_x, max_x, min_y, max_y FROM gpkg_contents",
            [],
            |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?)),
        )
        .unwrap();
    // The tie point is the upper-left corner; rows run south.
    assert_eq!((minx, maxx), (500_000.0, 500_000.0 + 16.0 * 30.0));
    assert_eq!((miny, maxy), (4_650_000.0 - 16.0 * 30.0, 4_650_000.0));

    // However the scene happens to segment, the polygons must tile it exactly
    // and the areas must follow from the 30 m pixel the tags declare.
    let (px_total, area): (i64, f64) = conn
        .query_row("SELECT sum(n_pixels), sum(area) FROM t", [], |r| {
            Ok((r.get(0)?, r.get(1)?))
        })
        .unwrap();
    assert_eq!(px_total, 256, "every pixel accounted for once");
    assert_eq!(area, 256.0 * 900.0, "30 m pixels, from ModelPixelScaleTag");
}
