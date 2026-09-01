//! Where the pixels sit on the ground.
//!
//! The raster outputs never had to answer this: a region map is the same shape
//! as its input, so carrying the georeferencing string through untouched was
//! enough. A polygon has to be *placed*, so both containers' georeferencing has
//! to be reduced to the same six numbers here.
//!
//! ENVI carries it as a `map info` string, GeoTIFF as a pair of tag arrays
//! (`src/io/tiff.rs` reads those). Neither is authoritative over the other; a
//! file that has neither produces polygons in pixel coordinates, which is a
//! usable answer and is flagged as such rather than guessed at.

/// Pixel-to-map affine, in GDAL's parameter order:
///
/// ```text
/// x = t[0] + col * t[1] + row * t[2]
/// y = t[3] + col * t[4] + row * t[5]
/// ```
///
/// `(col, row) = (0, 0)` is the upper-left **corner** of the upper-left pixel,
/// not its centre. Polygon vertices are corners, so this is the convention that
/// costs no half-pixel correction downstream -- and it is the one both ENVI and
/// GeoTIFF use for their tie points, so nothing has to be shifted on the way in
/// either.
pub type Transform = [f64; 6];

/// Pixel coordinates: one map unit per pixel, y increasing downwards. What a
/// file with no georeferencing gets.
pub const PIXEL_SPACE: Transform = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Map coordinate of a pixel *corner*.
#[inline]
pub fn apply(t: &Transform, col: f64, row: f64) -> (f64, f64) {
    (
        t[0] + col * t[1] + row * t[2],
        t[3] + col * t[4] + row * t[5],
    )
}

/// Area of one pixel in square map units, from the determinant of the linear
/// part. Signed area would carry the north-up flip; nobody wants a negative
/// stand, so take the magnitude.
pub fn pixel_area(t: &Transform) -> f64 {
    (t[1] * t[5] - t[2] * t[4]).abs()
}

/// Parse an ENVI `map info` string into a transform and, where it can be
/// worked out, an EPSG code.
///
/// The string is the comma-separated body of the header record, without the
/// surrounding braces:
///
/// ```text
/// UTM, 1, 1, 500000.0, 4649776.0, 30.0, 30.0, 10, North, WGS-84, units=Meters
/// ```
///
/// Fields 1 and 2 are the reference pixel, **1-based**, naming the upper-left
/// corner of that pixel; 3 and 4 are its easting and northing; 5 and 6 the
/// pixel size. Everything after that is projection identification and varies by
/// projection, so it is read defensively -- a `map info` we only half
/// understand still yields a usable transform, which is the part that matters.
pub fn from_envi_map_info(s: &str) -> Option<(Transform, Option<u32>)> {
    let f: Vec<&str> = s.split(',').map(str::trim).collect();
    if f.len() < 7 {
        return None;
    }
    let refx: f64 = f[1].parse().ok()?;
    let refy: f64 = f[2].parse().ok()?;
    let east: f64 = f[3].parse().ok()?;
    let north: f64 = f[4].parse().ok()?;
    let xsize: f64 = f[5].parse().ok()?;
    let ysize: f64 = f[6].parse().ok()?;
    if xsize == 0.0 || ysize == 0.0 || !xsize.is_finite() || !ysize.is_finite() {
        return None;
    }
    // Back the tie point out to the image's own upper-left corner. ENVI's y
    // pixel size is written positive and the rows run south, hence the sign
    // flip in the last term.
    let x0 = east - (refx - 1.0) * xsize;
    let y0 = north + (refy - 1.0) * ysize;
    Some(([x0, xsize, 0.0, y0, 0.0, -ysize], envi_epsg(&f)))
}

/// EPSG code from the projection fields of a `map info`, for the two cases that
/// cover essentially all forestry imagery. Anything else returns `None` and the
/// GeoPackage falls back to the header's `coordinate system string`.
fn envi_epsg(f: &[&str]) -> Option<u32> {
    fn wgs84(d: &str) -> bool {
        let d = d.to_ascii_uppercase();
        d.contains("WGS") && d.contains("84")
    }

    let proj = f[0];
    if proj.eq_ignore_ascii_case("UTM") {
        // UTM, refx, refy, east, north, dx, dy, zone, hemisphere, datum, ...
        let zone: u32 = f.get(7)?.trim().parse().ok()?;
        if !(1..=60).contains(&zone) {
            return None;
        }
        if !wgs84(f.get(9).copied().unwrap_or("")) {
            return None;
        }
        let south = f.get(8)?.trim().to_ascii_uppercase().starts_with('S');
        return Some(if south { 32700 + zone } else { 32600 + zone });
    }
    if proj.to_ascii_lowercase().contains("geographic") && wgs84(f.get(7).copied().unwrap_or("")) {
        return Some(4326);
    }
    None
}

/// A WKT1 definition for the EPSG codes we can generate one for.
///
/// GeoPackage wants a `definition` in `gpkg_spatial_ref_sys`. Readers in
/// practice resolve the CRS from `organization` plus
/// `organization_coordsys_id`, so a code we cannot spell out is still read
/// correctly -- but WGS 84 and its UTM zones are formulaic, they are what
/// Landsat Collection 2 ships in, and writing them out means the common case
/// needs no EPSG database at the other end.
pub fn wkt_for_epsg(epsg: u32) -> Option<String> {
    const GEOGCS: &str = "GEOGCS[\"WGS 84\",\
DATUM[\"WGS_1984\",\
SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],\
AUTHORITY[\"EPSG\",\"6326\"]],\
PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],\
UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],\
AUTHORITY[\"EPSG\",\"4326\"]]";

    if epsg == 4326 {
        return Some(GEOGCS.to_string());
    }
    let (zone, south) = match epsg {
        32601..=32660 => (epsg - 32600, false),
        32701..=32760 => (epsg - 32700, true),
        _ => return None,
    };
    let meridian = zone as i32 * 6 - 183;
    let false_northing = if south { 10000000 } else { 0 };
    let hemi = if south { 'S' } else { 'N' };
    Some(format!(
        "PROJCS[\"WGS 84 / UTM zone {zone}{hemi}\",{GEOGCS},\
PROJECTION[\"Transverse_Mercator\"],\
PARAMETER[\"latitude_of_origin\",0],\
PARAMETER[\"central_meridian\",{meridian}],\
PARAMETER[\"scale_factor\",0.9996],\
PARAMETER[\"false_easting\",500000],\
PARAMETER[\"false_northing\",{false_northing}],\
UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],\
AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],\
AUTHORITY[\"EPSG\",\"{epsg}\"]]"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn envi_utm_north() {
        let (t, epsg) = from_envi_map_info(
            "UTM, 1, 1, 500000.0, 4649776.0, 30.0, 30.0, 10, North, WGS-84, units=Meters",
        )
        .unwrap();
        assert_eq!(epsg, Some(32610));
        assert_eq!(t, [500000.0, 30.0, 0.0, 4649776.0, 0.0, -30.0]);
        // Corner (0,0) is the tie point; one pixel across and down is +30 east
        // and 30 *south*.
        assert_eq!(apply(&t, 0.0, 0.0), (500000.0, 4649776.0));
        assert_eq!(apply(&t, 1.0, 1.0), (500030.0, 4649746.0));
        assert_eq!(pixel_area(&t), 900.0);
    }

    #[test]
    fn envi_reference_pixel_is_one_based() {
        // A tie point on pixel (1,1) -- ENVI's own origin -- and on pixel (2,3)
        // describing the same grid must give the same transform.
        let a = from_envi_map_info("UTM, 1, 1, 500000, 4649776, 30, 30, 10, North, WGS-84")
            .unwrap()
            .0;
        let b = from_envi_map_info("UTM, 2, 3, 500030, 4649716, 30, 30, 10, North, WGS-84")
            .unwrap()
            .0;
        assert_eq!(a, b);
    }

    #[test]
    fn envi_utm_south() {
        let (_, epsg) =
            from_envi_map_info("UTM, 1, 1, 500000, 6000000, 30, 30, 55, South, WGS-84").unwrap();
        assert_eq!(epsg, Some(32755));
    }

    #[test]
    fn envi_geographic() {
        let (t, epsg) =
            from_envi_map_info("Geographic Lat/Lon, 1, 1, -123.0, 45.0, 0.001, 0.001, WGS-84")
                .unwrap();
        assert_eq!(epsg, Some(4326));
        assert_eq!(apply(&t, 0.0, 0.0), (-123.0, 45.0));
    }

    #[test]
    fn envi_unknown_projection_still_gives_a_transform() {
        let (t, epsg) =
            from_envi_map_info("Albers Conical Equal Area, 1, 1, 100, 200, 30, 30, Custom").unwrap();
        assert_eq!(epsg, None);
        assert_eq!(t[1], 30.0);
    }

    #[test]
    fn envi_rubbish_is_rejected_rather_than_guessed() {
        assert!(from_envi_map_info("UTM, 1, 1").is_none());
        assert!(from_envi_map_info("UTM, 1, 1, x, y, 30, 30, 10, North, WGS-84").is_none());
        // A zero pixel size would make every polygon degenerate.
        assert!(from_envi_map_info("UTM, 1, 1, 0, 0, 0, 0, 10, North, WGS-84").is_none());
    }

    #[test]
    fn utm_wkt_central_meridian() {
        let w = wkt_for_epsg(32610).unwrap();
        assert!(w.contains("\"central_meridian\",-123"), "{w}");
        assert!(w.contains("UTM zone 10N"));
        assert!(wkt_for_epsg(32755).unwrap().contains("false_northing\",10000000"));
        assert!(wkt_for_epsg(3857).is_none());
    }
}
