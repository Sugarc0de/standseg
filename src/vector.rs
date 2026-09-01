//! Region map to polygons.
//!
//! A region map is the honest output of the segmenter and a poor deliverable:
//! the people who commission a stand map work in polygons, and asking them to
//! polygonise a raster is asking them to install GDAL to finish our job.
//!
//! The method is boundary stitching, not contour tracing. Every side of every
//! pixel whose neighbour is in a different region is emitted as a directed unit
//! edge, oriented so the region always lies to the left of travel; the edges
//! then chain into closed rings. That orientation does all the work downstream:
//! shells and holes come out with opposite winding for free, and no ring can be
//! traced twice.
//!
//! Two properties are worth stating because they are what make the output
//! trustworthy rather than merely plausible:
//!
//! - Vertices land exactly on pixel corners. No smoothing, no simplification,
//!   no Douglas-Peucker. A polygon edge is where the region boundary is, and
//!   two adjacent stands share their vertices exactly, so the coverage has no
//!   slivers and no gaps.
//! - The polygon area of a region is `npix * pixel_area` exactly, holes already
//!   subtracted. That is an identity, not an approximation, which is why the
//!   `area` attribute is computed from the pixel count rather than by shoelace.

use std::collections::HashMap;

use smallvec::SmallVec;

use crate::geo::{self, Transform};

/// A corner of the pixel grid, packed so it can key a hash map in one word.
/// Columns run `0..=nsamps` and rows `0..=nlines`, so an image at this
/// program's 65536-per-axis ceiling still leaves both halves comfortable.
type Corner = u64;

#[inline]
fn corner(col: u32, row: u32) -> Corner {
    ((row as u64) << 32) | col as u64
}

#[inline]
fn unpack(c: Corner) -> (u32, u32) {
    ((c & 0xffff_ffff) as u32, (c >> 32) as u32)
}

/// One boundary side of one pixel, directed so the region is on the left.
#[derive(Clone, Copy)]
struct Edge {
    region: u32,
    from: Corner,
    to: Corner,
}

impl Edge {
    #[inline]
    fn dir(&self) -> (i32, i32) {
        let (fx, fy) = unpack(self.from);
        let (tx, ty) = unpack(self.to);
        (tx as i32 - fx as i32, ty as i32 - fy as i32)
    }
}

/// One region's geometry, in map coordinates.
pub struct RegionPolygons {
    pub id: u32,
    pub npix: u64,
    /// One entry per part. Element 0 of each is the shell, the rest its holes.
    /// Shells wind counter-clockwise and holes clockwise, as OGC asks.
    pub parts: Vec<Vec<Vec<(f64, f64)>>>,
}

impl RegionPolygons {
    /// Bounding box over every ring: `(min_x, min_y, max_x, max_y)`.
    pub fn envelope(&self) -> (f64, f64, f64, f64) {
        let (mut nx, mut ny) = (f64::INFINITY, f64::INFINITY);
        let (mut xx, mut xy) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for ring in self.parts.iter().flatten() {
            for &(x, y) in ring {
                nx = nx.min(x);
                ny = ny.min(y);
                xx = xx.max(x);
                xy = xy.max(y);
            }
        }
        (nx, ny, xx, xy)
    }
}

/// Turn every region in `rband` into polygons.
///
/// `skip` is the id that is not a region -- 0, the artificial region holding
/// masked and nodata pixels, whenever a mask was in play. Passing `None`
/// polygonises everything, which is what an unmasked run wants: there, 0 is a
/// perfectly ordinary stand.
///
/// Output is ordered by region id, so the same map always produces the same
/// file.
pub fn polygonize(
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
    skip: Option<u32>,
    t: &Transform,
) -> Vec<RegionPolygons> {
    let mut edges: Vec<Edge> = Vec::new();
    let mut npix: HashMap<u32, u64> = HashMap::new();

    for row in 0..nlines {
        for col in 0..nsamps {
            let p = row * nsamps + col;
            let r = rband[p];
            if skip == Some(r) {
                continue;
            }
            *npix.entry(r).or_insert(0) += 1;

            let (c, rw) = (col as u32, row as u32);
            // The four corners of this pixel's square, in grid coordinates.
            let (a, b) = (corner(c, rw), corner(c + 1, rw));
            let (cc, d) = (corner(c + 1, rw + 1), corner(c, rw + 1));

            // Each side whose neighbour is a different region, oriented so the
            // interior stays on the left. Chained for one pixel these give
            // b -> a -> d -> cc -> b.
            if row == 0 || rband[p - nsamps] != r {
                edges.push(Edge { region: r, from: b, to: a });
            }
            if col == 0 || rband[p - 1] != r {
                edges.push(Edge { region: r, from: a, to: d });
            }
            if row + 1 == nlines || rband[p + nsamps] != r {
                edges.push(Edge { region: r, from: d, to: cc });
            }
            if col + 1 == nsamps || rband[p + 1] != r {
                edges.push(Edge { region: r, from: cc, to: b });
            }
        }
    }

    // Group by region so each one can be stitched against a small local map
    // rather than one the size of the image.
    edges.sort_unstable_by_key(|e| e.region);

    let mut out = Vec::new();
    let mut i = 0;
    while i < edges.len() {
        let id = edges[i].region;
        let mut j = i;
        while j < edges.len() && edges[j].region == id {
            j += 1;
        }
        let rings = stitch(&edges[i..j]);
        let parts = assemble(rings, t);
        if !parts.is_empty() {
            out.push(RegionPolygons {
                id,
                npix: npix.get(&id).copied().unwrap_or(0),
                parts,
            });
        }
        i = j;
    }
    out
}

/// Chain one region's edges into closed rings, in grid coordinates.
fn stitch(edges: &[Edge]) -> Vec<Vec<Corner>> {
    let mut outgoing: HashMap<Corner, SmallVec<[u32; 4]>> = HashMap::with_capacity(edges.len());
    for (i, e) in edges.iter().enumerate() {
        outgoing.entry(e.from).or_default().push(i as u32);
    }

    let mut used = vec![false; edges.len()];
    let mut rings = Vec::new();

    for start in 0..edges.len() {
        if used[start] {
            continue;
        }
        let origin = edges[start].from;
        let mut ring = vec![origin];
        let mut cur = start;
        loop {
            used[cur] = true;
            let end = edges[cur].to;
            ring.push(end);
            if end == origin {
                break;
            }
            match next_edge(edges, &outgoing, &used, cur, end) {
                Some(n) => cur = n,
                // Unreachable for a well-formed region map: boundary edges of a
                // set of pixels always close. Stopping beats looping if one
                // ever does not.
                None => break,
            }
        }
        // A ring needs three distinct corners plus the repeated close.
        if ring.len() >= 4 && ring[0] == *ring.last().unwrap() {
            rings.push(ring);
        }
    }
    rings
}

/// Pick the next edge at a junction.
///
/// Four boundary edges of the same region meet at a corner where the region
/// touches itself diagonally -- a checkerboard pinch, which `-8` makes
/// reachable. Taking the sharpest turn in the ring's own rotational sense
/// closes off the part we are walking and leaves the other part to be traced
/// as its own ring, which is what makes the result a multipolygon of two
/// squares meeting at a point rather than one self-intersecting ring.
fn next_edge(
    edges: &[Edge],
    outgoing: &HashMap<Corner, SmallVec<[u32; 4]>>,
    used: &[bool],
    cur: usize,
    at: Corner,
) -> Option<usize> {
    let din = edges[cur].dir();
    outgoing
        .get(&at)?
        .iter()
        .map(|&i| i as usize)
        .filter(|&i| !used[i])
        .min_by_key(|&i| turn_score(din, edges[i].dir()))
}

/// Rank a turn: sharpest first, reversal last. Rows run downwards, so a
/// negative cross product is the sense every ring here turns in.
#[inline]
fn turn_score(din: (i32, i32), dout: (i32, i32)) -> u8 {
    if dout == din {
        return 1; // straight on
    }
    match din.0 * dout.1 - din.1 * dout.0 {
        c if c < 0 => 0, // sharpest, hugging the interior
        c if c > 0 => 2,
        _ => 3, // doubling back
    }
}

/// Twice the signed area of a ring in grid coordinates. Shells come out
/// negative and holes positive, which is the sole thing this is used for --
/// map-space orientation is fixed up separately, so a transform that flips an
/// axis cannot turn a stand inside out.
fn signed_area2(ring: &[Corner]) -> f64 {
    let mut s = 0.0;
    for w in ring.windows(2) {
        let (x0, y0) = unpack(w[0]);
        let (x1, y1) = unpack(w[1]);
        s += x0 as f64 * y1 as f64 - x1 as f64 * y0 as f64;
    }
    s
}

/// Sort rings into shells with their holes, and put them in map coordinates.
fn assemble(rings: Vec<Vec<Corner>>, t: &Transform) -> Vec<Vec<Vec<(f64, f64)>>> {
    let (mut shells, mut holes) = (Vec::new(), Vec::new());
    for r in rings {
        if signed_area2(&r) < 0.0 {
            shells.push(r);
        } else {
            holes.push(r);
        }
    }
    if shells.is_empty() {
        return Vec::new();
    }

    let mut parts: Vec<Vec<Vec<Corner>>> = shells.into_iter().map(|s| vec![s]).collect();
    for h in holes {
        let owner = if parts.len() == 1 {
            // The overwhelmingly common case, and it needs no geometry.
            0
        } else {
            let (px, py) = inside_hole(&h);
            parts
                .iter()
                .enumerate()
                .filter(|(_, p)| point_in_ring(&p[0], px, py))
                // Nested shells of one region are possible: an annulus with an
                // island of itself inside. The innermost containing shell owns
                // the hole.
                .min_by(|a, b| {
                    let (sa, sb) = (signed_area2(&a.1[0]).abs(), signed_area2(&b.1[0]).abs());
                    sa.partial_cmp(&sb).unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0)
        };
        parts[owner].push(h);
    }

    parts
        .into_iter()
        .map(|part| {
            part.into_iter()
                .enumerate()
                .map(|(i, ring)| to_map(&ring, t, i == 0))
                .collect()
        })
        .collect()
}

/// A point strictly inside a hole's cavity, in grid coordinates.
///
/// The topmost horizontal edge of a cavity always has cavity pixels directly
/// below it, so that pixel's centre is inside and -- being on a half-integer --
/// can never lie on a ring, which keeps the containment test exact.
fn inside_hole(ring: &[Corner]) -> (f64, f64) {
    let mut best: Option<(u32, u32)> = None;
    for w in ring.windows(2) {
        let (x0, y0) = unpack(w[0]);
        let (x1, y1) = unpack(w[1]);
        if y0 != y1 {
            continue; // vertical
        }
        let x = x0.min(x1);
        if best.is_none_or(|(_, by)| y0 < by) {
            best = Some((x, y0));
        }
    }
    let (x, y) = best.unwrap_or_else(|| unpack(ring[0]));
    (x as f64 + 0.5, y as f64 + 0.5)
}

/// Ray casting. The rings here are rectilinear and the test point is always on
/// a half-integer, so no vertex or edge case can be hit.
fn point_in_ring(ring: &[Corner], px: f64, py: f64) -> bool {
    let mut inside = false;
    for w in ring.windows(2) {
        let (x0, y0) = unpack(w[0]);
        let (x1, y1) = unpack(w[1]);
        let (x0, y0) = (x0 as f64, y0 as f64);
        let (x1, y1) = (x1 as f64, y1 as f64);
        if (y0 > py) != (y1 > py) {
            let x = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
            if px < x {
                inside = !inside;
            }
        }
    }
    inside
}

/// Grid ring to map coordinates, wound the way OGC wants: shells
/// counter-clockwise, holes clockwise, measured in the map's own frame.
fn to_map(ring: &[Corner], t: &Transform, is_shell: bool) -> Vec<(f64, f64)> {
    let mut pts: Vec<(f64, f64)> = ring
        .iter()
        .map(|&c| {
            let (x, y) = unpack(c);
            geo::apply(t, x as f64, y as f64)
        })
        .collect();

    let mut s = 0.0;
    for w in pts.windows(2) {
        s += w[0].0 * w[1].1 - w[1].0 * w[0].1;
    }
    if (s > 0.0) != is_shell {
        pts.reverse();
    }
    pts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo::PIXEL_SPACE;

    /// Shoelace area of a whole part, holes subtracted.
    fn part_area(part: &[Vec<(f64, f64)>]) -> f64 {
        part.iter()
            .map(|ring| {
                let mut s = 0.0;
                for w in ring.windows(2) {
                    s += w[0].0 * w[1].1 - w[1].0 * w[0].1;
                }
                s / 2.0
            })
            .sum()
    }

    fn go(rband: &[u32], nl: usize, ns: usize) -> Vec<RegionPolygons> {
        polygonize(rband, nl, ns, None, &PIXEL_SPACE)
    }

    #[test]
    fn single_pixel_is_a_unit_square() {
        let p = go(&[7], 1, 1);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].id, 7);
        assert_eq!(p[0].npix, 1);
        assert_eq!(p[0].parts.len(), 1);
        assert_eq!(p[0].parts[0].len(), 1, "no holes");
        // Closed, four corners plus the repeat.
        assert_eq!(p[0].parts[0][0].len(), 5);
        assert_eq!(part_area(&p[0].parts[0]), 1.0, "counter-clockwise unit square");
    }

    #[test]
    fn two_regions_share_their_boundary_exactly() {
        // 2x2, split down the middle.
        let p = go(&[1, 2, 1, 2], 2, 2);
        assert_eq!(p.len(), 2);
        for r in &p {
            assert_eq!(r.npix, 2);
            assert_eq!(part_area(&r.parts[0]), 2.0);
        }
        // The shared edge is the segment x=1 from y=0 to y=2, and it appears in
        // both rings, so the two polygons abut with no sliver.
        let has = |r: &RegionPolygons, pt: (f64, f64)| r.parts[0][0].contains(&pt);
        for r in &p {
            assert!(has(r, (1.0, 0.0)) && has(r, (1.0, 2.0)));
        }
    }

    #[test]
    fn a_ring_region_gets_a_hole() {
        // 3x3 of region 1 with a single pixel of region 2 at the centre.
        let p = go(&[1, 1, 1, 1, 2, 1, 1, 1, 1], 3, 3);
        assert_eq!(p.len(), 2);
        let one = p.iter().find(|r| r.id == 1).unwrap();
        assert_eq!(one.npix, 8);
        assert_eq!(one.parts.len(), 1, "one shell");
        assert_eq!(one.parts[0].len(), 2, "shell plus one hole");
        // Holes wind the other way, so the shoelace sum is already net of it.
        assert_eq!(part_area(&one.parts[0]), 8.0);
    }

    #[test]
    fn nested_island_inside_a_hole() {
        // Region 1 as a 5x5 annulus, region 2 as the ring inside it, and
        // region 1 again as the single centre pixel: 1's hole contains 1.
        #[rustfmt::skip]
        let m = [
            1, 1, 1, 1, 1,
            1, 2, 2, 2, 1,
            1, 2, 1, 2, 1,
            1, 2, 2, 2, 1,
            1, 1, 1, 1, 1,
        ];
        let p = go(&m, 5, 5);
        let one = p.iter().find(|r| r.id == 1).unwrap();
        assert_eq!(one.npix, 17);
        assert_eq!(one.parts.len(), 2, "outer annulus and the centre island");
        let total: f64 = one.parts.iter().map(|x| part_area(x)).sum();
        assert_eq!(total, 17.0, "area is the pixel count, holes subtracted");
        // The hole belongs to the annulus, not to the island.
        let big = one
            .parts
            .iter()
            .max_by(|a, b| part_area(a).partial_cmp(&part_area(b)).unwrap())
            .unwrap();
        assert_eq!(big.len(), 2);
    }

    #[test]
    fn diagonal_touch_is_two_parts_not_a_crossed_ring() {
        // The checkerboard pinch: region 1 meets itself at exactly one corner.
        let p = go(&[1, 2, 2, 1], 2, 2);
        let one = p.iter().find(|r| r.id == 1).unwrap();
        assert_eq!(one.npix, 2);
        assert_eq!(one.parts.len(), 2, "two squares, not one figure-of-eight");
        for part in &one.parts {
            assert_eq!(part.len(), 1);
            assert_eq!(part_area(part), 1.0);
        }
    }

    #[test]
    fn area_is_the_pixel_count_for_every_region() {
        // A map with enough awkwardness to be worth trusting: nested rings,
        // diagonal touches and single pixels.
        #[rustfmt::skip]
        let m = [
            1, 1, 1, 1, 2, 2,
            1, 3, 3, 1, 2, 4,
            1, 3, 5, 1, 2, 2,
            1, 1, 1, 1, 2, 2,
            6, 6, 2, 2, 2, 2,
            6, 2, 2, 7, 7, 2,
        ];
        for r in go(&m, 6, 6) {
            let a: f64 = r.parts.iter().map(|p| part_area(p)).sum();
            assert_eq!(a, r.npix as f64, "region {}", r.id);
        }
    }

    #[test]
    fn skip_omits_the_nodata_region_entirely() {
        let m = [0, 1, 0, 1];
        let p = polygonize(&m, 2, 2, Some(0), &PIXEL_SPACE);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].id, 1);
        assert_eq!(p[0].npix, 2);
        // Region 1 owns the right-hand column, so it is one 1x2 part -- and
        // crucially the nodata column beside it is absent rather than being a
        // polygon of its own.
        assert_eq!(p[0].parts.len(), 1);
        assert_eq!(part_area(&p[0].parts[0]), 2.0);
    }

    #[test]
    fn a_transform_places_and_scales_the_polygons() {
        // North-up, 30 m pixels, tie point at a UTM easting/northing.
        let t = [500000.0, 30.0, 0.0, 4650000.0, 0.0, -30.0];
        let p = polygonize(&[9], 1, 1, None, &t);
        let (nx, ny, xx, xy) = p[0].envelope();
        assert_eq!((nx, xx), (500000.0, 500030.0));
        assert_eq!((ny, xy), (4649970.0, 4650000.0));
        // Rows run south, but the ring must still wind counter-clockwise on
        // the ground.
        assert_eq!(part_area(&p[0].parts[0]), 900.0);
        assert_eq!(geo::pixel_area(&t) * p[0].npix as f64, 900.0);
    }

    #[test]
    fn output_is_ordered_by_region_id() {
        let p = go(&[5, 3, 9, 1], 2, 2);
        let ids: Vec<u32> = p.iter().map(|r| r.id).collect();
        assert_eq!(ids, vec![1, 3, 5, 9]);
    }
}
