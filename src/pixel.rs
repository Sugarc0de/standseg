//! Phase 0: the two sweeps over the raw image that happen before any region
//! exists, plus the construction of the initial region list.
//!
//! Faithful port of `pixel.c`. The purpose of the whole phase is to shrink the
//! initial region list -- on Case 1 it turns 62500 pixels into 55226 regions.

use crate::config::SegConfig;
use crate::contig::{CCLEAR, CMONO, E_EDGE, N_EDGE, S_EDGE, W_EDGE};
use crate::image::{Image, IntSample, Raster, RasterRef};
use crate::region::{merge_regions, RegionId, RegionList};

/// The C's `MAXLONG` sentinel, used when a neighbour is out of bounds or masked.
const MAXLONG: i64 = i64::MAX;

pub struct Bands {
    pub cband: Vec<u8>,
    pub rband: Vec<RegionId>,
    pub nreg: usize,
}

/// Squared distance between two pixel vectors, in integer arithmetic.
///
/// Note the asymmetry with `RegionList::dist2`: Phase 0 compares *pixels*
/// exactly in integers, everything later compares *centroids* in f32.
///
/// i64 is wide enough for every sample type we accept: the worst case is 16-bit
/// samples over the 255-band maximum, 255 * 65535^2 = 1.1e15, well inside i64.
#[inline]
fn pix_dist2<T: IntSample>(a: &[T], b: &[T]) -> i64 {
    let mut dist2: i64 = 0;
    for i in 0..a.len() {
        let diff = a[i].to_i64() - b[i].to_i64();
        dist2 += diff * diff;
    }
    dist2
}

/// For every pixel, flag every neighbour sitting at its minimum distance --
/// but only if that minimum is within the general tolerance.
fn pix_nnbr<T: IntSample>(
    img: &Raster<'_, T>,
    cfg: &SegConfig,
    mask: Option<&[u8]>,
    cband: &mut [u8],
    tg2: f32,
) {
    let (nl, ns) = (img.nlines, img.nsamps);
    let conn = &cfg.conn;
    let mut ndist2 = [MAXLONG; 8];

    for l in 0..nl {
        for s in 0..ns {
            let p = l * ns + s;
            cband[p] = CCLEAR;
            if let Some(m) = mask {
                if m[p] == 0 {
                    continue;
                }
            }
            let cpix = img.pixel(l, s);
            let mut mdist2 = MAXLONG;

            for (d, &(dx, dy)) in conn.deltas.iter().enumerate().take(conn.ncdir) {
                let (nx, ny) = (s as i32 + dx, l as i32 + dy);
                let in_bounds = nx >= 0 && ny >= 0 && (nx as usize) < ns && (ny as usize) < nl;
                let ok = in_bounds && mask.is_none_or(|m| m[ny as usize * ns + nx as usize] != 0);
                ndist2[d] = if ok {
                    pix_dist2(cpix, img.pixel(ny as usize, nx as usize))
                } else {
                    MAXLONG
                };
                mdist2 = mdist2.min(ndist2[d]);
            }

            // The C is `long <= float`, which promotes the long to float. We
            // do the comparison in f64 instead. For 8-bit input that is not a
            // change: mdist2 tops out at 255^2 * nbands, which is exact in f32
            // as well as f64, and tg2 promotes exactly -- so the boolean is
            // identical, and the golden fixtures still hold. For 16-bit input
            // f32 would start rounding, so f64 is the honest width.
            if (mdist2 as f64) <= tg2 as f64 {
                for (d, &dist) in ndist2.iter().enumerate().take(conn.ncdir) {
                    if dist == mdist2 {
                        conn.set(&mut cband[p], d);
                    }
                }
            }
        }
    }
}

/// Pair up mutually-nearest pixels, at most one merge per pixel, assigning
/// region ids as we go.
fn pix_merge(nl: usize, ns: usize, cfg: &SegConfig, mask: Option<&[u8]>, bands: &mut Bands) {
    let conn = &cfg.conn;
    let mut nregions: RegionId = 0;
    // A rotating start direction, carried across pixels. Easy to miss and it
    // changes which pairs form.
    let mut idir: usize = 0;

    for l in 0..nl {
        for s in 0..ns {
            let p = l * ns + s;
            if mask.is_some_and(|m| m[p] == 0) || bands.rband[p] > 0 {
                continue;
            }
            if bands.cband[p] == CCLEAR {
                nregions += 1;
                bands.rband[p] = nregions;
                continue;
            }

            let mut d = idir;
            let mut merged = false;
            // The C is `while (advance_dir(d) != idir)`, which pre-advances and
            // therefore never tries direction `idir` itself.
            loop {
                d += 1;
                if d == conn.ncdir {
                    d = 0;
                }
                if d == idir {
                    break;
                }
                if !conn.has(bands.cband[p], d) {
                    continue;
                }
                let (dx, dy) = conn.deltas[d];
                let np = (l as i32 + dy) as usize * ns + (s as i32 + dx) as usize;
                if bands.rband[np] > 0 {
                    continue;
                }
                if conn.has(bands.cband[np], conn.reverse(d)) {
                    nregions += 1;
                    bands.rband[p] = nregions;
                    bands.rband[np] = nregions;
                    bands.cband[p] = conn.flags[d];
                    bands.cband[np] = conn.flags[conn.reverse(d)];
                    idir = d;
                    merged = true;
                    break;
                }
            }

            if !merged {
                bands.cband[p] = CMONO;
                nregions += 1;
                bands.rband[p] = nregions;
            }
        }
    }

    bands.nreg = nregions as usize;
}

/// Mark directions that are out of bounds or masked as "contiguous", so the
/// region routines never look outside the image or into nodata.
///
/// This is what makes the unchecked neighbour reads in `merge_regions` and
/// `reg_nnbr` safe, and it is also what stops regions growing along a shoreline.
fn pix_check_bounds_and_mask(
    cfg: &SegConfig,
    mask: Option<&[u8]>,
    cband: &mut [u8],
    nl: usize,
    ns: usize,
    x: usize,
    y: usize,
) {
    let p = y * ns + x;
    if y == 0 {
        cband[p] |= N_EDGE;
    }
    if y == nl - 1 {
        cband[p] |= S_EDGE;
    }
    if x == 0 {
        cband[p] |= W_EDGE;
    }
    if x == ns - 1 {
        cband[p] |= E_EDGE;
    }

    if let Some(m) = mask {
        let conn = &cfg.conn;
        for d in 0..conn.ncdir {
            if !conn.has(cband[p], d) {
                let (dx, dy) = conn.deltas[d];
                let np = (y as i32 + dy) as usize * ns + (x as i32 + dx) as usize;
                if m[np] == 0 {
                    conn.set(&mut cband[p], d);
                }
            }
        }
    }
}

/// Turn the pixel pairings into the initial region list.
fn make_region_list<T: IntSample>(
    img: &Raster<'_, T>,
    cfg: &SegConfig,
    mask: Option<&[u8]>,
    bands: &mut Bands,
    rl: &mut RegionList,
) -> Result<(), String> {
    let (nl, ns) = (img.nlines, img.nsamps);
    let conn = &cfg.conn;
    let offs = conn.offsets(ns);
    let dummy: RegionId = bands.nreg as RegionId + 1;

    for l in 0..nl {
        for s in 0..ns {
            let p = l * ns + s;
            let rid = bands.rband[p];
            if rid == 0 {
                continue;
            }
            if rl.is_active(rid) {
                continue;
            }
            let cdf = bands.cband[p];
            if cdf != 0 {
                // Exactly one bit is set here (pix_merge guarantees it).
                let d = conn.flags.iter().position(|&f| f == cdf).ok_or_else(|| {
                    format!("contiguity byte {cdf:#x} at ({s},{l}) is not a single direction flag")
                })?;
                let (dx, dy) = conn.deltas[d];
                let (nx, ny) = ((s as i32 + dx) as usize, (l as i32 + dy) as usize);

                rl.from_pixel(rid, s as u16, l as u16, img.pixel(l, s));
                rl.from_pixel(dummy, nx as u16, ny as u16, img.pixel(ny, nx));
                pix_check_bounds_and_mask(cfg, mask, &mut bands.cband, nl, ns, s, l);
                pix_check_bounds_and_mask(cfg, mask, &mut bands.cband, nl, ns, nx, ny);
                merge_regions(
                    rl,
                    &mut bands.rband,
                    &mut bands.cband,
                    ns,
                    conn,
                    &offs,
                    rid,
                    dummy,
                )?;
            } else {
                rl.from_pixel(rid, s as u16, l as u16, img.pixel(l, s));
                pix_check_bounds_and_mask(cfg, mask, &mut bands.cband, nl, ns, s, l);
            }
        }
    }
    Ok(())
}

/// Run all of Phase 0. Returns the bands and the initial region list.
///
/// Dispatches once on the input sample width; everything below is monomorphic,
/// so the 8-bit path compiles to exactly what it did before 16-bit existed.
pub fn phase0(
    img: &Image,
    cfg: &SegConfig,
    mask: Option<&[u8]>,
) -> Result<(Bands, RegionList), String> {
    img.with_raster(|r| match r {
        RasterRef::U8(r) => phase0_typed(&r, cfg, mask),
        RasterRef::U16(r) => phase0_typed(&r, cfg, mask),
        RasterRef::I16(r) => phase0_typed(&r, cfg, mask),
        // Unreachable in practice -- main.rs refuses a float input with a fuller
        // message before we get here -- but stated rather than left to a panic.
        RasterRef::F32(_) => Err(
            "the first stage segments 8- and 16-bit integer imagery only; \
             32-bit float is accepted as the --stage2 image, not as the input"
                .to_string(),
        ),
    })
}

fn phase0_typed<T: IntSample>(
    img: &Raster<'_, T>,
    cfg: &SegConfig,
    mask: Option<&[u8]>,
) -> Result<(Bands, RegionList), String> {
    let npix = img.npixels();
    let tg2 = cfg.tols[0] * cfg.tols[0];

    let mut bands = Bands {
        cband: vec![0u8; npix],
        rband: vec![0 as RegionId; npix],
        nreg: 0,
    };

    pix_nnbr(img, cfg, mask, &mut bands.cband, tg2);
    pix_merge(img.nlines, img.nsamps, cfg, mask, &mut bands);

    // Two extra slots: id 0 for masked/nodata pixels, and the scratch region.
    let mut rl = RegionList::new(bands.nreg + 2, img.nbands);
    make_region_list(img, cfg, mask, &mut bands, &mut rl)?;

    Ok((bands, rl))
}
