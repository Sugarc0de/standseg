//! Segment development — the second phase of Ye et al. (2025).
//!
//! Takes the region map stage 1 produced and merges those micro-segments using
//! a **second image over the same grid**: forest structure, age, or species
//! rather than the spectral proxies stage 1 saw. See PLAN.md section 13 and
//! `tests/STAGE2.md`; the oracle is `tools/stage2_oracle/`.
//!
//! Deliberately not part of [`crate::segment::Segmenter`]. The rules differ at
//! almost every point — f64 centroids instead of f32, the *absorbing* region's
//! id surviving instead of the lower one, a forced-mutual write-back instead of
//! a mutual-nearest test, no contiguity band, no RNG — so folding it in would
//! mean threading a mode flag through the hot loop of a phase that is currently
//! byte-exact against 1992.
//!
//! ## Float layers and numpy's summation order
//!
//! The oracle's initial centroids are `b_images.mean(axis=1)`, and numpy
//! accumulates a float32 mean *in float32*. That is not a rounding detail to be
//! improved on: summing the same pixels in f64 instead moves 5.7 % of the output
//! pixels on a 1000^2 crop of `exp_150`. Reproducing the oracle means
//! reproducing its arithmetic, summation order included.
//!
//! Which order that is depends on something the oracle never intended to
//! choose. `b_images = image[:, coords[:, 0], coords[:, 1]]` is
//! **non-contiguous when there is more than one band** and contiguous when there
//! is exactly one, because the advanced-index result is laid out band-last. numpy
//! uses its pairwise summation only on contiguous input; on strided input it
//! accumulates sequentially. So:
//!
//! | bands | `b_images` | numpy sums |
//! |---|---|---|
//! | 1 | contiguous | pairwise (`npy_pairwise_sum`, `PW_BLOCKSIZE` = 128) |
//! | 2+ | strided | sequentially, in `coords` order |
//!
//! Verified against numpy 2.4.3 over 285 cases spanning 1-6 bands and regions of
//! 1 to 12 000 pixels: zero mismatches in either the sums or the means. Both
//! orders are implemented below and picked by band count, which is why a
//! single-band float layer costs a gather buffer and a multi-band one does not.

use std::collections::HashMap;

use crate::image::{Image, Raster, RasterRef, Sample};

#[derive(Debug, Clone)]
pub struct Stage2Config {
    /// Minimum region size. A region below this looks for a partner.
    pub nmin: u32,
    /// Absolute maximum. A merge whose *sum* would exceed this is refused.
    pub nmax: u32,
}

/// Per-pass counters, the stage-2 equivalent of `myseg.log`: enough to localise
/// a divergence to one pass and one rejection reason.
#[derive(Debug, Clone, Copy, Default)]
pub struct PassStats2 {
    pub considered: usize,
    pub busy: usize,
    pub no_cand: usize,
    pub inf: usize,
    pub over_max: usize,
    pub not_mutual: usize,
    pub merged: usize,
    pub nreg: usize,
}

#[derive(Debug, Clone)]
pub struct Stage2Result {
    /// Counted the oracle's way: passes executed, including the final one that
    /// merged nothing, plus one. This is the number in the output filename.
    pub passes: usize,
    pub nreg: usize,
    pub dropped_majority_nodata: usize,
    pub stats: Vec<PassStats2>,
}

/// Python's `math.isclose`, including its infinity rule.
///
/// Getting this wrong is a silent divergence: a naive
/// `|a-b| <= rel*max(|a|,|b|)` returns *true* when one side is infinite, and the
/// running best distance starts at infinity, so every first candidate would
/// register as a tie and be kept rather than taken.
#[inline]
fn isclose(a: f64, b: f64, rel_tol: f64) -> bool {
    if a == b {
        return true;
    }
    if a.is_infinite() || b.is_infinite() {
        return false;
    }
    let diff = (b - a).abs();
    diff <= (rel_tol * b).abs() || diff <= (rel_tol * a).abs()
}

/// Tolerance on the near-tie test when choosing among candidate neighbours.
const TIE_REL_TOL: f64 = 1e-6;
/// Tolerance on the "are these two still each other's partner" test.
const MUTUAL_REL_TOL: f64 = 1e-9;

struct Reg {
    id: u32,
    npix: u32,
    ulx: u32,
    uly: u32,
    lrx: u32,
    lry: u32,
    alive: bool,
    nnbr_id: u32,
    nnbr_d2: f64,
    /// 0 = free this pass, 1 = has absorbed, 2 = has been absorbed.
    state: u8,
}

struct Regions {
    regs: Vec<Reg>,
    /// `nbands` f64 per region, in the same order as `regs`.
    ctr: Vec<f64>,
    nbands: usize,
    idx: HashMap<u32, u32>,
}

impl Regions {
    #[inline]
    fn dist2(&self, i: usize, j: usize) -> f64 {
        let (oi, oj) = (i * self.nbands, j * self.nbands);
        let mut d = 0.0f64;
        for b in 0..self.nbands {
            let t = self.ctr[oi + b] - self.ctr[oj + b];
            d += t * t;
        }
        d
    }
}

/// Build the region list from the map, and take the stage-2 centroids.
///
/// For integer samples the means are exactly reproducible against numpy without
/// imitating its pairwise summation: every partial sum is an exact integer below
/// 2^53, so summation order cannot matter, and one `i64` sum with one `f64`
/// divide gives the same bits. For `f32` it very much can matter, so that path
/// reproduces numpy's order instead -- see `pairwise_sum_f32`.
//
/// numpy's pairwise summation for contiguous `float32` -- the single-band case.
///
/// `npy_pairwise_sum` (numpy `loops_utils.h`): naive below 8, eight interleaved
/// partial sums combined in a fixed tree up to `PW_BLOCKSIZE` = 128, and
/// split-in-half above it with the split rounded down to a multiple of 8.
fn pairwise_sum_f32(a: &[f32]) -> f32 {
    let n = a.len();
    if n < 8 {
        let mut r = 0.0f32;
        for &v in a {
            r += v;
        }
        return r;
    }
    if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - (n % 8) {
            for k in 0..8 {
                r[k] += a[i + k];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        return res;
    }
    let mut n2 = n / 2;
    n2 -= n2 % 8;
    pairwise_sum_f32(&a[..n2]) + pairwise_sum_f32(&a[n2..])
}

/// A float32 sum turned into a centroid the way numpy finishes `.mean()`: the
/// divide is performed in f64 and rounded back to f32
/// (`true_divide(..., out=<f32>, casting='unsafe')`), then widened to f64, which
/// is what `.tolist()` hands the Python.
fn mean_from_f32_sum(sum: f32, n: usize) -> f64 {
    ((sum as f64 / n as f64) as f32) as f64
}

fn build<T: Sample>(
    r: &Raster<'_, T>,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
) -> (Regions, Vec<u32>) {
    let nbands = r.nbands;
    let mut idx: HashMap<u32, u32> = HashMap::new();
    let mut regs: Vec<Reg> = Vec::new();
    // Exact for the integer widths: every partial sum is an integer well below
    // 2^53, so order cannot matter. Left unused when `T::ORDERED_MEAN`.
    let mut sums: Vec<i64> = Vec::new();
    // The multi-band float case, which numpy accumulates sequentially: that is
    // exactly a streaming sum in raster order, so it needs no buffer at all.
    let mut fsums: Vec<f32> = Vec::new();
    let mut zeros: Vec<u32> = Vec::new();
    // Single-band float is the one case numpy sums pairwise, which cannot be
    // streamed -- it needs the region's samples laid out together.
    let gather = T::ORDERED_MEAN && nbands == 1;

    for y in 0..nlines {
        for x in 0..nsamps {
            let p = y * nsamps + x;
            let id = rband[p];
            if id == 0 {
                continue;
            }
            // Insertion order is raster-first-occurrence order, which is the
            // order every later loop visits regions in.
            let i = match idx.get(&id) {
                Some(&i) => i as usize,
                None => {
                    let i = regs.len();
                    idx.insert(id, i as u32);
                    regs.push(Reg {
                        id,
                        npix: 0,
                        ulx: x as u32,
                        uly: y as u32,
                        lrx: x as u32,
                        lry: y as u32,
                        alive: true,
                        nnbr_id: 0,
                        nnbr_d2: f64::INFINITY,
                        state: 0,
                    });
                    if !T::ORDERED_MEAN {
                        sums.resize((i + 1) * nbands, 0);
                    } else if !gather {
                        fsums.resize((i + 1) * nbands, 0.0);
                    }
                    zeros.push(0);
                    i
                }
            };
            let g = &mut regs[i];
            g.npix += 1;
            g.ulx = g.ulx.min(x as u32);
            g.uly = g.uly.min(y as u32);
            g.lrx = g.lrx.max(x as u32);
            g.lry = g.lry.max(y as u32);

            let pix = r.pixel_at(p);
            let mut all_zero = true;
            let o = i * nbands;
            for b in 0..nbands {
                let s = pix[b];
                if !T::ORDERED_MEAN {
                    sums[o + b] += s.to_f64() as i64;
                } else if !gather {
                    fsums[o + b] += s.to_f64() as f32;
                }
                all_zero &= s.is_zero();
            }
            if all_zero {
                zeros[i] += 1;
            }
        }
    }

    let mut ctr = vec![0.0f64; regs.len() * nbands];
    if gather {
        // Single-band float: numpy sums this one pairwise, so the region's
        // samples have to sit together, in `coords` (raster) order. One extra
        // pass and one `npixels` float buffer, both dropped before the merge
        // loop starts. Nothing else pays for this.
        let nreg = regs.len();
        let mut off = vec![0usize; nreg + 1];
        for i in 0..nreg {
            off[i + 1] = off[i] + regs[i].npix as usize;
        }
        let mut gath = vec![0.0f32; off[nreg]];
        let mut fill = vec![0usize; nreg];
        for y in 0..nlines {
            for x in 0..nsamps {
                let p = y * nsamps + x;
                let id = rband[p];
                if id == 0 {
                    continue;
                }
                let i = idx[&id] as usize;
                gath[off[i] + fill[i]] = r.pixel_at(p)[0].to_f64() as f32;
                fill[i] += 1;
            }
        }
        for i in 0..nreg {
            let n = regs[i].npix as usize;
            ctr[i] = mean_from_f32_sum(pairwise_sum_f32(&gath[off[i]..off[i + 1]]), n);
        }
    } else if T::ORDERED_MEAN {
        for i in 0..regs.len() {
            let n = regs[i].npix as usize;
            for b in 0..nbands {
                ctr[i * nbands + b] = mean_from_f32_sum(fsums[i * nbands + b], n);
            }
        }
    } else {
        for i in 0..regs.len() {
            let n = regs[i].npix as f64;
            for b in 0..nbands {
                ctr[i * nbands + b] = sums[i * nbands + b] as f64 / n;
            }
        }
    }

    // "More than half the pixels are nodata" -> the region is excluded, and its
    // pixels become region 0. This is where non-treed area enters: the stage-1
    // map need not carry a mask at all.
    //
    // Compared as `2*zeros > npix` rather than `zeros/npix > 0.5`; identical for
    // any region small enough to exist, and not at the mercy of a rounding.
    let mut dropped = Vec::new();
    for i in 0..regs.len() {
        if (zeros[i] as u64) * 2 > regs[i].npix as u64 {
            regs[i].alive = false;
            dropped.push(i as u32);
        }
    }

    (
        Regions {
            regs,
            ctr,
            nbands,
            idx,
        },
        dropped,
    )
}

/// Rewrite every pixel of region `id` inside `bbox` to `to`.
fn relabel(rband: &mut [u32], nsamps: usize, g: (u32, u32, u32, u32), id: u32, to: u32) {
    let (ulx, uly, lrx, lry) = g;
    for y in uly..=lry {
        let row = y as usize * nsamps;
        for x in ulx..=lrx {
            let p = row + x as usize;
            if rband[p] == id {
                rband[p] = to;
            }
        }
    }
}

/// Nearest neighbour of region `i` among the regions its pixels touch.
///
/// Candidates are visited in **ascending region id** and a near-tie **keeps the
/// incumbent**, so the smallest id among near-equal candidates wins. That is the
/// rule `tests/STAGE2.md` pins; the oracle's coin flip is not reachable under
/// it, so this phase consumes no randomness at all.
fn find_nearest(
    rl: &Regions,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
    i: usize,
    cands: &mut Vec<u32>,
) -> (u32, f64) {
    cands.clear();
    let g = &rl.regs[i];
    let (id, ulx, uly, lrx, lry) = (g.id, g.ulx, g.uly, g.lrx, g.lry);

    for y in uly..=lry {
        let row = y as usize * nsamps;
        for x in ulx..=lrx {
            let p = row + x as usize;
            if rband[p] != id {
                continue;
            }
            let mut push = |q: usize| {
                let v = rband[q];
                if v != 0 && v != id {
                    cands.push(v);
                }
            };
            if y > 0 {
                push(p - nsamps);
            }
            if (x as usize) + 1 < nsamps {
                push(p + 1);
            }
            if (y as usize) + 1 < nlines {
                push(p + nsamps);
            }
            if x > 0 {
                push(p - 1);
            }
        }
    }
    cands.sort_unstable();
    cands.dedup();

    let mut best_id = 0u32;
    let mut best_d2 = f64::INFINITY;
    for &c in cands.iter() {
        let Some(&j) = rl.idx.get(&c) else { continue };
        let j = j as usize;
        if !rl.regs[j].alive {
            continue;
        }
        let d = rl.dist2(i, j);
        if isclose(d, best_d2, TIE_REL_TOL) {
            // Tie: keep the incumbent, which is the lower id.
        } else if d < best_d2 {
            best_d2 = d;
            best_id = c;
        }
    }
    (best_id, best_d2)
}

/// Run segment development over `rband`, in place.
pub fn run(
    rband: &mut [u32],
    nlines: usize,
    nsamps: usize,
    img: &Image,
    cfg: &Stage2Config,
) -> Result<Stage2Result, String> {
    if cfg.nmin == 0 {
        return Err("stage 2 minimum region size must be > 0".into());
    }
    if cfg.nmax < cfg.nmin {
        return Err("stage 2 maximum region size must be >= the minimum".into());
    }
    if rband.len() != nlines * nsamps {
        return Err(format!(
            "region map is {} pixels, expected {}x{}",
            rband.len(),
            nlines,
            nsamps
        ));
    }
    if img.nlines != nlines || img.nsamps != nsamps {
        return Err(format!(
            "the second image is {}x{} but the region map is {}x{}; they must be \
             the same grid",
            img.nlines, img.nsamps, nlines, nsamps
        ));
    }

    let (mut rl, dropped) = img.with_raster(|r| match r {
        RasterRef::U8(r) => build(&r, rband, nlines, nsamps),
        RasterRef::U16(r) => build(&r, rband, nlines, nsamps),
        RasterRef::I16(r) => build(&r, rband, nlines, nsamps),
        RasterRef::F32(r) => build(&r, rband, nlines, nsamps),
    });

    for &i in &dropped {
        let g = &rl.regs[i as usize];
        relabel(rband, nsamps, (g.ulx, g.uly, g.lrx, g.lry), g.id, 0);
    }

    let mut nreg = rl.regs.iter().filter(|g| g.alive).count();
    let mut stats: Vec<PassStats2> = Vec::new();
    let mut cands: Vec<u32> = Vec::new();

    // The oracle's loop, including how it counts: `passes` ends up one greater
    // than the number executed, and that number is what names the output file.
    let mut old_nreg = 0usize;
    let mut passes = 1usize;
    while nreg != old_nreg {
        let mut st = PassStats2::default();

        // Only the distance is reset. The oracle resets `nearest_region_dist`
        // and leaves `nearest_region_id` standing from the previous pass, so a
        // region that finds no candidate this time still carries a stale
        // partner id -- which lands it in the `inf` bucket rather than
        // `no_cand`. That is only a counter difference, but it is the counter
        // that would otherwise send someone hunting the wrong pass.
        //
        // A stale id is always safe to follow: a distance is finite only if it
        // was written this pass, and both writers set the id with it.
        for g in rl.regs.iter_mut() {
            g.nnbr_d2 = f64::INFINITY;
            g.state = 0;
        }

        // Pass 1: every undersized region picks a partner, and *forces* the
        // partner to point back at it if it is the closest claimant so far.
        // This is the paper's relaxation -- A may merge with B even when A is
        // not B's nearest neighbour.
        for i in 0..rl.regs.len() {
            if !rl.regs[i].alive || rl.regs[i].npix >= cfg.nmin {
                continue;
            }
            let (best_id, best_d2) = find_nearest(&rl, rband, nlines, nsamps, i, &mut cands);
            if best_id == 0 {
                continue;
            }
            rl.regs[i].nnbr_id = best_id;
            rl.regs[i].nnbr_d2 = best_d2;
            let j = rl.idx[&best_id] as usize;
            if best_d2 < rl.regs[j].nnbr_d2 {
                rl.regs[j].nnbr_id = rl.regs[i].id;
                rl.regs[j].nnbr_d2 = best_d2;
            }
        }

        // Pass 2: merge, in the same order, at most one merge per region.
        //
        // A region absorbed *earlier in this pass* is still visited: the oracle
        // does not remove it from its region dict until the pass ends, so it
        // reaches the `busy` test and is counted there. Skipping it here would
        // give the same map but different per-pass counters, and the counters
        // are how a divergence gets localised.
        for i in 0..rl.regs.len() {
            if !rl.regs[i].alive && rl.regs[i].state != 2 {
                continue;
            }
            if rl.regs[i].npix >= cfg.nmin {
                continue;
            }
            st.considered += 1;
            if rl.regs[i].state != 0 {
                st.busy += 1;
                continue;
            }
            let pid = rl.regs[i].nnbr_id;
            if pid == 0 {
                st.no_cand += 1;
                continue;
            }
            let j = rl.idx[&pid] as usize;
            if rl.regs[j].state != 0 {
                st.busy += 1;
                continue;
            }
            let d = rl.regs[i].nnbr_d2;
            if d.is_infinite() {
                st.inf += 1;
                continue;
            }
            if u64::from(rl.regs[i].npix) + u64::from(rl.regs[j].npix) > u64::from(cfg.nmax) {
                st.over_max += 1;
                continue;
            }
            if !isclose(d, rl.regs[j].nnbr_d2, MUTUAL_REL_TOL) {
                st.not_mutual += 1;
                continue;
            }

            // i absorbs j, and keeps i's id -- the smaller region survives here,
            // which is the opposite of stage 1's lower-id rule.
            let (n1, n2) = (rl.regs[i].npix as f64, rl.regs[j].npix as f64);
            let nb = rl.nbands;
            let (oi, oj) = (i * nb, j * nb);
            for b in 0..nb {
                rl.ctr[oi + b] = (n1 * rl.ctr[oi + b] + n2 * rl.ctr[oj + b]) / (n1 + n2);
            }
            let jg = (
                rl.regs[j].ulx,
                rl.regs[j].uly,
                rl.regs[j].lrx,
                rl.regs[j].lry,
            );
            let (jid, iid) = (rl.regs[j].id, rl.regs[i].id);
            relabel(rband, nsamps, jg, jid, iid);

            rl.regs[i].ulx = rl.regs[i].ulx.min(jg.0);
            rl.regs[i].uly = rl.regs[i].uly.min(jg.1);
            rl.regs[i].lrx = rl.regs[i].lrx.max(jg.2);
            rl.regs[i].lry = rl.regs[i].lry.max(jg.3);
            rl.regs[i].npix += rl.regs[j].npix;
            rl.regs[i].state = 1;
            rl.regs[j].state = 2;
            rl.regs[j].alive = false;
            st.merged += 1;
            nreg -= 1;
        }

        st.nreg = nreg;
        stats.push(st);
        old_nreg = nreg + st.merged; // the count this pass started with
        passes += 1;
    }

    Ok(Stage2Result {
        passes,
        nreg,
        dropped_majority_nodata: dropped.len(),
        stats,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The same deterministic sequence the numpy side used to produce the
    /// constants below. `0.017` and `2.1` are deliberately *not* exact binary
    /// fractions -- an earlier version of this test used `0.375` and `5.0`,
    /// which are, so every summation order agreed and the test proved nothing.
    fn gen(n: usize) -> Vec<f32> {
        (0..n)
            .map(|i| ((i * i) % 97) as f32 * 0.017f32 - 2.1f32)
            .collect()
    }

    fn sequential(a: &[f32]) -> f32 {
        let mut r = 0.0f32;
        for &v in a {
            r += v;
        }
        r
    }

    /// numpy 2.4.3's own `np.add.reduce`, `mean(axis=1)` and sequential sums,
    /// recorded as raw f32 bit patterns so the comparison is exact. Sizes
    /// straddle every branch of the pairwise algorithm: below 8, the 8-way
    /// block, exactly 128, just past it, and the recursive split.
    const EXPECT: &[(usize, u32, u32, u32)] = &[
        (1, 0xc0066666, 0xc0066666, 0xc0066666),
        (2, 0xc085db22, 0xc005db22, 0xc085db22),
        (7, 0xc15272af, 0xbff08311, 0xc15272af),
        (8, 0xc166b850, 0xbfe6b850, 0xc166b851),
        (9, 0xc176e977, 0xbfdb7a31, 0xc176e978),
        (15, 0xc1c11cab, 0xbfcdfc72, 0xc1c11cab),
        (16, 0xc1cdb22c, 0xbfcdb22c, 0xc1cdb22c),
        (63, 0xc29f5603, 0xbfa1dd79, 0xc29f5602),
        (127, 0xc325c9b9, 0xbfa717e9, 0xc325c9ba),
        (128, 0xc3276dd2, 0xbfa76dd2, 0xc3276dd3),
        (129, 0xc3280872, 0xbfa6bafc, 0xc3280873),
        (130, 0xc3293709, 0xbfa69c97, 0xc329370b),
        (200, 0xc3826168, 0xbfa6e314, 0xc3826169),
        (255, 0xc3a31915, 0xbfa3bcd2, 0xc3a31916),
        (256, 0xc3a3d9ba, 0xbfa3d9ba, 0xc3a3d9ba),
        (257, 0xc3a461ca, 0xbfa3be0c, 0xc3a461ca),
        (1000, 0xc4a0d70a, 0xbfa4b33d, 0xc4a0d711),
        (4096, 0xc5a47e6c, 0xbfa47e6c, 0xc5a47e71),
        (5000, 0xc5c892f0, 0xbfa44f69, 0xc5c8930f),
    ];

    /// `pairwise_sum_f32` is numpy's contiguous (single-band) summation, bit for
    /// bit, and `mean_from_f32_sum` finishes it the way `.mean()` does.
    #[test]
    fn the_pairwise_sum_is_numpys_contiguous_sum() {
        for &(n, sum_bits, mean_bits, _) in EXPECT {
            let a = gen(n);
            let s = pairwise_sum_f32(&a);
            assert_eq!(s.to_bits(), sum_bits, "pairwise sum of {n} samples");
            let m = mean_from_f32_sum(s, n);
            assert_eq!((m as f32).to_bits(), mean_bits, "mean of {n} samples");
            // The centroid is the f32 mean widened, never an f64 mean.
            assert_eq!(m, f32::from_bits(mean_bits) as f64);
        }
    }

    /// The streaming accumulation the multi-band path uses is numpy's *strided*
    /// summation, which is plain sequential -- and genuinely different from the
    /// pairwise one. If these two ever agree everywhere, the band-count switch in
    /// `build` has stopped mattering and something is wrong.
    #[test]
    fn the_sequential_sum_is_numpys_strided_sum() {
        let mut differ = 0;
        for &(n, sum_bits, _, seq_bits) in EXPECT {
            let a = gen(n);
            assert_eq!(sequential(&a).to_bits(), seq_bits, "sequential sum of {n}");
            if seq_bits != sum_bits {
                differ += 1;
            }
        }
        assert!(
            differ > 0,
            "pairwise and sequential agreed on every size; the constants no \
             longer discriminate between numpy's two summation orders"
        );
    }

    /// The reason any of this exists: f32 accumulation is not an f64 sum, and on
    /// real data that difference moves the segmentation. If this starts passing
    /// trivially, the float path has quietly become an f64 accumulation.
    #[test]
    fn f32_and_f64_accumulation_genuinely_differ() {
        let a = gen(5000);
        let exact: f64 = a.iter().map(|&v| v as f64).sum::<f64>() / a.len() as f64;
        let oracle = mean_from_f32_sum(pairwise_sum_f32(&a), a.len());
        assert_ne!(oracle, exact, "f32 and f64 accumulation agreed");
        assert!(
            (oracle - exact).abs() < 1e-3,
            "they should differ in the last bits, not wildly: {oracle} vs {exact}"
        );
    }
}
