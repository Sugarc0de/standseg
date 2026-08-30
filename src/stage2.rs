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

use std::collections::HashMap;

use crate::image::{Image, RasterRef, Raster, Sample};

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
/// The means are exactly reproducible against numpy without imitating its
/// pairwise summation: the samples are integers, so every partial sum is an
/// exact integer below 2^53 and summation order cannot matter. One `i64` sum
/// and one `f64` divide gives the same bits.
fn build<T: Sample>(
    r: &Raster<'_, T>,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
) -> (Regions, Vec<u32>) {
    let nbands = r.nbands;
    let mut idx: HashMap<u32, u32> = HashMap::new();
    let mut regs: Vec<Reg> = Vec::new();
    let mut sums: Vec<i64> = Vec::new();
    let mut zeros: Vec<u32> = Vec::new();

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
                    sums.resize((i + 1) * nbands, 0);
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
                let v = pix[b].to_i64();
                sums[o + b] += v;
                all_zero &= v == 0;
            }
            if all_zero {
                zeros[i] += 1;
            }
        }
    }

    let mut ctr = vec![0.0f64; regs.len() * nbands];
    for i in 0..regs.len() {
        let n = regs[i].npix as f64;
        for b in 0..nbands {
            ctr[i * nbands + b] = sums[i * nbands + b] as f64 / n;
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

    (Regions { regs, ctr, nbands, idx }, dropped)
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

        for g in rl.regs.iter_mut() {
            g.nnbr_id = 0;
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
        for i in 0..rl.regs.len() {
            if !rl.regs[i].alive || rl.regs[i].npix >= cfg.nmin {
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
            let jg = (rl.regs[j].ulx, rl.regs[j].uly, rl.regs[j].lrx, rl.regs[j].lry);
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
