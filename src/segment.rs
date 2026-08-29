//! The segmentation driver: nearest-neighbour search, the two pass kinds, and
//! region-list compaction. Faithful port of `region.c` plus `segment.c`'s
//! `main_loop` / `wind_up`.

use crate::config::SegConfig;
use crate::contig::Connectivity;
use crate::nbrset::NbrSet;
use crate::pixel::Bands;
use crate::region::{
    merge_regions, RegionId, RegionList, RF_ACTIVE, RF_MERGE, RF_SPECIAL,
};
use crate::rng::GlibcRandom;

const N_DHISTBINS: usize = 1000;
const MAXFLOAT: f32 = f32::MAX;

/// Nearest-neighbour record, one per region id.
///
/// `id` doubles as the old->new translation table during compaction, exactly as
/// the C reuses `nnbrlist[].nbr_id`. At 15000^2 a separate table would be
/// another 800 MB. PLAN.md section 4, trick 5.
#[derive(Clone, Copy, Debug)]
pub struct Nbr {
    pub id: RegionId,
    pub d2: f32,
}

impl Default for Nbr {
    fn default() -> Self {
        Nbr { id: 0, d2: 0.0 }
    }
}

/// Per-pass statistics, mirroring the counters `myseg.log` prints.
#[derive(Debug, Default, Clone, Copy)]
pub struct PassStats {
    pub nreg: usize,
    pub dmin2: f32,
    pub maxpix: u32,
    pub merge_attempts: u64,
    pub special_merge_attempts: u64,
    pub nnbr_gone: u64,
    pub wrong_partner: u64,
    pub nnbr_d2_big: u64,
    pub both_viable: u64,
    pub npix_big: u64,
    pub merging: u64,
    pub no_nbr: u64,
    pub norminpix: u32,
    pub absminpix: u32,
    pub tp2: f32,
}

pub struct Segmenter<'a> {
    pub cfg: &'a SegConfig,
    pub bands: Bands,
    pub rl: RegionList,
    pub nnbr: Vec<Nbr>,
    pub nlines: usize,
    pub nsamps: usize,
    pub nreg: usize,
    pub maxreg: usize,
    pub rng: GlibcRandom,
    set: NbrSet,
    nbr_offsets: [isize; 8],
    dhbin: Vec<i64>,
    binwidth2: f32,
    use_hist: bool,
    tg2: f32,
    tp2: f32,
    /// Auxiliary region-map mask (`-A`).
    pub aband: Option<Vec<u8>>,
}

impl<'a> Segmenter<'a> {
    pub fn new(cfg: &'a SegConfig, bands: Bands, rl: RegionList, nlines: usize, nsamps: usize) -> Self {
        let nreg = bands.nreg;
        Self {
            cfg,
            bands,
            rl,
            nnbr: vec![Nbr::default(); nreg + 1],
            nlines,
            nsamps,
            nreg,
            maxreg: nreg,
            rng: GlibcRandom::new(),
            set: NbrSet::new(),
            nbr_offsets: cfg.conn.offsets(nsamps),
            dhbin: vec![0; N_DHISTBINS + 1],
            binwidth2: 0.0,
            use_hist: cfg.cm < 1.0,
            tg2: 0.0,
            tp2: 0.0,
            aband: None,
        }
    }

    #[inline]
    fn conn(&self) -> &Connectivity {
        &self.cfg.conn
    }

    /// Find region `rid`'s nearest neighbour and record it.
    ///
    /// The bounding-box scan is the algorithm, not a C workaround: no region
    /// ever stores a pixel list, so the box plus a `regid == rid` test is the
    /// only spatial index there is.
    fn reg_nnbr(&mut self, rid: RegionId) -> Result<(), String> {
        let conn = *self.conn();
        let b = self.rl.bbox[rid as usize];
        let ncdir = conn.ncdir;
        let offs = self.nbr_offsets;
        // Copy the direction table into a fixed-size array: `conn.flags` is a
        // slice of runtime length, so every `flags[d]` in the inner loop is a
        // bounds check that LLVM cannot hoist.
        let mut flags = [0u8; 8];
        flags[..ncdir].copy_from_slice(&conn.flags[..ncdir]);
        self.set.clear();

        for l in b.uly as usize..=b.lry as usize {
            let row = l * self.nsamps;
            let (lo, hi) = (row + b.ulx as usize, row + b.lrx as usize);
            // Walk the row as a pair of slices so the per-pixel reads are
            // bounds-checked once for the row rather than once per pixel.
            let rrow = &self.bands.rband[lo..=hi];
            let crow = &self.bands.cband[lo..=hi];
            for (i, (&pid, &cmap)) in rrow.iter().zip(crow.iter()).enumerate() {
                if pid != rid || cmap == conn.internal {
                    continue;
                }
                let p = lo + i;
                for d in 0..ncdir {
                    if cmap & flags[d] == 0 {
                        // Safe without a bounds test for the same reason the C
                        // is: out-of-bounds directions are marked contiguous.
                        let np = (p as isize + offs[d]) as usize;
                        let nbr = self.bands.rband[np];
                        if !self.set.add(nbr) {
                            return Err(format!(
                                "more than {} neighbors of region {rid}",
                                crate::nbrset::MAX_NEIGHBORS
                            ));
                        }
                    }
                }
            }
        }

        // Select the minimum, breaking exact ties with the glibc RNG. The tie
        // branch is the only consumer of randomness in a normal pass, and the
        // number of draws depends on the insertion order above.
        let mut mdist2 = MAXFLOAT;
        let mut nnbr: RegionId = 0;
        // Move the set out rather than collecting into a fresh Vec: this runs
        // once per region per pass, which is over 100M times per pass at
        // 15000^2, and an allocation there dominates everything else.
        let set = std::mem::replace(&mut self.set, NbrSet::empty());
        for &nbr in set.as_slice() {
            let ndist2 = self.rl.dist2(rid, nbr);
            if ndist2 > mdist2 {
                continue;
            } else if ndist2 < mdist2 {
                mdist2 = ndist2;
                nnbr = nbr;
            } else if self.rng.flip() {
                nnbr = nbr;
            }
        }

        self.set = set;
        self.nnbr[rid as usize] = Nbr { id: nnbr, d2: mdist2 };
        Ok(())
    }

    fn clear_d2hist(&mut self) {
        self.dhbin.iter_mut().for_each(|b| *b = 0);
    }

    fn hit_d2hist(&mut self, dist2: f32) {
        // C truncates a float into an int here. When dist2 is MAXFLOAT (a region
        // with no neighbours) that conversion is out of range; both the x86 and
        // aarch64 results land in the overflow bin, and Rust's saturating `as`
        // lands in the same place.
        let idx = (dist2 / self.binwidth2) as i32;
        let idx = if idx > N_DHISTBINS as i32 || idx < 0 {
            N_DHISTBINS
        } else {
            idx as usize
        };
        self.dhbin[idx] += 1;
    }

    fn get_tp2(&mut self) {
        let maxmerge = (self.nreg as f32 * self.cfg.cm) as i64;
        let mut index = 0usize;
        let mut cfreq = 0i64;
        while cfreq <= maxmerge && index <= N_DHISTBINS {
            cfreq += self.dhbin[index];
            index += 1;
        }
        if index > N_DHISTBINS {
            self.tp2 = self.tg2;
            self.use_hist = false;
        } else {
            self.tp2 = self.binwidth2 * index as f32;
        }
    }

    /// A normal pass. Returns the statistics `myseg.log` prints.
    pub fn seg_pass(&mut self) -> Result<PassStats, String> {
        let d2hist = self.use_hist;
        if d2hist {
            self.clear_d2hist();
        }

        for r in 1..=self.maxreg as RegionId {
            if !self.rl.is_active(r) {
                continue;
            }
            self.rl.flags[r as usize] &= !RF_MERGE;
            self.reg_nnbr(r)?;
            if d2hist {
                let d2 = self.nnbr[r as usize].d2;
                self.hit_d2hist(d2);
            }
        }
        if d2hist {
            self.get_tp2();
        }

        let mut st = PassStats {
            dmin2: MAXFLOAT,
            tp2: self.tp2,
            ..Default::default()
        };

        for r in 1..=self.maxreg as RegionId {
            if !self.rl.is_active(r) {
                continue;
            }
            let nnbr_id = self.nnbr[r as usize].id;
            if nnbr_id == 0 {
                st.no_nbr += 1;
                continue;
            }
            let nnbr_d2 = self.nnbr[r as usize].d2;
            st.dmin2 = st.dmin2.min(nnbr_d2);
            st.merge_attempts += 1;

            let nflags = self.rl.flags[nnbr_id as usize];
            let npix_r = self.rl.npix[r as usize] as u32;
            if nflags & RF_ACTIVE == 0 || nflags & RF_MERGE != 0 {
                st.maxpix = st.maxpix.max(npix_r);
                st.nnbr_gone += 1;
                continue;
            }

            if nnbr_d2 > self.tp2 {
                st.nnbr_d2_big += 1;
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }
            // C: fabs(float - float) <= FLT_EPSILON, evaluated as doubles.
            let other_d2 = self.nnbr[nnbr_id as usize].d2;
            if ((other_d2 - nnbr_d2) as f64).abs() > f32::EPSILON as f64 {
                st.wrong_partner += 1;
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }
            let npix_n = self.rl.npix[nnbr_id as usize] as u32;
            if npix_r >= self.cfg.nviable && npix_n >= self.cfg.nviable {
                st.both_viable += 1;
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }
            if npix_r + npix_n > self.cfg.nmax {
                st.npix_big += 1;
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }

            st.merging += 1;
            let (lo, hi) = if r < nnbr_id { (r, nnbr_id) } else { (nnbr_id, r) };
            merge_regions(
                &mut self.rl,
                &mut self.bands.rband,
                &mut self.bands.cband,
                self.nsamps,
                &self.cfg.conn,
                &self.nbr_offsets,
                lo,
                hi,
            )?;
            self.rl.flags[lo as usize] |= RF_MERGE;
            st.maxpix = st.maxpix.max(self.rl.npix[lo as usize] as u32);
            self.nreg -= 1;
        }

        st.nreg = self.nreg;
        Ok(st)
    }

    /// Garbage-collect the region list so active ids are contiguous from 1.
    pub fn compact_region_list(&mut self) {
        let nb = self.rl.nbands;
        let mut nrid: RegionId = 1;
        for crid in 0..=self.maxreg as RegionId {
            if !self.rl.is_active(crid) {
                // The neighbour list is now a translation table.
                self.nnbr[crid as usize].id = 0;
                continue;
            }
            if crid != nrid {
                let (c, n) = (crid as usize, nrid as usize);
                self.rl.bbox[n] = self.rl.bbox[c];
                self.rl.npix[n] = self.rl.npix[c];
                self.rl.flags[n] = self.rl.flags[c];
                for b in 0..nb {
                    self.rl.ctr[n * nb + b] = self.rl.ctr[c * nb + b];
                }
            }
            self.nnbr[crid as usize].id = nrid;
            nrid += 1;
        }
        let new_nreg = (nrid - 1) as usize;

        for p in 0..self.bands.rband.len() {
            self.bands.rband[p] = self.nnbr[self.bands.rband[p] as usize].id;
        }

        let slots = new_nreg + 1;
        self.rl.bbox.truncate(slots);
        self.rl.npix.truncate(slots);
        self.rl.flags.truncate(slots);
        self.rl.ctr.truncate(slots * nb);
        self.rl.bbox.shrink_to_fit();
        self.rl.npix.shrink_to_fit();
        self.rl.flags.shrink_to_fit();
        self.rl.ctr.shrink_to_fit();
        self.nnbr.truncate(slots);
        self.nnbr.shrink_to_fit();

        self.maxreg = new_nreg;
    }

    /// Set the current tolerance and reset the histogram for it.
    pub fn set_tolerance(&mut self, tg: f32) {
        self.tg2 = tg * tg;
        self.tp2 = self.tg2;
        if self.cfg.cm < 1.0 {
            self.use_hist = true;
            self.binwidth2 = self.tg2 / N_DHISTBINS as f32;
        }
    }

    /// The smallest pixel size that can hold every region id, as
    /// `GDAL_write_image` computes it.
    pub fn region_map_nbytes(&self) -> usize {
        let mut nbits = 0;
        let mut n = self.nreg;
        while n != 0 {
            nbits += 1;
            n >>= 1;
        }
        match nbits {
            0..=8 => 1,
            9..=16 => 2,
            _ => 4,
        }
    }
}

/// Mark every pixel of region `rid` in a byte band. Used by `-A`.
fn mark_reg_in_image(
    rl: &RegionList,
    rband: &[RegionId],
    band: &mut [u8],
    nsamps: usize,
    rid: RegionId,
    val: u8,
) {
    let b = rl.bbox[rid as usize];
    for l in b.uly as usize..=b.lry as usize {
        for s in b.ulx as usize..=b.lrx as usize {
            let p = l * nsamps + s;
            if rband[p] == rid {
                band[p] = val;
            }
        }
    }
}

impl<'a> Segmenter<'a> {
    /// An auxiliary pass: force undersized regions to merge, with no distance
    /// ceiling and the relaxed `nabsmax` size cap.
    pub fn seg_apass(&mut self) -> Result<PassStats, String> {
        // nbr_d2 is reused here as "closest approach by any subminimal region",
        // so it has to be reset for every slot, not just active ones.
        for r in 1..=self.maxreg {
            self.nnbr[r].d2 = MAXFLOAT;
        }

        let use_norb = self.cfg.norm_band.is_some();
        for r in 1..=self.maxreg as RegionId {
            if !self.rl.is_active(r) {
                continue;
            }
            self.rl.flags[r as usize] &= !RF_MERGE;

            if let Some(nb) = self.cfg.norm_band {
                let v = self.rl.ctr(r)[nb];
                if v < self.cfg.nblow || v > self.cfg.nbhigh {
                    self.rl.flags[r as usize] |= RF_SPECIAL;
                } else {
                    self.rl.flags[r as usize] &= !RF_SPECIAL;
                }
            }

            let special = self.rl.flags[r as usize] & RF_SPECIAL != 0;
            let npix = self.rl.npix[r as usize] as u32;
            let undersized = (special && npix < self.cfg.nabsmin)
                || (!special && npix < self.cfg.nnormin);
            if !undersized {
                continue;
            }

            self.reg_nnbr(r)?;
            let nnbr_id = self.nnbr[r as usize].id;
            if nnbr_id == 0 {
                continue;
            }
            let mine = self.nnbr[r as usize].d2;
            let theirs = self.nnbr[nnbr_id as usize].d2;
            self.nnbr[nnbr_id as usize].d2 = theirs.min(mine);
        }

        let mut st = PassStats {
            dmin2: MAXFLOAT,
            norminpix: u32::MAX,
            absminpix: u32::MAX,
            ..Default::default()
        };

        for r in 1..=self.maxreg as RegionId {
            if !self.rl.is_active(r) {
                continue;
            }
            let special = self.rl.flags[r as usize] & RF_SPECIAL != 0;
            let npix_r = self.rl.npix[r as usize] as u32;
            let floor = if special { self.cfg.nabsmin } else { self.cfg.nnormin };

            // `track` records the smallest surviving region of this kind.
            macro_rules! track {
                () => {
                    if special {
                        st.absminpix = st.absminpix.min(npix_r);
                    } else {
                        st.norminpix = st.norminpix.min(npix_r);
                    }
                };
            }

            if npix_r >= floor {
                track!();
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }

            let nnbr_id = self.nnbr[r as usize].id;
            if nnbr_id == 0 {
                st.no_nbr += 1;
                continue;
            }
            let nnbr_d2 = self.nnbr[r as usize].d2;
            st.dmin2 = st.dmin2.min(nnbr_d2);
            if special {
                st.special_merge_attempts += 1;
            } else {
                st.merge_attempts += 1;
            }

            let nflags = self.rl.flags[nnbr_id as usize];
            if nflags & RF_ACTIVE == 0 || nflags & RF_MERGE != 0 {
                track!();
                st.maxpix = st.maxpix.max(npix_r);
                st.nnbr_gone += 1;
                continue;
            }

            let other_d2 = self.nnbr[nnbr_id as usize].d2;
            if ((other_d2 - nnbr_d2) as f64).abs() > f32::EPSILON as f64 {
                st.wrong_partner += 1;
                track!();
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }

            let npix_n = self.rl.npix[nnbr_id as usize] as u32;
            if npix_r + npix_n > self.cfg.nabsmax {
                st.npix_big += 1;
                track!();
                st.maxpix = st.maxpix.max(npix_r);
                continue;
            }

            st.merging += 1;
            if self.cfg.armm {
                // Mark the region that is being absorbed. Ties consume a draw.
                let loser = if npix_r > npix_n || (npix_r == npix_n && self.rng.flip()) {
                    nnbr_id
                } else {
                    r
                };
                if let Some(ab) = self.aband.as_mut() {
                    mark_reg_in_image(&self.rl, &self.bands.rband, ab, self.nsamps, loser, 0);
                }
            }

            let (lo, hi) = if r < nnbr_id { (r, nnbr_id) } else { (nnbr_id, r) };
            merge_regions(
                &mut self.rl,
                &mut self.bands.rband,
                &mut self.bands.cband,
                self.nsamps,
                &self.cfg.conn,
                &self.nbr_offsets,
                lo,
                hi,
            )?;
            self.rl.flags[lo as usize] |= RF_MERGE;
            let merged = self.rl.npix[lo as usize] as u32;
            if special {
                st.absminpix = st.absminpix.min(merged);
            } else {
                st.norminpix = st.norminpix.min(merged);
            }
            st.maxpix = st.maxpix.max(merged);
            self.nreg -= 1;
        }

        let _ = use_norb;
        st.nreg = self.nreg;
        Ok(st)
    }
}
