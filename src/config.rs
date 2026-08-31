//! Segmentation parameters, validated exactly as `main.c` validates them.

use crate::contig::{Connectivity, EIGHT_WAY, FOUR_WAY};
use crate::region::MAX_REGION_PIXELS;

#[derive(Debug, Clone)]
pub struct SegConfig {
    /// Final tolerances (`-t`), consumed left to right.
    pub tols: Vec<f32>,
    /// Merge coefficient (`-m`). The C reads this uninitialised when `-m` is
    /// absent; we default it to 1.0, which is the documented "no restriction".
    pub cm: f32,
    /// `-n Nabsmin,Nnormin,Nviable,Nmax,Nabsmax`
    pub nabsmin: u32,
    pub nnormin: u32,
    pub nviable: u32,
    pub nmax: u32,
    pub nabsmax: u32,
    /// Normality band and interval (`-B`, `-N`). A region whose centroid in
    /// `norm_band` falls outside `[nblow, nbhigh]` is *special*, and is held to
    /// `nabsmin` pixels in Phase 2 rather than `nnormin`. Band index is
    /// zero-based, as it is in the C.
    pub norm_band: Option<usize>,
    pub nblow: f32,
    pub nbhigh: f32,
    /// Auxiliary region map mask (`-A`): record which side of each Phase 2
    /// merge was absorbed, and write it out beside the armap.
    pub armm: bool,
    pub conn: Connectivity,
    /// Region count above which the nearest-neighbour sweep goes parallel.
    /// Below it the fan-out costs more than the scan. 0 forces parallel on.
    pub par_threshold: usize,
    /// 0 = rayon's default (one per core), 1 = force the serial path.
    pub threads: usize,
}

impl Default for SegConfig {
    fn default() -> Self {
        Self {
            tols: vec![],
            cm: 1.0,
            nabsmin: 1,
            nnormin: 1,
            // "No limit". The C spelled this 65535 because `npix` was an
            // `unsigned short`; it is not a limit any more (section 12.2).
            nviable: MAX_REGION_PIXELS,
            nmax: MAX_REGION_PIXELS,
            nabsmax: MAX_REGION_PIXELS,
            norm_band: None,
            nblow: 0.0,
            nbhigh: 255.0,
            armm: false,
            conn: FOUR_WAY,
            par_threshold: 200_000,
            threads: 0,
        }
    }
}

impl SegConfig {
    pub fn eight_way(mut self, yes: bool) -> Self {
        self.conn = if yes { EIGHT_WAY } else { FOUR_WAY };
        self
    }

    /// Apply `-B band` and `-N low,high`. The C requires them together and
    /// checks `0 <= low < high`; its further `high <= 255` is a uint8-era
    /// bound, so the caller range-checks against the actual input instead.
    pub fn with_normality(mut self, band: usize, low: f32, high: f32) -> Result<Self, String> {
        // Deliberately a negated `<` rather than `>=`: if either bound is NaN
        // every comparison is false, and this way that lands in the error arm
        // instead of sailing through. clippy would rather see `partial_cmp`.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(low < high) {
            return Err("normality interval (-N low,high) must have low < high".into());
        }
        if low < 0.0 {
            return Err("normality interval (-N low,high) must have low >= 0".into());
        }
        self.norm_band = Some(band);
        self.nblow = low;
        self.nbhigh = high;
        Ok(self)
    }

    /// Apply `-n`, left to right, with the C's "0 means default" rule.
    pub fn with_n(mut self, n: &[u32]) -> Result<Self, String> {
        if !n.is_empty() {
            self.nabsmin = if n[0] == 0 { 1 } else { n[0] };
        }
        if n.len() >= 2 {
            self.nnormin = if n[1] == 0 { self.nabsmin } else { n[1] };
            if self.nnormin < self.nabsmin {
                return Err("Nnormin (-n Nabsmin,Nnormin) must be >= Nabsmin".into());
            }
        }
        if n.len() >= 3 {
            self.nviable = if n[2] == 0 { MAX_REGION_PIXELS } else { n[2] };
            if self.nviable < self.nnormin {
                return Err("Nviable must be >= Nnormin".into());
            }
        }
        if n.len() >= 4 {
            self.nmax = if n[3] == 0 { MAX_REGION_PIXELS } else { n[3] };
            if self.nmax < self.nviable {
                return Err("Nmax must be >= Nviable".into());
            }
        }
        if n.len() >= 5 {
            self.nabsmax = if n[4] == 0 { MAX_REGION_PIXELS } else { n[4] };
            if self.nabsmax < self.nmax {
                return Err("Nabsmax must be >= Nmax".into());
            }
        }
        Ok(self)
    }
}
