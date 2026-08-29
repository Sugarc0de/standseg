//! Segmentation parameters, validated exactly as `main.c` validates them.

use crate::contig::{Connectivity, EIGHT_WAY, FOUR_WAY};
use crate::region::MAX_USHORT;

#[derive(Debug, Clone)]
pub struct SegConfig {
    /// Final tolerances (`-t`), consumed left to right.
    pub tols: Vec<f32>,
    /// Merge coefficient (`-m`). The C reads this uninitialised when `-m` is
    /// absent; we default it to 1.0, which is the documented "no restriction".
    pub cm: f32,
    /// Log threshold base and increment (`-l`).
    pub lthr: f32,
    pub lincr: f32,
    /// `-n Nabsmin,Nnormin,Nviable,Nmax,Nabsmax`
    pub nabsmin: u32,
    pub nnormin: u32,
    pub nviable: u32,
    pub nmax: u32,
    pub nabsmax: u32,
    /// Normality band and interval (`-B`, `-N`).
    pub norm_band: Option<usize>,
    pub nblow: f32,
    pub nbhigh: f32,
    /// Log band (`-b`).
    pub log_band: Option<usize>,
    /// Auxiliary region map mask (`-A`).
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
            lthr: 0.0,
            lincr: 0.0,
            nabsmin: 1,
            nnormin: 1,
            nviable: MAX_USHORT,
            nmax: MAX_USHORT,
            nabsmax: MAX_USHORT,
            norm_band: None,
            nblow: 0.0,
            nbhigh: 255.0,
            log_band: None,
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
            self.nviable = if n[2] == 0 { MAX_USHORT } else { n[2] };
            if self.nviable > MAX_USHORT || self.nviable < self.nnormin {
                return Err("Nviable must be Nnormin <= Nviable <= 65535".into());
            }
        }
        if n.len() >= 4 {
            self.nmax = if n[3] == 0 { MAX_USHORT } else { n[3] };
            if self.nmax < self.nviable || self.nmax > MAX_USHORT {
                return Err("Nmax must be Nviable <= Nmax <= 65535".into());
            }
        }
        if n.len() >= 5 {
            self.nabsmax = if n[4] == 0 { MAX_USHORT } else { n[4] };
            if self.nabsmax < self.nmax || self.nabsmax > MAX_USHORT {
                return Err("Nabsmax must be Nmax <= Nabsmax <= 65535".into());
            }
        }
        Ok(self)
    }
}
