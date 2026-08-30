//! Top-level control flow: `segment.c`'s `main_loop` and `wind_up`, with the
//! two `exit(0)`s removed so both phases always run and both maps are written.

use crate::config::SegConfig;
use crate::image::Image;
use crate::pixel::phase0;
use crate::segment::{PassStats, Segmenter};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Normal passes -- produce `.rmap.<n>`.
    Normal,
    /// Auxiliary passes -- produce `.armap.<n>`.
    Auxiliary,
}

/// Hooks so the caller can log passes and write maps without the driver ever
/// holding a second copy of the region band (which is 900 MB at 15000^2).
pub trait Observer {
    fn on_start(&mut self, _nreg: usize, _npix: usize) {}
    fn on_memory(&mut self, _report: &MemReport) {}
    fn on_pass(&mut self, _phase: Phase, _pass: usize, _stats: &PassStats) {}
    fn on_no_merges(&mut self, _pass: usize) {}
    fn on_map(&mut self, _phase: Phase, _pass: usize, _seg: &Segmenter) -> Result<(), String> {
        Ok(())
    }
}

/// The C computes this from `getpagesize()`, which is 4096 on the Linux box that
/// produced the fixtures but 16384 on Apple silicon. Compaction only renumbers
/// active regions in place, order-preserving, so its timing cannot change the
/// result -- but pin the Linux value anyway so behaviour does not drift by host.
const PAGESIZE: usize = 4096;
const MIN_RECLAIM_PAGES: usize = 8;

fn reclaim_trigger(nbands: usize) -> usize {
    // sizeof(region) + sizeof(neighbor) + sizeof(float) * nbands, as the C does.
    (PAGESIZE * MIN_RECLAIM_PAGES) / (12 + 8 + 4 * nbands) + 1
}

/// Exact sizes of the arrays that dominate the footprint, so the section 6
/// budget can be checked against reality rather than sampled with `ps`.
#[derive(Debug, Clone, Copy)]
pub struct MemReport {
    pub npixels: usize,
    pub nreg: usize,
    pub nbands: usize,
    pub image: usize,
    pub cband: usize,
    pub rband: usize,
    pub rlist: usize,
    pub ctrlist: usize,
    pub nnbrlist: usize,
}

impl MemReport {
    /// Peak while the image is still resident (during Phase 0).
    pub fn peak_phase0(&self) -> usize {
        self.image + self.cband + self.rband + self.rlist + self.ctrlist
    }
    /// Peak once the image is freed and the neighbour list exists.
    pub fn peak_phase1(&self) -> usize {
        self.cband + self.rband + self.rlist + self.ctrlist + self.nnbrlist
    }
    pub fn peak(&self) -> usize {
        self.peak_phase0().max(self.peak_phase1())
    }
}

pub struct RunResult {
    pub normal_passes: usize,
    pub aux_passes: usize,
    pub final_nreg: usize,
}

/// Takes the image **by value** and drops it as soon as Phase 0 is done, which
/// is what `segment.c` does with `free_image()`. Nothing after Phase 0 reads
/// pixels -- the centroids carry everything -- and on a 15000^2 x 6 scene the
/// buffer is 1.35 GB, so holding it would inflate peak memory by a fifth.
pub fn run(
    img: Image,
    cfg: &SegConfig,
    mask: Option<&[u8]>,
    obs: &mut dyn Observer,
) -> Result<RunResult, String> {
    if cfg.tols.is_empty() {
        return Err("at least one final tolerance (-t tol) required".into());
    }

    let (nlines, nsamps, nbands, npixels) = (img.nlines, img.nsamps, img.nbands, img.npixels());

    let (bands, rl) = phase0(&img, cfg, mask)?;
    obs.on_start(bands.nreg, npixels);
    let nreg0 = bands.nreg;
    obs.on_memory(&MemReport {
        npixels,
        nreg: nreg0,
        nbands,
        image: img.nbytes(),
        cband: npixels,
        rband: npixels * 4,
        // 13 bytes per region: BBox(8) + npix(4) + flags(1).
        rlist: (nreg0 + 2) * 13,
        ctrlist: (nreg0 + 2) * nbands * 4,
        nnbrlist: (nreg0 + 1) * 8,
    });
    drop(img);

    let mut seg = Segmenter::new(cfg, bands, rl, nlines, nsamps);
    let trigger = reclaim_trigger(nbands);

    // --- Phase 1: normal passes, one run per tolerance -------------------
    let mut pass = 0usize;
    for &tol in &cfg.tols {
        seg.set_tolerance(tol);
        loop {
            let old_nreg = seg.nreg;
            pass += 1;
            let stats = seg.seg_pass()?;
            // The C logs every pass including the final no-merge one; it just
            // does it after the break test (`log_pass(Spr, FALSE, FALSE)`).
            obs.on_pass(Phase::Normal, pass, &stats);
            if old_nreg == seg.nreg {
                obs.on_no_merges(pass);
                break;
            }
            if seg.maxreg - seg.nreg >= trigger {
                seg.compact_region_list();
            }
        }
        seg.compact_region_list();
        obs.on_map(Phase::Normal, pass, &seg)?;
    }

    // --- Phase 2: auxiliary passes ---------------------------------------
    // The C bails out here when nnormin == 1, and again under `-S breakpoint`.
    // Both are removed: we always run Phase 2 and always write the armap.
    if cfg.armm {
        seg.aband = Some(vec![1u8; npixels]);
    }

    let mut apass = 0usize;
    let mut old_nreg = usize::MAX;
    while old_nreg != seg.nreg {
        old_nreg = seg.nreg;
        apass += 1;
        let stats = seg.seg_apass()?;
        obs.on_pass(Phase::Auxiliary, apass, &stats);
        if seg.maxreg - seg.nreg >= trigger {
            seg.compact_region_list();
        }
    }

    seg.compact_region_list();
    obs.on_map(Phase::Auxiliary, apass, &seg)?;

    Ok(RunResult {
        normal_passes: pass,
        aux_passes: apass,
        final_nreg: seg.nreg,
    })
}
