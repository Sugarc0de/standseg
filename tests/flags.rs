//! The options that were half-present: `-B`/`-N` (normality band) and `-A`
//! (auxiliary region map mask).
//!
//! All three had working logic in the engine and no way to reach it from the
//! command line, and `-A` computed its mask and then threw it away. `-b` and
//! `-l` -- the per-pass single-band debug files -- were never implemented at
//! all, and their config fields have been deleted rather than left looking
//! like features.

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, Observer, Phase};
use fast_segment::image::{Image, Samples};
use fast_segment::region::RegionId;
use fast_segment::segment::Segmenter;

#[derive(Default)]
struct Capture {
    armap: Option<Vec<RegionId>>,
    armask: Option<Vec<u8>>,
    nreg: usize,
}

impl Observer for Capture {
    fn on_map(&mut self, phase: Phase, _pass: usize, seg: &Segmenter) -> Result<(), String> {
        if phase == Phase::Auxiliary {
            self.armap = Some(seg.bands.rband.clone());
            self.armask = seg.aband.clone();
            self.nreg = seg.nreg;
        }
        Ok(())
    }
}

/// A bright field with 25 isolated dark pixels. Each dark pixel is its own
/// region -- 190 DN from its neighbours, far past the tolerance.
fn speckled() -> Image {
    let (nl, ns) = (20usize, 20usize);
    let mut v = vec![200u8; nl * ns];
    for l in (1..nl).step_by(4) {
        for s in (1..ns).step_by(4) {
            v[l * ns + s] = 10;
        }
    }
    Image::from_samples(nl, ns, 1, Samples::U8(v))
}

fn base_cfg() -> SegConfig {
    // Nabsmin 1, Nnormin 25: a normal region under 25 pixels gets force-merged
    // in Phase 2, a special one under 1 pixel never does.
    SegConfig {
        tols: vec![10.0],
        cm: 1.0,
        ..Default::default()
    }
    .with_n(&[1, 25, 0, 0, 0])
    .unwrap()
}

fn segment(img: Image, cfg: &SegConfig) -> Capture {
    let mut cap = Capture::default();
    run(img, cfg, None, &mut cap).expect("segmentation failed");
    cap
}

/// Without `-B`/`-N` the dark specks are ordinary undersized regions and
/// Phase 2 absorbs them. With them they are *special*, held to Nabsmin instead
/// of Nnormin, and they survive. That is the whole point of the option, and it
/// was unreachable from the command line.
#[test]
fn a_normality_band_keeps_special_regions_that_would_be_absorbed() {
    let plain = segment(speckled(), &base_cfg());
    let norb = segment(
        speckled(),
        &base_cfg().with_normality(0, 50.0, 255.0).unwrap(),
    );

    assert!(
        norb.nreg > plain.nreg,
        "normality band should preserve regions: {} with, {} without",
        norb.nreg,
        plain.nreg
    );
    assert_ne!(
        norb.armap, plain.armap,
        "the normality band did not change the segmentation at all"
    );

    // The dark specks specifically: each should still be its own region.
    let armap = norb.armap.unwrap();
    let ns = 20usize;
    // row * width + col; clippy would rather see the folded constant.
    #[allow(clippy::identity_op)]
    let speck = armap[1 * ns + 1];
    let field = armap[0];
    assert_ne!(
        speck, field,
        "a dark speck was absorbed into the bright field"
    );
}

/// An interval that contains every centroid marks nothing special, so the run
/// must match the one with no `-B` at all.
#[test]
fn an_all_inclusive_interval_changes_nothing() {
    let plain = segment(speckled(), &base_cfg());
    let wide = segment(
        speckled(),
        &base_cfg().with_normality(0, 0.0, 255.0).unwrap(),
    );
    assert_eq!(plain.armap, wide.armap);
    assert_eq!(plain.nreg, wide.nreg);
}

/// `-A` used to allocate the mask, fill it in during Phase 2, and drop it.
#[test]
fn the_auxiliary_mask_records_what_phase_2_absorbed() {
    let cfg = SegConfig {
        armm: true,
        ..base_cfg()
    };
    let cap = segment(speckled(), &cfg);

    let mask = cap.armask.expect("-A should produce a mask");
    assert_eq!(mask.len(), 20 * 20);
    assert!(
        mask.iter().all(|&v| v == 0 || v == 1),
        "the mask must be binary"
    );
    assert!(
        mask.contains(&0),
        "Phase 2 absorbed regions here, so the mask cannot be all ones"
    );
    assert!(mask.contains(&1), "the mask cannot be all zeros");

    // Without -A there is no mask at all -- no wasted allocation.
    assert!(segment(speckled(), &base_cfg()).armask.is_none());
}

#[test]
fn normality_interval_is_validated() {
    assert!(SegConfig::default().with_normality(0, 200.0, 50.0).is_err());
    assert!(SegConfig::default().with_normality(0, 50.0, 50.0).is_err());
    assert!(SegConfig::default().with_normality(0, -1.0, 50.0).is_err());
    let c = SegConfig::default().with_normality(2, 50.0, 200.0).unwrap();
    assert_eq!(c.norm_band, Some(2));
    assert_eq!((c.nblow, c.nbhigh), (50.0, 200.0));
}
