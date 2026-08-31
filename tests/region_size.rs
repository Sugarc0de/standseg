//! Region sizes above the C's 65535-pixel ceiling.
//!
//! `npix` was an `unsigned short` in 1992, and the "no limit" settings for
//! `-n Nviable,Nmax,Nabsmax` were spelled 65535 for that reason. At 1 m
//! resolution 65535 pixels is a 256 m square, which is smaller than plenty of
//! real forest stands, so the ceiling had become a silent cap on the answer.

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, Observer, Phase};
use fast_segment::image::{Image, Samples};
use fast_segment::region::{RegionId, MAX_REGION_PIXELS};
use fast_segment::segment::Segmenter;

#[derive(Default)]
struct Capture {
    armap: Option<Vec<RegionId>>,
    npix: Vec<u32>,
    nreg: usize,
}

impl Observer for Capture {
    fn on_map(&mut self, phase: Phase, _pass: usize, seg: &Segmenter) -> Result<(), String> {
        if phase == Phase::Auxiliary {
            self.armap = Some(seg.bands.rband.clone());
            self.npix = seg.rl.npix.clone();
            self.nreg = seg.nreg;
        }
        Ok(())
    }
}

/// A uniform 300 x 300 field is one stand: 90000 pixels, well past 65535.
/// Before the widening, `Nmax`/`Nabsmax` defaulted to 65535 and the merges
/// simply stopped there, leaving the field split into arbitrary pieces.
#[test]
fn a_uniform_field_becomes_one_region_of_90000_pixels() {
    let (nl, ns) = (300usize, 300usize);
    let img = Image::from_samples(nl, ns, 1, Samples::U8(vec![100u8; nl * ns]));

    let cfg = SegConfig {
        tols: vec![10.0],
        cm: 1.0,
        ..Default::default()
    };
    let mut cap = Capture::default();
    run(img, &cfg, None, &mut cap).expect("segmentation failed");

    let armap = cap.armap.expect("no armap");
    assert_eq!(cap.nreg, 1, "a uniform field should end as a single region");
    let id = armap[0];
    assert!(armap.iter().all(|&r| r == id), "field did not fully merge");
    assert_eq!(
        cap.npix[id as usize],
        (nl * ns) as u32,
        "region pixel count did not survive past 65535"
    );
    assert!(cap.npix[id as usize] > 65535);
}

/// The old ceiling is still available as an explicit `-n` setting, and still
/// binds: asking for Nmax = Nabsmax = 65535 caps the same field below 90000.
#[test]
fn an_explicit_65535_ceiling_still_caps() {
    let (nl, ns) = (300usize, 300usize);
    let img = Image::from_samples(nl, ns, 1, Samples::U8(vec![100u8; nl * ns]));

    let cfg = SegConfig {
        tols: vec![10.0],
        cm: 1.0,
        ..Default::default()
    }
    .with_n(&[1, 1, 65535, 65535, 65535])
    .unwrap();
    let mut cap = Capture::default();
    run(img, &cfg, None, &mut cap).expect("segmentation failed");

    assert!(
        cap.nreg > 1,
        "an explicit 65535 ceiling should prevent one region"
    );
    assert!(
        cap.npix.iter().all(|&n| n <= 65535),
        "a region exceeded the ceiling that was asked for"
    );
}

/// `-n` no longer refuses values above 65535, and "0 means no limit" now means
/// no limit rather than 65535.
#[test]
fn n_accepts_values_above_65535() {
    let cfg = SegConfig::default()
        .with_n(&[15, 15, 100, 250_000, 250_000])
        .unwrap();
    assert_eq!(cfg.nmax, 250_000);
    assert_eq!(cfg.nabsmax, 250_000);

    let cfg = SegConfig::default().with_n(&[15, 15, 0, 0, 0]).unwrap();
    assert_eq!(cfg.nviable, MAX_REGION_PIXELS);
    assert_eq!(cfg.nmax, MAX_REGION_PIXELS);
    assert_eq!(cfg.nabsmax, MAX_REGION_PIXELS);

    // Ordering is still enforced.
    assert!(SegConfig::default().with_n(&[15, 15, 100, 50, 50]).is_err());
}

/// More than 5000 neighbours on one region.
///
/// The C carried `MAX_NEIGHBORS = 5000` and aborted the entire run on the
/// 5001st -- `add_to_set` returned FALSE and `reg_nnbr` gave up. A long thin
/// stand bordered by noise reaches that easily: here a 2501-pixel row is
/// flanked above and below by 5002 singleton regions.
#[test]
fn a_region_with_more_than_5000_neighbours_completes() {
    let (nl, ns) = (3usize, 2501usize);
    let mut v = vec![0u16; nl * ns];
    // Rows 0 and 2: every pixel its own region, 20 DN apart so a tolerance of
    // 10 never merges them with each other or with the row between.
    for s in 0..ns {
        let x = 10_000 + 20 * s as u16;
        v[s] = x;
        v[2 * ns + s] = x;
    }
    // Row 1 is uniform 0: one stand, 2501 pixels, 5002 neighbours.
    let img = Image::from_samples(nl, ns, 1, Samples::U16(v));

    let cfg = SegConfig {
        tols: vec![10.0],
        cm: 1.0,
        ..Default::default()
    };
    let mut cap = Capture::default();
    run(img, &cfg, None, &mut cap).expect("a region with >5000 neighbours must not abort the run");

    let armap = cap.armap.expect("no armap");
    let mid = armap[ns];
    assert!(
        armap[ns..2 * ns].iter().all(|&r| r == mid),
        "the middle row should be a single region"
    );
    assert_eq!(cap.npix[mid as usize], ns as u32);

    // Count its distinct neighbours to show the old ceiling really was crossed.
    let mut nbrs: std::collections::HashSet<RegionId> = std::collections::HashSet::new();
    for s in 0..ns {
        nbrs.insert(armap[s]);
        nbrs.insert(armap[2 * ns + s]);
    }
    assert!(
        nbrs.len() > 5000,
        "expected more than 5000 neighbouring regions, got {}",
        nbrs.len()
    );
}
