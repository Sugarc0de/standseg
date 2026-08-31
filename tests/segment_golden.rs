//! M3 + M4 gate: the whole segmenter, both phases, against both golden cases.
//!
//! This is the definition of done from PLAN.md. If it passes, the port is
//! faithful down to the last bit -- including the glibc `random()` tie-break,
//! since a single desynced draw would diverge the region map.

use std::path::{Path, PathBuf};

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, Observer, Phase};
use fast_segment::segment::{PassStats, Segmenter};

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(rel)
}

fn test_config() -> SegConfig {
    SegConfig {
        tols: vec![10.0],
        cm: 0.1,
        ..Default::default()
    }
    .with_n(&[15, 15, 100, 2500, 2500])
    .unwrap()
}

#[derive(Default)]
struct Capture {
    rmap: Option<(usize, Vec<u8>)>,
    armap: Option<(usize, Vec<u8>)>,
    passes: Vec<(Phase, usize, PassStats)>,
}

impl Observer for Capture {
    fn on_pass(&mut self, phase: Phase, pass: usize, s: &PassStats) {
        self.passes.push((phase, pass, *s));
    }
    fn on_map(&mut self, phase: Phase, pass: usize, seg: &Segmenter) -> Result<(), String> {
        let nbytes = seg.region_map_nbytes();
        let mut out = Vec::with_capacity(seg.nlines * seg.nsamps * nbytes);
        for &v in &seg.bands.rband {
            out.extend_from_slice(&v.to_le_bytes()[..nbytes]);
        }
        match phase {
            Phase::Normal => self.rmap = Some((pass, out)),
            Phase::Auxiliary => self.armap = Some((pass, out)),
            Phase::Stage2 => unreachable!("stage 2 was not asked for"),
        }
        Ok(())
    }
}

fn segment_case(input: &str) -> Capture {
    let img = fast_segment::io::read(&golden(input)).expect("read input");
    let mut cap = Capture::default();
    run(img, &test_config(), None, &mut cap).expect("segmentation failed");
    cap
}

/// The IPW containers carry a text header ahead of the payload; the `proof/`
/// files are the same payload raw. Either works as the comparison target.
fn ipw_payload(rel: &str) -> Vec<u8> {
    let raw = std::fs::read(golden(rel)).unwrap();
    raw[raw.len() - 125_000..].to_vec()
}

#[test]
fn case1_matches_golden_both_phases() {
    let cap = segment_case("misc/temp_byte_bip");

    let (rpass, rmap) = cap.rmap.expect("no rmap written");
    let (apass, armap) = cap.armap.expect("no armap written");

    assert_eq!(rpass, 51, "normal phase converged on the wrong pass");
    assert_eq!(apass, 58, "auxiliary phase converged on the wrong pass");

    assert_eq!(
        rmap,
        std::fs::read(golden("test_3456/expected/proof/regmap.rmap.51")).unwrap(),
        "Case 1 rmap does not match golden"
    );
    assert_eq!(
        armap,
        std::fs::read(golden("test_3456/expected/proof/regmap.armap.58")).unwrap(),
        "Case 1 armap does not match golden"
    );
}

#[test]
fn case2_matches_golden_both_phases() {
    let cap = segment_case("LC80220492014083LGN00/input/LC80220492014083LGN00_stack.ipw");

    let (rpass, rmap) = cap.rmap.expect("no rmap written");
    let (apass, armap) = cap.armap.expect("no armap written");

    assert_eq!(rpass, 17);
    assert_eq!(apass, 1);

    assert_eq!(
        rmap,
        ipw_payload("LC80220492014083LGN00/expected/t10-m1-n15_15_100_2500_2500_myseg.rmap.17"),
        "Case 2 rmap does not match golden"
    );
    assert_eq!(
        armap,
        ipw_payload("LC80220492014083LGN00/expected/t10-m1-n15_15_100_2500_2500_myseg.armap.1"),
        "Case 2 armap does not match golden"
    );
}

/// Spot-check the per-pass statistics against the values `myseg.log` records,
/// so a regression reports *which* pass broke rather than just "125000 bytes
/// differ".
#[test]
fn case1_pass_statistics_match_log() {
    let cap = segment_case("misc/temp_byte_bip");
    let normal: Vec<_> = cap
        .passes
        .iter()
        .filter(|(p, _, _)| *p == Phase::Normal)
        .collect();

    assert_eq!(normal.len(), 51);

    // Pass 1, straight from myseg.log.
    let (_, _, s1) = normal[0];
    assert_eq!(s1.nreg, 52664);
    assert_eq!(s1.maxpix, 4);
    assert_eq!(s1.merge_attempts, 52678);
    assert_eq!(s1.nnbr_gone, 2479);
    assert_eq!(s1.wrong_partner, 386);
    assert_eq!(s1.nnbr_d2_big, 47251);
    assert_eq!(s1.both_viable, 0);
    assert_eq!(s1.npix_big, 0);
    assert_eq!(s1.merging, 2562);

    // The final pass, which must produce no merges at all.
    let (_, _, s51) = normal[50];
    assert_eq!(s51.nreg, 22199);
    assert_eq!(s51.merging, 0);
    assert_eq!(s51.merge_attempts, 22199);
    assert_eq!(s51.wrong_partner, 13);
    assert_eq!(s51.both_viable, 4);
}

/// The parallel nearest-neighbour sweep must be bit-identical to the serial one.
///
/// The golden cases are far below the production threshold, so this forces the
/// parallel path on to exercise it. If the out-of-order collect ever perturbed
/// the `flip()` stream, this is where it would show.
#[test]
fn parallel_sweep_matches_golden_exactly() {
    let img = fast_segment::io::read(&golden("misc/temp_byte_bip")).expect("read");
    let cfg = SegConfig {
        par_threshold: 0, // force parallel
        ..test_config()
    };
    let mut cap = Capture::default();
    run(img, &cfg, None, &mut cap).expect("segmentation");

    let (rpass, rmap) = cap.rmap.expect("rmap");
    let (apass, armap) = cap.armap.expect("armap");
    assert_eq!((rpass, apass), (51, 58));
    assert_eq!(
        rmap,
        std::fs::read(golden("test_3456/expected/proof/regmap.rmap.51")).unwrap(),
        "parallel sweep diverged on rmap"
    );
    assert_eq!(
        armap,
        std::fs::read(golden("test_3456/expected/proof/regmap.armap.58")).unwrap(),
        "parallel sweep diverged on armap"
    );
}

/// Serial and parallel must agree on Case 2 as well, which has 8 bands.
#[test]
fn parallel_and_serial_agree_on_case2() {
    let path = golden("LC80220492014083LGN00/input/LC80220492014083LGN00_stack.ipw");
    let mut ser = Capture::default();
    run(
        fast_segment::io::read(&path).unwrap(),
        &SegConfig {
            threads: 1,
            ..test_config()
        },
        None,
        &mut ser,
    )
    .unwrap();
    let mut par = Capture::default();
    run(
        fast_segment::io::read(&path).unwrap(),
        &SegConfig {
            par_threshold: 0,
            ..test_config()
        },
        None,
        &mut par,
    )
    .unwrap();
    assert_eq!(ser.rmap, par.rmap);
    assert_eq!(ser.armap, par.armap);
}
