//! M2 gate: Phase 0 must produce the initial region counts recorded in the
//! golden `myseg.log` -- "N of a possible 62500 regions are required".

use std::path::{Path, PathBuf};

use fast_segment::config::SegConfig;

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(rel)
}

/// The invocation both golden cases were produced with.
fn test_config() -> SegConfig {
    SegConfig {
        tols: vec![10.0],
        cm: 0.1,
        ..Default::default()
    }
    .with_n(&[15, 15, 100, 2500, 2500])
    .unwrap()
}

fn nreg_for(input: &str) -> usize {
    let img = fast_segment::io::read(&golden(input)).expect("read input");
    let (bands, _rl) = fast_segment::pixel::phase0(&img, &test_config(), None).expect("phase0");
    bands.nreg
}

#[test]
fn case1_initial_region_count() {
    assert_eq!(
        nreg_for("misc/temp_byte_bip"),
        55226,
        "Case 1 Phase 0 region count does not match myseg.log"
    );
}

#[test]
fn case2_initial_region_count() {
    assert_eq!(
        nreg_for("LC80220492014083LGN00/input/LC80220492014083LGN00_stack.ipw"),
        31609,
        "Case 2 Phase 0 region count does not match myseg.log"
    );
}

/// Whichever container Case 1 is read from, Phase 0 must agree.
#[test]
fn case1_same_from_either_container() {
    assert_eq!(
        nreg_for("misc/temp_byte_bip"),
        nreg_for("test_3456/input/test_3456.bip.ipw")
    );
}
