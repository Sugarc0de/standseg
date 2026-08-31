//! `-t` is a distance in raw DN units, so parameters tuned on 8-bit imagery are
//! silently wrong on 16-bit imagery: the spectral phase merges nothing, the size
//! rules force merges anyway, and the result is a plausible-looking map shaped by
//! region size rather than by the image. That is the worst kind of failure -- it
//! produces output. These tests pin the warning that catches it, and, just as
//! importantly, pin that it stays quiet on the real cases.

use std::path::{Path, PathBuf};
use std::process::Command;

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(rel)
}

fn outdir() -> PathBuf {
    let d = std::env::temp_dir().join("fast_segment_tolwarn");
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Returns stderr from a run of the CLI.
fn run(input: &Path, tol: &str, base: &str) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_fast_segment"))
        .arg(input)
        .args([
            "-t",
            tol,
            "-m",
            ".1",
            "-n",
            "15,15,100,2500,2500",
            "-o",
            base,
        ])
        .arg("--outdir")
        .arg(outdir())
        .output()
        .expect("run fast_segment");
    assert!(out.status.success(), "segmentation failed");
    String::from_utf8_lossy(&out.stderr).into_owned()
}

const WARNING: &str = "merged almost nothing";

/// The 16-bit Landsat stack runs 0..8990. A tolerance of 10 is an 8-bit value:
/// the first pass leaves essentially every pixel its own region.
#[test]
fn warns_when_an_8bit_tolerance_meets_16bit_data() {
    let stack = golden("LC80220492014083LGN00/input/LC80220492014083LGN00_stack");
    let err = run(&stack, "10", "tolwarn_bad");
    assert!(err.contains(WARNING), "expected a warning, got:\n{err}");
    // The message has to be actionable, not just alarming: it should say what
    // units -t is in and what range the image actually occupies.
    assert!(err.contains("DN units"), "no units in the message:\n{err}");
    assert!(
        err.contains("8990"),
        "no observed range in the message:\n{err}"
    );
}

/// The same scene with the tolerance scaled to its range must be silent.
#[test]
fn stays_quiet_when_the_tolerance_matches_the_data() {
    let stack = golden("LC80220492014083LGN00/input/LC80220492014083LGN00_stack");
    let err = run(&stack, "350", "tolwarn_good");
    assert!(
        !err.contains(WARNING),
        "false positive on scaled tolerance:\n{err}"
    );
}

/// The case that matters most: the golden 8-bit input at the golden tolerance,
/// which is the invocation the reference outputs were produced with. A warning
/// here would be a false positive on the one run we know is correct.
#[test]
fn stays_quiet_on_the_golden_case() {
    let err = run(&golden("misc/temp_byte_bip"), "10", "tolwarn_golden");
    assert!(
        !err.contains(WARNING),
        "false positive on the golden case:\n{err}"
    );
}
