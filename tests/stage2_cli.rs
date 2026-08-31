//! The CLI side of the two-input variant: that `--stage2` composes with stage 1
//! in one invocation, that it leaves stage 1 alone, and that the argument rules
//! actually refuse what they say they refuse.
//!
//! Output goes to `build/out/`, never into `tests/` -- the original program kept
//! inputs and outputs in one directory, which is how an oracle gets clobbered by
//! its own program.

use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_fast_segment");

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn outdir(name: &str) -> PathBuf {
    let d = root().join("build/out/stage2_cli").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(BIN)
        .current_dir(root())
        .args(args)
        .output()
        .expect("run fast_segment")
}

fn ok(args: &[&str]) -> std::process::Output {
    let o = run(args);
    assert!(
        o.status.success(),
        "fast_segment {args:?} failed:\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
    o
}

/// `--rmap` + `--stage2` reproduces the fixture through the command line, not
/// just through the library.
#[test]
fn the_cli_reproduces_a_fixture_byte_for_byte() {
    let d = outdir("rmap_only");
    ok(&[
        "--rmap",
        "tests/stage2/e2e_gsv/input/rmap",
        "--stage2",
        "tests/stage2/e2e_gsv/input/layer",
        "--n2",
        "50,8000",
        "-o",
        "gsv",
        "--outdir",
        d.to_str().unwrap(),
    ]);
    let ours = std::fs::read(d.join("gsv.armap.39")).expect("no gsv.armap.39 written");
    let theirs = std::fs::read(root().join("tests/stage2/e2e_gsv/expected/armap.39")).unwrap();
    assert_eq!(ours.len(), 80_000, "wrong payload size");
    assert_eq!(ours, theirs, "CLI output differs from the oracle");
}

/// One invocation must equal two. This is the composition the whole option
/// exists for: segment with `-t`, then develop the result against a second
/// image, without an intermediate file.
///
/// The second image here is the test scene itself. Feeding stage 2 the same
/// pixels stage 1 saw is not a meaningful segmentation, but it is a perfectly
/// good second image for testing that the two halves are wired together -- and
/// it needs no fixture that is not already in the repo.
#[test]
fn one_invocation_equals_two() {
    const SCENE: &str = "tests/golden/misc/temp_byte_bip";
    const TOLS: [&str; 6] = ["-t", "10", "-m", ".1", "-n", "15,15,100,2500,2500"];

    // (a) stage 1 alone, the golden invocation.
    let a = outdir("plain");
    let mut args: Vec<&str> = TOLS.to_vec();
    let ad = a.to_str().unwrap().to_string();
    args.extend(["-o", "a", "--outdir", &ad, SCENE]);
    ok(&args);

    // (b) the same run, with segment development bolted on.
    let b = outdir("composed");
    let mut args: Vec<&str> = TOLS.to_vec();
    let bd = b.to_str().unwrap().to_string();
    args.extend([
        "--stage2", SCENE, "--n2", "100,8000", "-o", "b", "--outdir", &bd, SCENE,
    ]);
    ok(&args);

    // Stage 1 is untouched by asking for stage 2 -- same pass count, same bytes.
    let a_rmap = std::fs::read(a.join("a.rmap.51")).expect("no a.rmap.51");
    let b_rmap = std::fs::read(b.join("b.rmap.51")).expect("--stage2 changed the stage-1 output");
    assert_eq!(a_rmap, b_rmap, "--stage2 perturbed the stage-1 region map");

    // (c) stage 2 on its own, fed the map stage 1 just wrote.
    let c = outdir("split");
    ok(&[
        "--rmap",
        a.join("a.rmap.51").to_str().unwrap(),
        "--stage2",
        SCENE,
        "--n2",
        "100,8000",
        "-o",
        "c",
        "--outdir",
        c.to_str().unwrap(),
    ]);

    let composed = std::fs::read_dir(&b)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .find(|n| n.starts_with("b.armap.") && !n.ends_with(".hdr"))
        .expect("no armap from the composed run");
    let split = composed.replace("b.armap.", "c.armap.");
    assert!(
        c.join(&split).exists(),
        "composed run converged as {composed}, split run did not produce {split}"
    );
    assert_eq!(
        std::fs::read(b.join(&composed)).unwrap(),
        std::fs::read(c.join(&split)).unwrap(),
        "one invocation and two produced different maps"
    );
}

/// The argument rules. Each of these would otherwise be a run that silently did
/// something other than what was asked.
#[test]
fn the_argument_rules_refuse_what_they_claim_to() {
    let d = outdir("refusals");
    let od = d.to_str().unwrap();
    let cases: [(&str, Vec<&str>); 6] = [
        (
            "--stage2 given but no size rules",
            vec![
                "-t",
                "10",
                "-o",
                "x",
                "--outdir",
                od,
                "--stage2",
                "tests/stage2/e2e_gsv/input/layer",
                "tests/golden/misc/temp_byte_bip",
            ],
        ),
        (
            "--n2 given but no second-stage image",
            vec![
                "-t",
                "10",
                "-o",
                "x",
                "--outdir",
                od,
                "--n2",
                "50,8000",
                "tests/golden/misc/temp_byte_bip",
            ],
        ),
        (
            "--rmap skips stage 1, so there is nothing to do without --stage2",
            vec![
                "-o",
                "x",
                "--outdir",
                od,
                "--rmap",
                "tests/stage2/e2e_gsv/input/rmap",
            ],
        ),
        (
            "--rmap skips stage 1, so -t has nothing to apply to",
            vec![
                "-t",
                "10",
                "-o",
                "x",
                "--outdir",
                od,
                "--rmap",
                "tests/stage2/e2e_gsv/input/rmap",
                "--stage2",
                "tests/stage2/e2e_gsv/input/layer",
                "--n2",
                "50,8000",
            ],
        ),
        (
            "no input image",
            vec!["-t", "10", "-o", "x", "--outdir", od],
        ),
        (
            "cannot be combined",
            vec![
                "-t",
                "10",
                "-o",
                "x",
                "--outdir",
                od,
                "-A",
                "--n2",
                "50,8000",
                "--stage2",
                "tests/golden/misc/temp_byte_bip",
                "tests/golden/misc/temp_byte_bip",
            ],
        ),
    ];
    for (want, args) in cases {
        let o = run(&args);
        assert!(!o.status.success(), "expected a refusal for {args:?}");
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(
            err.contains(want),
            "expected {want:?} in stderr, got:\n{err}"
        );
    }
}

/// A second image over a different grid is the one mistake that would otherwise
/// produce a plausible-looking wrong answer, so it must fail -- and fail before
/// stage 1 spends any time.
#[test]
fn a_mismatched_grid_is_refused() {
    let d = outdir("grid");
    let o = run(&[
        "-t",
        "10",
        "-o",
        "x",
        "--outdir",
        d.to_str().unwrap(),
        "--stage2",
        "tests/stage2/e2e_gsv/input/layer", // 200x200
        "--n2",
        "50,8000",
        "tests/golden/misc/temp_byte_bip", // 250x250
    ]);
    assert!(
        !o.status.success(),
        "a 200x200 second image was accepted for a 250x250 scene"
    );
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("same"), "unhelpful error: {err}");
}
