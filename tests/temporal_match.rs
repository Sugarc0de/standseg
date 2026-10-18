//! The byte gate on `standseg temporal`: run each `tests/temporal` case step by
//! step and compare every artifact the oracle kept -- the stage-1 map, the
//! change mask, the combined map and the final map -- with what the reference
//! pipeline (`tools/temporal_oracle`) wrote. Byte-exact, pass counts included:
//! the expected file names carry them.
//!
//! Output from the CLI tests goes to `build/out/temporal_test/`, never into
//! `tests/`.

use std::path::{Path, PathBuf};
use std::process::Command;

use standseg::config::SegConfig;
use standseg::driver::Observer;
use standseg::image::Image;
use standseg::io;
use standseg::stage2::Stage2Config;
use standseg::temporal::{change_mask, normalize, Series};

const BIN: &str = env!("CARGO_BIN_EXE_standseg");

struct Case {
    name: &'static str,
    t: &'static str,
    m: &'static str,
    n: &'static str,
    n2: (u32, u32),
    labels: &'static [&'static str],
}

const PAPER_N: &str = "18,18,18,36,36";
const SEVEN_01: &[&str] = &["2001", "2002", "2003", "2004", "2005", "2006", "2007"];

const CASES: &[Case] = &[
    Case {
        name: "tiny_synthetic",
        t: "50",
        m: "1.0",
        n: "16,16,32,200,200",
        n2: (6, 200),
        labels: &["1", "2", "3"],
    },
    Case {
        name: "fire1_crop",
        t: "50",
        m: "0.2",
        n: PAPER_N,
        n2: (80, 8000),
        labels: SEVEN_01,
    },
    Case {
        name: "fire2_crop",
        t: "50",
        m: "0.2",
        n: PAPER_N,
        n2: (80, 8000),
        labels: &["1996", "1997", "1998", "1999", "2000", "2001", "2002"],
    },
    Case {
        name: "fire3_crop",
        t: "50",
        m: "0.2",
        n: PAPER_N,
        n2: (80, 8000),
        labels: SEVEN_01,
    },
];

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn case_dir(c: &Case) -> PathBuf {
    root().join("tests/temporal").join(c.name)
}

fn outdir(name: &str) -> PathBuf {
    let d = root().join("build/out/temporal_test").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

struct Quiet;
impl Observer for Quiet {}

fn pack(band: &[u32], nbytes: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(band.len() * nbytes);
    for &v in band {
        out.extend_from_slice(&v.to_le_bytes()[..nbytes]);
    }
    out
}

/// The single expected file for `<label>.<kind>[.<passes>]`, with its pass count.
fn expected(c: &Case, label: &str, kind: &str) -> Option<(PathBuf, Option<usize>)> {
    let dir = case_dir(c).join("expected");
    let prefix = format!("{label}.{kind}");
    let mut hits: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| !n.ends_with(".hdr") && (*n == prefix || n.starts_with(&format!("{prefix}."))))
        .collect();
    assert!(
        hits.len() <= 1,
        "{}: several expected {prefix}: {hits:?}",
        c.name
    );
    let n = hits.pop()?;
    let passes = n
        .strip_prefix(&format!("{prefix}."))
        .map(|p| p.parse().unwrap());
    Some((dir.join(n), passes))
}

fn load(p: &Path) -> Image {
    let (img, nd) = io::read_with_nodata(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    normalize(&img, nd).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn input(c: &Case, label: &str, what: &str) -> PathBuf {
    case_dir(c).join("input").join(label).join(what)
}

fn cfg(c: &Case) -> SegConfig {
    SegConfig {
        tols: c.t.split(',').map(|v| v.parse().unwrap()).collect(),
        cm: c.m.parse().unwrap(),
        ..Default::default()
    }
    .with_n(
        &c.n.split(',')
            .map(|v| v.parse().unwrap())
            .collect::<Vec<u32>>(),
    )
    .unwrap()
}

fn assert_bytes(c: &Case, what: &str, ours: &[u8], theirs: &[u8]) {
    assert_eq!(
        ours.len(),
        theirs.len(),
        "{} {what}: {} bytes, oracle {}",
        c.name,
        ours.len(),
        theirs.len()
    );
    if let Some(i) = ours.iter().zip(theirs).position(|(a, b)| a != b) {
        panic!(
            "{} {what}: differs from the oracle at byte {i}; {} of {} bytes differ",
            c.name,
            ours.iter().zip(theirs).filter(|(a, b)| a != b).count(),
            ours.len()
        );
    }
}

/// Run a case through the library, optionally with a different stage-2 rule.
/// Returns the final map's bytes per step.
fn run_library(c: &Case, n2: (u32, u32), check: bool) -> Vec<Vec<u8>> {
    let cfg = cfg(c);
    let mut series: Option<Series> = None;
    let mut prev_classes: Option<Image> = None;
    let mut finals = vec![];
    for label in c.labels {
        let image = load(&input(c, label, "image"));
        let layer = load(&input(c, label, "stage2"));
        let classes = io::read(&input(c, label, "classes")).unwrap();
        let change = prev_classes
            .as_ref()
            .map(|p| change_mask(p, &classes).unwrap());
        prev_classes = Some(classes);
        let (nl, ns) = (image.nlines, image.nsamps);
        let s = series.get_or_insert_with(|| {
            Series::new(
                &cfg,
                Stage2Config {
                    nmin: n2.0,
                    nmax: n2.1,
                },
                nl,
                ns,
            )
        });
        let r = s
            .step(image, change.as_deref(), &layer, &mut Quiet)
            .unwrap();
        let fin = pack(&r.map, r.nbytes);
        if check {
            let what = |k: &str| format!("step {label} {k}");
            let (p, n) = expected(c, label, "rmap").expect("expected rmap");
            assert_bytes(
                c,
                &what("stage-1 map"),
                &pack(&r.stage1, r.stage1_nbytes),
                &std::fs::read(p).unwrap(),
            );
            assert_eq!(
                n,
                Some(r.stage1_passes),
                "{} {}: stage-1 passes",
                c.name,
                what("rmap")
            );

            match (&change, expected(c, label, "change")) {
                (None, None) => {}
                (Some(m), Some((p, _))) => {
                    assert_bytes(c, &what("change mask"), m, &std::fs::read(p).unwrap())
                }
                (m, e) => panic!(
                    "{} {}: change mask ours {} oracle {}",
                    c.name,
                    label,
                    m.is_some(),
                    e.is_some()
                ),
            }
            match (&r.combined, expected(c, label, "combined")) {
                (None, None) => {}
                (Some(m), Some((p, _))) => assert_bytes(
                    c,
                    &what("combined map"),
                    &pack(m, r.nbytes),
                    &std::fs::read(p).unwrap(),
                ),
                (m, e) => panic!(
                    "{} {}: combined ours {} oracle {}",
                    c.name,
                    label,
                    m.is_some(),
                    e.is_some()
                ),
            }
            let (p, n) = expected(c, label, "armap").expect("expected armap");
            assert_bytes(c, &what("final map"), &fin, &std::fs::read(p).unwrap());
            assert_eq!(
                n,
                Some(r.stage2.passes),
                "{} {}: stage-2 passes",
                c.name,
                what("armap")
            );
        }
        finals.push(fin);
    }
    finals
}

/// The whole point: every artifact of every step of every case, byte for byte.
#[test]
fn every_step_of_every_case_reproduces_the_oracle() {
    for c in CASES {
        run_library(c, c.n2, true);
    }
}

fn cli(c: &Case, d: &Path, extra: &[&str]) {
    let n2 = format!("{},{}", c.n2.0, c.n2.1);
    let mut args = vec![
        "temporal",
        "--normalize",
        "--intermediates",
        "-t",
        c.t,
        "-m",
        c.m,
        "-n",
        c.n,
        "--n2",
        &n2,
        "-o",
        "case",
        "--outdir",
        d.to_str().unwrap(),
    ];
    args.extend(extra);
    let o = Command::new(BIN)
        .current_dir(root())
        .args(&args)
        .output()
        .expect("run standseg");
    assert!(
        o.status.success(),
        "standseg {args:?} failed:\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

/// Every expected file, as the command line writes it: `<base>.<expected name>`.
fn compare_dir(c: &Case, d: &Path) -> usize {
    let mut n = 0;
    for e in std::fs::read_dir(case_dir(c).join("expected")).unwrap() {
        let name = e.unwrap().file_name().to_string_lossy().into_owned();
        if name.ends_with(".hdr") {
            continue;
        }
        let ours = std::fs::read(d.join(format!("case.{name}")))
            .unwrap_or_else(|_| panic!("{}: CLI wrote no case.{name}", c.name));
        let theirs = std::fs::read(case_dir(c).join("expected").join(&name)).unwrap();
        assert_bytes(c, &format!("CLI {name}"), &ours, &theirs);
        n += 1;
    }
    n
}

/// The same through the command line, driven by each case's manifest.
#[test]
fn the_cli_reproduces_every_case_from_its_manifest() {
    let mut files = 0;
    for c in CASES {
        let d = outdir(&format!("manifest_{}", c.name));
        let m = case_dir(c).join("input/manifest.csv");
        cli(c, &d, &["--manifest", m.to_str().unwrap()]);
        files += compare_dir(c, &d);
    }
    assert_eq!(files, 88, "expected-file count changed");
}

/// Change masks given directly (`--change`) must give what the class layers
/// they were derived from give (`--classes`). The oracle's own masks are used.
#[test]
fn change_masks_and_class_layers_agree() {
    for c in [&CASES[0], &CASES[1]] {
        let d = outdir(&format!("change_{}", c.name));
        let mut args: Vec<String> = vec![];
        for (t, label) in c.labels.iter().enumerate() {
            args.push("--image".into());
            args.push(input(c, label, "image").to_string_lossy().into());
            args.push("--stage2".into());
            args.push(input(c, label, "stage2").to_string_lossy().into());
            if t > 0 {
                args.push("--change".into());
                args.push(
                    expected(c, label, "change")
                        .unwrap()
                        .0
                        .to_string_lossy()
                        .into(),
                );
            }
        }
        args.push("--labels".into());
        args.push(c.labels.join(","));
        let a: Vec<&str> = args.iter().map(String::as_str).collect();
        cli(c, &d, &a);
        assert!(compare_dir(c, &d) > 0);
    }
}

/// A negative control: if the gate passed because nothing was compared, or
/// because the steps ignored their parameters, this would pass too. Raising
/// stage 2's Nmin on the synthetic case to 60 puts its 56-pixel right-hand
/// segment under the floor, so the final map of every step must move.
#[test]
fn the_comparison_can_fail() {
    let c = &CASES[0];
    let base = run_library(c, c.n2, false);
    let moved = run_library(c, (60, c.n2.1), false);
    for (i, (a, b)) in base.iter().zip(&moved).enumerate() {
        assert_ne!(
            a, b,
            "step {i}: Nmin 60 left the final map identical; the gate is inert"
        );
    }
}

/// A step where nothing changed re-segments nothing: the result is the
/// previous step's segments, developed against the new layer. With an
/// unchanged layer too, the map must come back as it was.
#[test]
fn a_step_with_no_change_carries_everything_forward() {
    let c = &CASES[0];
    let cfg = cfg(c);
    let image = load(&input(c, "1", "image"));
    let layer = load(&input(c, "1", "stage2"));
    let (nl, ns) = (image.nlines, image.nsamps);
    let mut s = Series::new(
        &cfg,
        Stage2Config {
            nmin: c.n2.0,
            nmax: c.n2.1,
        },
        nl,
        ns,
    );
    let r0 = s.step(image.clone(), None, &layer, &mut Quiet).unwrap();
    let none = vec![0u8; nl * ns];
    let r1 = s.step(image, Some(&none), &layer, &mut Quiet).unwrap();
    assert_eq!(r1.changed, 0);
    assert!(
        r1.stage1.iter().all(|&v| v == 0),
        "stage 1 labelled a masked pixel"
    );
    assert_eq!(r1.map, r0.map, "an unchanged step altered the map");
}
