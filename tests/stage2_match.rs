//! The byte gate on stage 2: run our segment-development phase over each
//! `tests/stage2` case and compare the region-map payload with what Elaine's
//! Python produced.
//!
//! This is what makes `tests/stage2/` a verifier rather than inert data. The
//! pass condition is the same as `tests/golden/`: byte-exact, and the pass count
//! in the expected filename matching the one we converge on. Region ids carry
//! meaning here -- the surviving id is the *absorbing* region's -- so "same
//! partition, different numbering" is a failure.

use std::path::{Path, PathBuf};

use standseg::io;
use standseg::stage2::{self, Stage2Config};

struct Case {
    name: &'static str,
    nmin: u32,
    nmax: u32,
    /// The expected filename encodes the pass count the oracle stopped at.
    expected: &'static str,
    regions_out: usize,
}

const CASES: &[Case] = &[
    Case {
        name: "tiny_synthetic",
        nmin: 4,
        nmax: 9,
        expected: "armap.4",
        regions_out: 5,
    },
    Case {
        name: "p95_250",
        nmin: 80,
        nmax: 8000,
        expected: "armap.71",
        regions_out: 272,
    },
    Case {
        name: "species_250",
        nmin: 80,
        nmax: 8000,
        expected: "armap.78",
        regions_out: 382,
    },
    Case {
        name: "age_capped",
        nmin: 60,
        nmax: 200,
        expected: "armap.40",
        regions_out: 1907,
    },
    Case {
        name: "e2e_gsv",
        nmin: 50,
        nmax: 8000,
        expected: "armap.39",
        regions_out: 291,
    },
    Case {
        name: "e2e_masked",
        nmin: 50,
        nmax: 8000,
        expected: "armap.39",
        regions_out: 468,
    },
];

fn stage2_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/stage2")
        .join(rel)
}

/// Serialise a region band the way the fixture stores it: little-endian, one
/// band, `nbytes` per sample.
fn pack(rband: &[u32], nbytes: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(rband.len() * nbytes);
    for &v in rband {
        out.extend_from_slice(&v.to_le_bytes()[..nbytes]);
    }
    out
}

/// Where the two maps first differ, as `(index, ours, theirs)`.
fn first_diff(a: &[u8], b: &[u8]) -> Option<(usize, u8, u8)> {
    a.iter()
        .zip(b)
        .enumerate()
        .find(|(_, (x, y))| x != y)
        .map(|(i, (x, y))| (i, *x, *y))
}

fn run_case(c: &Case) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let rm = io::read_region_map(&stage2_path(&format!("{}/input/rmap", c.name)))
        .unwrap_or_else(|e| panic!("{}: read rmap: {e}", c.name));
    let layer = io::read(&stage2_path(&format!("{}/input/layer", c.name)))
        .unwrap_or_else(|e| panic!("{}: read layer: {e}", c.name));

    let mut rband = rm.rband.clone();
    let res = stage2::run(
        &mut rband,
        rm.nlines,
        rm.nsamps,
        &layer,
        &Stage2Config {
            nmin: c.nmin,
            nmax: c.nmax,
        },
    )
    .unwrap_or_else(|e| panic!("{}: stage 2: {e}", c.name));

    let expected = std::fs::read(stage2_path(&format!("{}/expected/{}", c.name, c.expected)))
        .unwrap_or_else(|e| panic!("{}: read expected: {e}", c.name));

    (pack(&rband, rm.nbytes), expected, res.passes, res.nreg)
}

/// The whole point. Every case, byte for byte.
#[test]
fn every_case_reproduces_the_oracle_byte_for_byte() {
    for c in CASES {
        let (ours, theirs, passes, nreg) = run_case(c);
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{}: produced {} bytes, expected {}",
            c.name,
            ours.len(),
            theirs.len()
        );
        if let Some((i, a, b)) = first_diff(&ours, &theirs) {
            panic!(
                "{}: region map differs from the oracle at byte {i} (ours {a:#04x}, \
                 oracle {b:#04x}); {} of {} bytes differ",
                c.name,
                ours.iter().zip(&theirs).filter(|(x, y)| x != y).count(),
                ours.len()
            );
        }

        // The pass count is part of the answer: the oracle names its output
        // after the pass it stopped on, so a map that matches after a different
        // number of passes means the convergence rule diverged.
        let want: usize = c.expected.rsplit('.').next().unwrap().parse().unwrap();
        assert_eq!(
            passes, want,
            "{}: converged in {passes} passes, oracle {want}",
            c.name
        );
        assert_eq!(nreg, c.regions_out, "{}: region count", c.name);
    }
}

/// A negative control. If the comparison above passed because both sides were
/// empty, or because `run` quietly returned its input, this would pass too --
/// so perturb one parameter and require the bytes to move. `tiny_synthetic` is
/// 25 pixels and its `Nmax` of 9 is what stops two of its merges.
#[test]
fn the_comparison_can_fail() {
    let c = &CASES[0];
    let rm = io::read_region_map(&stage2_path(&format!("{}/input/rmap", c.name))).unwrap();
    let layer = io::read(&stage2_path(&format!("{}/input/layer", c.name))).unwrap();
    let expected =
        std::fs::read(stage2_path(&format!("{}/expected/{}", c.name, c.expected))).unwrap();

    let mut moved = 0;
    for (nmin, nmax) in [(c.nmin, c.nmax + 90), (c.nmin + 2, c.nmax)] {
        let mut rband = rm.rband.clone();
        stage2::run(
            &mut rband,
            rm.nlines,
            rm.nsamps,
            &layer,
            &Stage2Config { nmin, nmax },
        )
        .unwrap();
        if pack(&rband, rm.nbytes) != expected {
            moved += 1;
        }
    }
    assert_eq!(
        moved, 2,
        "changing Nmin/Nmax left the output identical; the gate is inert"
    );
}

/// Stage 2 must be a pure function of its inputs -- no RNG, no map iteration
/// order leaking into the result. Under the pinned tie-break rule it consumes no
/// randomness at all, unlike stage 1.
#[test]
fn the_phase_is_deterministic() {
    for c in CASES
        .iter()
        .filter(|c| c.name != "p95_250" && c.name != "species_250")
    {
        let (a, _, pa, _) = run_case(c);
        let (b, _, pb, _) = run_case(c);
        assert_eq!(a, b, "{}: two runs disagree", c.name);
        assert_eq!(pa, pb, "{}: two runs took different pass counts", c.name);
    }
}

/// Every pass, not just the last one.
///
/// Byte equality of the final map is the pass condition, but it says nothing
/// about *where* a future change went wrong. `case.json` records the oracle's
/// per-pass merge and rejection counts, so a divergence names the pass and the
/// reason instead of a byte offset.
#[test]
fn the_per_pass_counters_match_the_oracle() {
    for c in CASES {
        let json = std::fs::read_to_string(stage2_path(&format!("{}/case.json", c.name)))
            .unwrap_or_else(|e| panic!("{}: read case.json: {e}", c.name));
        let want = parse_per_pass(&json);

        let rm = io::read_region_map(&stage2_path(&format!("{}/input/rmap", c.name))).unwrap();
        let layer = io::read(&stage2_path(&format!("{}/input/layer", c.name))).unwrap();
        let mut rband = rm.rband.clone();
        let got = stage2::run(
            &mut rband,
            rm.nlines,
            rm.nsamps,
            &layer,
            &Stage2Config {
                nmin: c.nmin,
                nmax: c.nmax,
            },
        )
        .unwrap();

        assert_eq!(got.stats.len(), want.len(), "{}: number of passes", c.name);
        for (i, (g, w)) in got.stats.iter().zip(&want).enumerate() {
            let mine = [
                g.considered,
                g.no_cand,
                g.busy,
                g.inf,
                g.over_max,
                g.not_mutual,
                g.merged,
            ];
            assert_eq!(
                mine,
                *w,
                "{}: pass {} counters differ \
                 [considered,no_cand,busy,inf,over_max,not_mutual,merged]",
                c.name,
                i + 1
            );
        }
        // The loop stops on a pass that merges nothing -- that is the whole
        // termination rule, so make it explicit rather than implied.
        assert_eq!(
            got.stats.last().unwrap().merged,
            0,
            "{}: last pass merged",
            c.name
        );
    }
}

/// Pull the `per_pass` array out of a `case.json` as
/// `[considered, no_cand, busy, inf, over_max, not_mutual, merged]` per entry.
///
/// A hand-rolled scanner rather than a JSON crate: this is a test fixture we
/// wrote ourselves in a shape we control, and a dependency here would be a
/// dependency in the shipped binary's tree.
fn parse_per_pass(json: &str) -> Vec<[usize; 7]> {
    const KEYS: [&str; 7] = [
        "considered",
        "no_cand",
        "busy",
        "inf",
        "over_max",
        "not_mutual",
        "merged",
    ];
    let start = json
        .find("\"per_pass\"")
        .expect("case.json has no per_pass");
    let body = &json[start..];
    let end = body.find(']').expect("per_pass array is unterminated");
    let mut out = Vec::new();
    for obj in body[..end].split('{').skip(1) {
        let obj = &obj[..obj.find('}').expect("per_pass entry is unterminated")];
        let mut row = [usize::MAX; 7];
        for (k, slot) in KEYS.iter().zip(&mut row) {
            let pat = format!("\"{k}\":");
            let at = obj
                .find(&pat)
                .unwrap_or_else(|| panic!("per_pass entry has no {k}: {obj}"))
                + pat.len();
            let digits: String = obj[at..]
                .trim_start()
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect();
            *slot = digits.parse().unwrap_or_else(|e| panic!("{k}: {e}"));
        }
        out.push(row);
    }
    assert!(!out.is_empty(), "per_pass was empty");
    out
}
