//! Gate on `tests/stage2` — the oracle for the two-input segment-development
//! phase (PLAN.md section 13, `tests/STAGE2.md`).
//!
//! The byte comparison lives next door in `stage2_match.rs`. This file does the
//! other half: it re-derives properties that any correct stage-2 output must
//! have, directly from the bytes on disk, using neither the Python that produced
//! them nor the Rust that has to reproduce them. So if someone regenerates the
//! fixtures with a broken oracle — and the Rust dutifully reproduces the broken
//! result — the byte test still passes and these do not.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

struct Case {
    name: &'static str,
    /// (lines, samples) and the number of stage-2 bands.
    shape: (usize, usize),
    nbands: usize,
    /// `-n` equivalents for stage 2: minimum and absolute maximum region size.
    nmin: u32,
    nmax: u32,
    expected: &'static str,
    /// ENVI data type of the region map: 12 = uint16, 13 = uint32.
    rmap_dtype: u32,
    regions_in: usize,
    regions_out: usize,
    /// Fraction of the *input* map that is already masked, x1000.
    masked_in_permille: usize,
}

/// Pinned from `tests/stage2/cases.json`. Kept as a literal table rather than
/// parsed so that a fixture regenerated with different numbers has to be
/// noticed and typed in here deliberately.
const CASES: &[Case] = &[
    Case {
        name: "tiny_synthetic",
        shape: (5, 5),
        nbands: 1,
        nmin: 4,
        nmax: 9,
        expected: "armap.4",
        rmap_dtype: 13,
        regions_in: 8,
        regions_out: 5,
        masked_in_permille: 0,
    },
    Case {
        name: "p95_250",
        shape: (250, 250),
        nbands: 1,
        nmin: 80,
        nmax: 8000,
        expected: "armap.71",
        rmap_dtype: 13,
        regions_in: 14712,
        regions_out: 272,
        masked_in_permille: 0,
    },
    Case {
        name: "species_250",
        shape: (250, 250),
        nbands: 6,
        nmin: 80,
        nmax: 8000,
        expected: "armap.78",
        rmap_dtype: 13,
        regions_in: 18154,
        regions_out: 382,
        masked_in_permille: 0,
    },
    Case {
        name: "age_capped",
        shape: (250, 250),
        nbands: 1,
        nmin: 60,
        nmax: 200,
        expected: "armap.40",
        rmap_dtype: 13,
        regions_in: 15091,
        regions_out: 1907,
        masked_in_permille: 0,
    },
    Case {
        name: "e2e_gsv",
        shape: (200, 200),
        nbands: 1,
        nmin: 50,
        nmax: 8000,
        expected: "armap.39",
        rmap_dtype: 12,
        regions_in: 6494,
        regions_out: 291,
        masked_in_permille: 0,
    },
    Case {
        name: "e2e_masked",
        shape: (200, 200),
        nbands: 1,
        nmin: 50,
        nmax: 8000,
        expected: "armap.39",
        rmap_dtype: 12,
        regions_in: 3171,
        regions_out: 468,
        masked_in_permille: 448,
    },
];

fn stage2(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/stage2")
        .join(rel)
}

struct Raster {
    nlines: usize,
    nsamps: usize,
    nbands: usize,
    data_type: u32,
    /// Band-sequential, widened to u32 whatever it was stored as.
    data: Vec<u32>,
}

/// Read one of our own ENVI fixtures. Deliberately a separate, minimal parser:
/// `io::envi` refuses data type 13, and using it here would make this test agree
/// with the reader rather than with the file.
fn read_envi(path: &Path) -> Raster {
    let hdr = std::fs::read_to_string(path.with_extension("hdr").with_file_name(format!(
        "{}.hdr",
        path.file_name().unwrap().to_string_lossy()
    )))
    .unwrap_or_else(|e| panic!("read {}.hdr: {e}", path.display()));

    let field = |k: &str| -> String {
        hdr.lines()
            .find(|l| l.trim_start().to_lowercase().starts_with(k))
            .and_then(|l| l.split('=').nth(1))
            .unwrap_or_else(|| panic!("{}: no '{k}' in header", path.display()))
            .trim()
            .to_string()
    };
    let num = |k: &str| field(k).parse::<usize>().unwrap();
    let (nlines, nsamps, nbands) = (num("lines"), num("samples"), num("bands"));
    let data_type = num("data type") as u32;
    assert_eq!(
        field("interleave"),
        "bsq",
        "{}: expected bsq",
        path.display()
    );
    assert_eq!(
        field("byte order"),
        "0",
        "{}: expected little-endian",
        path.display()
    );

    let raw = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let width = match data_type {
        1 => 1,
        12 => 2,
        13 => 4,
        d => panic!("{}: unexpected ENVI data type {d}", path.display()),
    };
    let n = nlines * nsamps * nbands;
    assert_eq!(
        raw.len(),
        n * width,
        "{}: wrong file length",
        path.display()
    );
    let data = raw
        .chunks_exact(width)
        .map(|c| {
            let mut b = [0u8; 4];
            b[..width].copy_from_slice(c);
            u32::from_le_bytes(b)
        })
        .collect();
    Raster {
        nlines,
        nsamps,
        nbands,
        data_type,
        data,
    }
}

fn sizes(map: &[u32]) -> HashMap<u32, usize> {
    let mut m = HashMap::new();
    for &v in map {
        if v != 0 {
            *m.entry(v).or_insert(0) += 1;
        }
    }
    m
}

/// The fixtures are shaped the way `case.json` says, in the container each
/// producer actually emitted — which is where the requirement that the region
/// map reader accept uint16 *and* uint32 comes from.
#[test]
fn every_case_reads_at_its_declared_shape() {
    for c in CASES {
        let rmap = read_envi(&stage2(&format!("{}/input/rmap", c.name)));
        let layer = read_envi(&stage2(&format!("{}/input/layer", c.name)));
        let exp = read_envi(&stage2(&format!("{}/expected/{}", c.name, c.expected)));

        assert_eq!(
            (rmap.nlines, rmap.nsamps),
            c.shape,
            "{}: rmap shape",
            c.name
        );
        assert_eq!(
            (layer.nlines, layer.nsamps),
            c.shape,
            "{}: layer shape",
            c.name
        );
        assert_eq!(
            (exp.nlines, exp.nsamps),
            c.shape,
            "{}: expected shape",
            c.name
        );
        assert_eq!(rmap.nbands, 1, "{}: region map must be single band", c.name);
        assert_eq!(exp.nbands, 1, "{}: region map must be single band", c.name);
        assert_eq!(layer.nbands, c.nbands, "{}: stage-2 band count", c.name);
        assert_eq!(rmap.data_type, c.rmap_dtype, "{}: rmap data type", c.name);
        assert_eq!(
            exp.data_type, c.rmap_dtype,
            "{}: expected data type",
            c.name
        );
        assert_eq!(layer.data_type, 1, "{}: stage-2 layers are uint8", c.name);

        assert_eq!(
            sizes(&rmap.data).len(),
            c.regions_in,
            "{}: regions in",
            c.name
        );
        assert_eq!(
            sizes(&exp.data).len(),
            c.regions_out,
            "{}: regions out",
            c.name
        );
        let masked = rmap.data.iter().filter(|&&v| v == 0).count() * 1000 / rmap.data.len();
        assert_eq!(
            masked, c.masked_in_permille,
            "{}: input mask fraction",
            c.name
        );
    }
}

/// Stage 2 only ever *merges*. Two pixels that shared a stage-1 region must
/// still share a stage-2 region — or have been masked out together. A fixture
/// that splits an input region did not come from this algorithm.
#[test]
fn no_case_splits_a_stage1_region() {
    for c in CASES {
        let rmap = read_envi(&stage2(&format!("{}/input/rmap", c.name)));
        let exp = read_envi(&stage2(&format!("{}/expected/{}", c.name, c.expected)));
        let mut seen: HashMap<u32, u32> = HashMap::new();
        for (&r, &e) in rmap.data.iter().zip(&exp.data) {
            if r == 0 {
                continue;
            }
            match seen.get(&r) {
                None => {
                    seen.insert(r, e);
                }
                Some(&first) => assert_eq!(
                    first, e,
                    "{}: stage-1 region {r} was split across stage-2 regions {first} and {e}",
                    c.name
                ),
            }
        }
    }
}

/// Merging cannot invent a region id (the absorbing region keeps its own), and
/// masking only ever grows: a pixel masked at stage 1 stays masked.
#[test]
fn output_ids_come_from_the_input_and_masking_only_grows() {
    for c in CASES {
        let rmap = read_envi(&stage2(&format!("{}/input/rmap", c.name)));
        let exp = read_envi(&stage2(&format!("{}/expected/{}", c.name, c.expected)));
        let have: HashSet<u32> = rmap.data.iter().copied().collect();
        for (&r, &e) in rmap.data.iter().zip(&exp.data) {
            if e != 0 {
                assert!(
                    have.contains(&e),
                    "{}: output invents region id {e}",
                    c.name
                );
            }
            if r == 0 {
                assert_eq!(
                    e, 0,
                    "{}: a masked stage-1 pixel came back as region {e}",
                    c.name
                );
            }
        }
    }
}

/// `Nmax` is checked on the *sum* before merging, so no region that actually
/// merged can exceed it. A region already over `Nmax` at stage 1 is left alone,
/// which is the one legal way to be larger.
#[test]
fn no_merged_region_exceeds_its_maximum() {
    for c in CASES {
        let rmap = read_envi(&stage2(&format!("{}/input/rmap", c.name)));
        let exp = read_envi(&stage2(&format!("{}/expected/{}", c.name, c.expected)));
        let before = sizes(&rmap.data);
        for (id, n) in sizes(&exp.data) {
            if n > c.nmax as usize {
                assert_eq!(
                    before.get(&id).copied(),
                    Some(n),
                    "{}: region {id} grew to {n} pixels, past Nmax {}",
                    c.name,
                    c.nmax
                );
            }
        }
    }
}

/// `age_capped` exists to make the maximum-region-size rejection bind, and
/// `species_250` to carry more than one stage-2 band. If either stops doing its
/// job the set has lost coverage even though every case still passes.
#[test]
fn the_set_still_covers_what_it_was_built_to_cover() {
    let capped = CASES.iter().find(|c| c.name == "age_capped").unwrap();
    let exp = read_envi(&stage2(&format!(
        "{}/expected/{}",
        capped.name, capped.expected
    )));
    let at_cap = sizes(&exp.data)
        .values()
        .filter(|&&n| n as u32 > capped.nmax / 2)
        .count();
    assert!(
        at_cap > 0,
        "age_capped no longer produces regions near Nmax; the cap is not being exercised"
    );

    assert!(
        CASES.iter().any(|c| c.nbands > 1),
        "no multi-band stage-2 case left in the set"
    );
    // Every case must actually reach its minimum somewhere -- a case where no
    // region ever gets to Nmin is not testing the phase, it is testing nothing.
    for c in CASES {
        assert!(c.nmin <= c.nmax, "{}: Nmin above Nmax", c.name);
        let exp = read_envi(&stage2(&format!("{}/expected/{}", c.name, c.expected)));
        assert!(
            sizes(&exp.data).values().any(|&n| n as u32 >= c.nmin),
            "{}: no output region reached Nmin {}",
            c.name,
            c.nmin
        );
    }
    assert!(
        CASES.iter().any(|c| c.masked_in_permille > 0),
        "no case left with a pre-masked stage-1 region map"
    );
    assert!(
        CASES.iter().any(|c| c.rmap_dtype == 12) && CASES.iter().any(|c| c.rmap_dtype == 13),
        "the set must exercise both uint16 and uint32 region maps"
    );
}
