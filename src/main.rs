use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;

use standseg::config::SegConfig;
use standseg::driver::{run_with_stage2, MemReport, Observer, Phase, Stage2Spec};
use standseg::image::Image;
use standseg::io;
use standseg::segment::{PassStats, Segmenter};
use standseg::stage2::{self, Stage2Config, Stage2Result};

/// Segment an image by region growing (Harward & Woodcock 1992).
#[derive(Parser, Debug)]
#[command(name = "standseg", version)]
struct Cli {
    /// Segmentation tolerances, comma separated
    #[arg(short = 't', value_delimiter = ',')]
    tols: Vec<f32>,

    /// Basename for output files
    #[arg(short = 'o')]
    base: String,

    /// Merge coefficient, 0 < cm <= 1
    #[arg(short = 'm', default_value_t = 1.0)]
    cm: f32,

    /// Nabsmin,Nnormin,Nviable,Nmax,Nabsmax
    #[arg(short = 'n', value_delimiter = ',')]
    n: Vec<u32>,

    /// Consider 8-way neighbours
    #[arg(short = '8', default_value_t = false)]
    eight: bool,

    /// Mask image; pixels valued 0 are excluded
    #[arg(short = 'M')]
    mask: Option<PathBuf>,

    /// Zero-based index of the band carrying the normality criterion. Requires -N.
    #[arg(short = 'B')]
    norm_band: Option<usize>,

    /// Normality interval low,high. A region whose -B band centroid falls
    /// outside it is "special" and is held to Nabsmin rather than Nnormin.
    #[arg(short = 'N', value_delimiter = ',', num_args = 1..=2, value_name = "LOW,HIGH")]
    norm_interval: Vec<f32>,

    /// Also write the auxiliary region map mask (<base>.armask.<pass>)
    #[arg(short = 'A', default_value_t = false)]
    armask: bool,

    /// Treat pixels with this value as nodata (water, cloud, non-treed area)
    #[arg(long, allow_negative_numbers = true)]
    nodata: Option<i64>,

    /// A pixel is nodata if ANY band matches, rather than all bands
    #[arg(long, default_value_t = false)]
    nodata_any: bool,

    /// Directory for output files
    #[arg(long, default_value = ".")]
    outdir: PathBuf,

    /// Output format for the region maps
    #[arg(long, value_enum, default_value_t = OutFormat::Envi)]
    format: OutFormat,

    /// Worker threads for the nearest-neighbour sweep. 0 = one per core,
    /// 1 = fully serial. Output is identical either way.
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Second-stage image: a *different* image over the same grid (forest
    /// structure, age, species). Replaces the auxiliary phase with segment
    /// development, Ye et al. 2025. Requires --n2.
    #[arg(long, value_name = "IMAGE")]
    stage2: Option<PathBuf>,

    /// Nmin,Nmax for stage 2: the minimum region size it merges up to, and the
    /// absolute ceiling a merge may not cross. Requires --stage2.
    #[arg(long = "n2", value_delimiter = ',', value_name = "NMIN,NMAX")]
    n2: Vec<u32>,

    /// Skip stage 1 and take its region map from this file, so stage 2 can be
    /// re-run against a different second image without re-segmenting. Requires
    /// --stage2, and then no input image or -t is needed.
    #[arg(long, value_name = "FILE")]
    rmap: Option<PathBuf>,

    /// Input image
    image: Option<PathBuf>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
enum OutFormat {
    /// Raw binary plus a .hdr sidecar -- byte-compatible with the originals.
    Envi,
    Tiff,
}

/// Prints the same running commentary as the C, so our log can be diffed
/// against the golden `myseg.log` line for line.
struct CLog {
    base: String,
    outdir: PathBuf,
    nbands: usize,
    masked: bool,
    geo: standseg::image::GeoRef,
    format: OutFormat,
    prov: io::Provenance,
    /// Set when -B/-N are in play: the C prints three extra lines then.
    norb: bool,
    /// Region-map width as of the end of stage 1. Stage 2 keeps stage-1 ids, so
    /// its own region count would understate how wide they are.
    stage2_nbytes: Option<usize>,
    /// Smallest tolerance and the image's observed sample range, kept only to
    /// warn when the two are wildly mismatched. `None` when a mask is in play,
    /// because then the region-to-pixel ratio no longer means what it should.
    tol_check: Option<(f32, f64, f64)>,
}

/// Fraction of pixels that may remain their own region after the first pass
/// before we say something. The real cases sit far below it: 88 % on Case 1,
/// 51 % on Case 2, 88 % on 16-bit Landsat at a correctly scaled tolerance.
/// An 8-bit tolerance on 16-bit data gives 99.0 - 100.0 %.
const LONELY_FRACTION: f64 = 0.98;

impl Observer for CLog {
    fn on_start(&mut self, nreg: usize, npix: usize) {
        println!("Completed the calculation of pixel nearest neighbors");
        println!("Initial pass over image completed");
        println!("{nreg} of a possible {npix} regions are required");
        if let Some((tol, lo, hi)) = self.tol_check {
            if npix > 0 && nreg as f64 / npix as f64 > LONELY_FRACTION {
                let pct = 100.0 * nreg as f64 / npix as f64;
                eprintln!(
                    "segment: warning: the first pass merged almost nothing -- {pct:.1} % of \
                     pixels are still their own region."
                );
                eprintln!(
                    "  -t is in DN units. This image runs {lo:.0} to {hi:.0}, and the smallest \
                     tolerance given is {tol}."
                );
                eprintln!(
                    "  The size rules will still force merges, so you will get a region map, but \
                     it will be shaped by region size rather than by spectral similarity. If this \
                     is 16-bit imagery, -t probably needs scaling up from an 8-bit value."
                );
            }
        }
        println!();
        println!("Creating region list");
        println!("Region list created");
        println!();
        println!("About to perform first general pass over region list");
        println!();
    }

    fn on_memory(&mut self, m: &MemReport) {
        let gb = |b: usize| b as f64 / 1e9;
        println!("Array sizes at peak:");
        println!(
            "\timage:        {:>8.3} GB (freed after this phase)",
            gb(m.image)
        );
        println!("\tcontiguity:   {:>8.3} GB", gb(m.cband));
        println!("\tregion band:  {:>8.3} GB", gb(m.rband));
        println!("\tregion list:  {:>8.3} GB", gb(m.rlist));
        println!("\tcentroids:    {:>8.3} GB", gb(m.ctrlist));
        println!("\tneighbours:   {:>8.3} GB", gb(m.nnbrlist));
        println!("\tpredicted peak: {:.3} GB", gb(m.peak()));
        println!();
    }

    fn on_pass(&mut self, phase: Phase, pass: usize, s: &PassStats) {
        match phase {
            Phase::Normal => {
                println!("Pass {pass} completed");
                println!(
                    "Tolerance for pass was {:.3}, (Tg = {:.3})",
                    s.tp2.sqrt(),
                    0.0
                );
                println!("{} regions remain after this pass", s.nreg);
                if s.no_nbr > 0 {
                    println!("{} regions possess no neighbors", s.no_nbr);
                }
                println!(
                    "The minimum nearest neighbor distance on this pass was {:.3}",
                    s.dmin2.sqrt()
                );
                println!(
                    "The largest region generated on this pass contained {} pixels",
                    s.maxpix
                );
                println!("Merges:\tattempted={}", s.merge_attempts);
                println!("\tnnbr_gone={}", s.nnbr_gone);
                println!("\twrong_partner={}", s.wrong_partner);
                println!("\tnnbr_d2_big={}", s.nnbr_d2_big);
                println!("\tboth_viable={}", s.both_viable);
                println!("\tnpix_big={}", s.npix_big);
                println!("\tmerging={}", s.merging);
            }
            Phase::Auxiliary => {
                println!("Auxiliary pass {pass} completed");
                println!("{} regions remain after this pass", s.nreg);
                if s.no_nbr > 0 {
                    println!("{} regions possess no neighbors", s.no_nbr);
                }
                if s.merge_attempts > 0 || s.special_merge_attempts > 0 {
                    println!(
                        "The minimum nearest neighbor distance on this pass was {:.3}",
                        s.dmin2.sqrt()
                    );
                }
                println!(
                    "The largest region generated on this pass contained {} pixels",
                    s.maxpix
                );
                println!(
                    "The smallest normal region remaining after this pass contained {} pixels",
                    s.norminpix
                );
                if self.norb {
                    println!(
                        "The smallest special region remaining after this pass contained {} pixels",
                        s.absminpix
                    );
                }
                println!("Normal merges:\tattempted={}", s.merge_attempts);
                if self.norb {
                    println!("Special merges:\tattempted={}", s.special_merge_attempts);
                }
                println!("\tnnbr_gone={}", s.nnbr_gone);
                println!("\twrong_partner={}", s.wrong_partner);
                println!("\tnpix_big={}", s.npix_big);
                println!("\tmerging={}", s.merging);
            }
            // Segment development keeps its own counters, which do not line up
            // with these -- it reports through `on_stage2` instead.
            Phase::Stage2 => return,
        }
        println!();
    }

    fn on_no_merges(&mut self, pass: usize) {
        println!("Pass {pass} resulted in no merges");
        println!();
    }

    fn on_stage2(&mut self, res: &Stage2Result, nbytes: usize) {
        self.stage2_nbytes = Some(nbytes);
        print_stage2(res);
    }

    fn on_map(&mut self, phase: Phase, pass: usize, seg: &Segmenter) -> Result<(), String> {
        let kind = match phase {
            Phase::Normal => "rmap",
            Phase::Auxiliary | Phase::Stage2 => "armap",
        };
        println!("Writing region map image");
        let nbytes = match phase {
            Phase::Stage2 => self
                .stage2_nbytes
                .unwrap_or_else(|| seg.region_map_nbytes()),
            _ => seg.region_map_nbytes(),
        };
        let path = match self.format {
            OutFormat::Envi => {
                let p = self.outdir.join(format!("{}.{kind}.{pass}", self.base));
                io::envi::write_region_map(
                    &p,
                    &seg.bands.rband,
                    seg.nlines,
                    seg.nsamps,
                    nbytes,
                    &self.geo,
                    self.masked,
                    &self.prov,
                )
                .map_err(|e| e.to_string())?;
                p
            }
            OutFormat::Tiff => {
                let p = self.outdir.join(format!("{}.{kind}.{pass}.tif", self.base));
                io::tiff::write_region_map(
                    &p,
                    &seg.bands.rband,
                    seg.nlines,
                    seg.nsamps,
                    nbytes,
                    &self.prov,
                )
                .map_err(|e| e.to_string())?;
                p
            }
        };
        let _ = &path;
        println!(
            "{}.{kind}.{pass} contains the final region map image",
            self.base
        );
        println!();

        // -A: the mask of which side of each Phase 2 merge was absorbed. The C
        // wrote it as <base>.armask.<pass>, one uint8 band.
        if phase == Phase::Auxiliary {
            if let Some(ab) = seg.aband.as_deref() {
                println!("Writing auxiliary region map mask");
                let p = self.outdir.join(format!("{}.armask.{pass}", self.base));
                let band: Vec<u32> = ab.iter().map(|&v| u32::from(v)).collect();
                io::envi::write_region_map(
                    &p, &band, seg.nlines, seg.nsamps, 1, &self.geo, false, &self.prov,
                )
                .map_err(|e| e.to_string())?;
                println!(
                    "{}.armask.{pass} contains the auxiliary region map mask",
                    self.base
                );
                println!();
            }
        }
        let _ = self.nbands;
        Ok(())
    }
}

/// The stage-2 equivalent of the C's per-pass commentary. Not a format anything
/// else parses -- the Python oracle prints nothing -- so it is shaped for
/// localising a divergence: the pass a merge count stops matching on is where to
/// look.
fn print_stage2(res: &Stage2Result) {
    println!("Segment development completed in {} passes", res.passes);
    if res.dropped_majority_nodata > 0 {
        println!(
            "{} regions were more than half nodata and were excluded",
            res.dropped_majority_nodata
        );
    }
    for (i, s) in res.stats.iter().enumerate() {
        println!("Development pass {} completed", i + 1);
        println!("{} regions remain after this pass", s.nreg);
        println!("Merges:\tattempted={}", s.considered);
        println!("\talready_merged={}", s.busy);
        println!("\tno_neighbor={}", s.no_cand);
        println!("\tno_distance={}", s.inf);
        println!("\tnpix_big={}", s.over_max);
        println!("\tnot_mutual={}", s.not_mutual);
        println!("\tmerging={}", s.merged);
    }
    println!();
}

/// Read the second image and check it against the grid stage 1 worked on.
fn read_stage2_image(path: &Path, nlines: usize, nsamps: usize) -> Result<Image, String> {
    let img = io::read(path).map_err(|e| e.to_string())?;
    if img.nlines != nlines || img.nsamps != nsamps {
        return Err(format!(
            "--stage2 image is {}x{} but the region map is {nlines}x{nsamps}; the \
             second stage reads a different image over the *same* grid",
            img.nlines, img.nsamps
        ));
    }
    println!(
        "Second-stage image has {} bands, {} samples, and {} lines ({} samples)",
        img.nbands,
        img.nsamps,
        img.nlines,
        img.data.kind()
    );
    Ok(img)
}

/// `--rmap`: stage 1 has already been run and its map is on disk, so run only
/// segment development. This is the path that lets a second image be swapped
/// without paying for the micro-segmentation again.
fn run_stage2_only(
    cli: &Cli,
    scfg: &Stage2Config,
    rmap_path: &Path,
    stage2_path: &Path,
) -> Result<(), String> {
    let rm = io::read_region_map(rmap_path).map_err(|e| e.to_string())?;
    println!(
        "Region map has {} samples and {} lines ({}-byte region ids)",
        rm.nsamps, rm.nlines, rm.nbytes
    );
    let img = read_stage2_image(stage2_path, rm.nlines, rm.nsamps)?;
    println!();

    let mut rband = rm.rband;
    let res = stage2::run(&mut rband, rm.nlines, rm.nsamps, &img, scfg)?;
    print_stage2(&res);

    std::fs::create_dir_all(&cli.outdir).map_err(|e| e.to_string())?;
    let prov = io::Provenance::from_args(std::env::args());
    println!("Writing region map image");
    match cli.format {
        OutFormat::Envi => {
            let p = cli
                .outdir
                .join(format!("{}.armap.{}", cli.base, res.passes));
            io::envi::write_region_map(
                &p, &rband, rm.nlines, rm.nsamps, rm.nbytes, &rm.geo, true, &prov,
            )
            .map_err(|e| e.to_string())?;
        }
        OutFormat::Tiff => {
            let p = cli
                .outdir
                .join(format!("{}.armap.{}.tif", cli.base, res.passes));
            io::tiff::write_region_map(&p, &rband, rm.nlines, rm.nsamps, rm.nbytes, &prov)
                .map_err(|e| e.to_string())?;
        }
    }
    println!(
        "{}.armap.{} contains the final region map image",
        cli.base, res.passes
    );
    Ok(())
}

/// Build a mask from an explicit mask image and/or a nodata value.
fn build_mask(
    img: &Image,
    mask_path: Option<&PathBuf>,
    nodata: Option<i64>,
    nodata_any: bool,
) -> Result<Option<Vec<u8>>, String> {
    let mut mask: Option<Vec<u8>> = None;

    if let Some(p) = mask_path {
        let m = io::read(p).map_err(|e| e.to_string())?;
        if m.nlines != img.nlines || m.nsamps != img.nsamps {
            return Err("input and mask images have different dimensions".into());
        }
        if m.nbands != 1 {
            return Err(format!("mask image has {} bands, expected 1", m.nbands));
        }
        // Any nonzero sample is valid data, whatever width the mask was stored at.
        mask = Some(m.to_mask());
    }

    if img.data.is_float() {
        return Err(format!(
            "the input image is {}; the first stage segments 8- and 16-bit \
             integer imagery only. A float image is accepted as the *second*-stage \
             layer -- pass it with --stage2, over an integer input or a --rmap.",
            img.data.kind()
        ));
    }

    if let Some(nd) = nodata {
        let (lo, hi) = img.data.value_range();
        if nd < lo || nd > hi {
            return Err(format!(
                "nodata value {nd} is outside the range of {} input ({lo} to {hi})",
                img.data.kind()
            ));
        }
        let mut v = mask.unwrap_or_else(|| vec![1u8; img.npixels()]);
        img.apply_nodata(nd, nodata_any, &mut v);
        mask = Some(v);
    }

    Ok(mask)
}

fn real_main() -> Result<(), String> {
    let cli = Cli::parse();

    // --- the two-input variant's own argument rules ----------------------
    // Deliberately strict. Every one of these is a case where the run would
    // otherwise silently do something other than what was asked.
    let scfg = match (cli.stage2.as_ref(), cli.n2.as_slice()) {
        (Some(_), [nmin, nmax]) => Some(Stage2Config {
            nmin: *nmin,
            nmax: *nmax,
        }),
        (Some(_), []) => return Err("--stage2 given but no size rules (--n2 Nmin,Nmax)".into()),
        (Some(_), _) => return Err("stage-2 size rules are --n2 Nmin,Nmax (two values)".into()),
        (None, []) => None,
        (None, _) => return Err("--n2 given but no second-stage image (--stage2)".into()),
    };
    if cli.rmap.is_some() && scfg.is_none() {
        return Err("--rmap skips stage 1, so there is nothing to do without --stage2".into());
    }
    if cli.rmap.is_some() && !cli.tols.is_empty() {
        return Err("--rmap skips stage 1, so -t has nothing to apply to".into());
    }
    if cli.rmap.is_some() && cli.image.is_some() {
        return Err("--rmap skips stage 1, so it takes no input image".into());
    }
    if cli.rmap.is_none() && cli.image.is_none() {
        return Err("no input image (or --rmap with --stage2)".into());
    }
    if cli.rmap.is_none() && cli.tols.is_empty() {
        return Err("at least one final tolerance (-t tol) required".into());
    }
    if scfg.is_some() && cli.armask {
        return Err(
            "-A records which side of each auxiliary-phase merge was absorbed, and              --stage2 replaces that phase; the two cannot be combined"
                .into(),
        );
    }
    if scfg.is_some() && cli.norm_band.is_some() {
        eprintln!(
            "segment: -B/-N only affect the auxiliary phase, which --stage2 replaces; \
             they will have no effect"
        );
    }

    // Stage 1 has nothing to contribute here: read the map and develop it.
    if let (Some(rp), Some(sp), Some(sc)) = (cli.rmap.as_ref(), cli.stage2.as_ref(), scfg.as_ref())
    {
        return run_stage2_only(&cli, sc, rp, sp);
    }
    let image = cli.image.clone().expect("checked above");

    if cli.cm <= 0.0 || cli.cm > 1.0 {
        return Err("merge coefficient must be > 0 and <= 1".into());
    }
    for t in &cli.tols {
        if *t < 0.0 {
            return Err("segmentation tolerance must be > 0".into());
        }
    }

    let cfg = SegConfig {
        tols: cli.tols.clone(),
        cm: cli.cm,
        ..Default::default()
    }
    .eight_way(cli.eight)
    .with_n(&cli.n)?;
    let cfg = SegConfig {
        threads: cli.threads,
        armm: cli.armask,
        ..cfg
    };
    // The C requires -B and -N together and refuses either alone.
    let cfg = match (cli.norm_band, cli.norm_interval.as_slice()) {
        (Some(b), [low, high]) => cfg.with_normality(b, *low, *high)?,
        (Some(_), []) => {
            return Err("normality band (-B) specified but no normality interval (-N)".into())
        }
        (None, [_, _]) => {
            return Err("normality interval (-N) but no normality band (-B) specified".into())
        }
        (None, []) => cfg,
        _ => return Err("normality interval is -N low,high (two values)".into()),
    };
    if cli.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .map_err(|e| e.to_string())?;
    }

    let (img, file_nodata) = io::read_with_nodata(&image).map_err(|e| e.to_string())?;
    println!(
        "Input image has {} bands, {} samples, and {} lines ({} samples)",
        img.nbands,
        img.nsamps,
        img.nlines,
        img.data.kind()
    );
    if cfg.tols.len() == 1 {
        println!("The segmentation tolerance is {:.6}", cfg.tols[0]);
    } else {
        println!("There are {} segmentation tolerances", cfg.tols.len());
    }
    println!("The merge coefficient is {:.6}", cfg.cm);
    println!();

    // An explicit --nodata wins; otherwise honour whatever the file declares
    // (ENVI `data ignore value`, GeoTIFF GDAL_NODATA).
    let (lo, hi) = img.data.value_range();
    let nodata = cli.nodata.or_else(|| {
        file_nodata.and_then(|v| {
            let r = v.round();
            if (v - r).abs() < 1e-9 && r >= lo as f64 && r <= hi as f64 {
                let r = r as i64;
                println!("Using nodata value {r} declared by the input file");
                Some(r)
            } else {
                eprintln!(
                    "segment: input declares nodata {v}, which is not a valid {} value; ignoring",
                    img.data.kind()
                );
                None
            }
        })
    });
    // The C caps the normality interval at 255 because pixels were uint8.
    // Check it against the input we actually have instead.
    if let Some(b) = cfg.norm_band {
        if b >= img.nbands {
            return Err(format!(
                "normality band (-B {b}) is out of range: the image has {} bands, \
                 indexed from 0",
                img.nbands
            ));
        }
        let (_, hi) = img.data.value_range();
        if cfg.nbhigh > hi as f32 {
            return Err(format!(
                "normality interval high ({}) exceeds the maximum {} value ({hi})",
                cfg.nbhigh,
                img.data.kind()
            ));
        }
        println!(
            "Special regions are those with band {b} outside ({:.6}, {:.6})",
            cfg.nblow, cfg.nbhigh
        );
    }

    let mask = build_mask(&img, cli.mask.as_ref(), nodata, cli.nodata_any)?;

    // Read the second image before segmenting, so a mismatched grid or an
    // unreadable file fails in a second rather than after the whole of stage 1.
    let stage2_img = match cli.stage2.as_ref() {
        Some(p) => Some(read_stage2_image(p, img.nlines, img.nsamps)?),
        None => None,
    };

    std::fs::create_dir_all(&cli.outdir).map_err(|e| e.to_string())?;
    let mut obs = CLog {
        base: cli.base.clone(),
        outdir: cli.outdir.clone(),
        nbands: img.nbands,
        // Segment development sets every pixel of a majority-non-treed region
        // to 0, so its output has nodata whether the input did or not.
        masked: mask.is_some() || stage2_img.is_some(),
        geo: img.geo.clone(),
        format: cli.format,
        // Recorded in the output header, the way IPW's `history` record was --
        // that line is how the invocation behind the golden fixtures was
        // recovered eleven years later.
        prov: io::Provenance::from_args(std::env::args()),
        norb: cfg.norm_band.is_some(),
        stage2_nbytes: None,
        tol_check: match (mask.is_some(), img.data.observed_range()) {
            (false, Some((lo, hi))) => cfg
                .tols
                .iter()
                .cloned()
                .fold(None::<f32>, |acc, t| Some(acc.map_or(t, |a: f32| a.min(t))))
                .map(|t| (t, lo, hi)),
            _ => None,
        },
    };

    let spec = match (stage2_img.as_ref(), scfg) {
        (Some(image), Some(cfg)) => Some(Stage2Spec { image, cfg }),
        _ => None,
    };
    let two_stage = spec.is_some();
    let r = run_with_stage2(img, &cfg, mask.as_deref(), &mut obs, spec)?;
    println!(
        "Normal segmentation completed in {} passes",
        r.normal_passes
    );
    if two_stage {
        println!("Segment development complete in {} passes", r.aux_passes);
    } else {
        println!("Auxiliary segmentation complete in {} passes", r.aux_passes);
    }
    Ok(())
}

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("segment: {e}");
            ExitCode::FAILURE
        }
    }
}
