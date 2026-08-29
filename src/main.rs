use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;

use fast_segment::config::SegConfig;
use fast_segment::driver::{run, MemReport, Observer, Phase};
use fast_segment::image::Image;
use fast_segment::io;
use fast_segment::segment::{PassStats, Segmenter};

/// Segment an image by region growing (Harward & Woodcock 1992).
#[derive(Parser, Debug)]
#[command(name = "segment", version)]
struct Cli {
    /// Segmentation tolerances, comma separated
    #[arg(short = 't', value_delimiter = ',', required = true)]
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

    /// Treat pixels with this value as nodata (water, non-treed area)
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

    /// Input image
    image: PathBuf,
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
    geo: fast_segment::image::GeoRef,
    format: OutFormat,
}

impl Observer for CLog {
    fn on_start(&mut self, nreg: usize, npix: usize) {
        println!("Completed the calculation of pixel nearest neighbors");
        println!("Initial pass over image completed");
        println!("{nreg} of a possible {npix} regions are required");
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
        println!("\timage:        {:>8.3} GB (freed after this phase)", gb(m.image));
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
                println!("Tolerance for pass was {:.3}, (Tg = {:.3})", s.tp2.sqrt(), 0.0);
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
                println!("Normal merges:\tattempted={}", s.merge_attempts);
                println!("\tnnbr_gone={}", s.nnbr_gone);
                println!("\twrong_partner={}", s.wrong_partner);
                println!("\tnpix_big={}", s.npix_big);
                println!("\tmerging={}", s.merging);
            }
        }
        println!();
    }

    fn on_no_merges(&mut self, pass: usize) {
        println!("Pass {pass} resulted in no merges");
        println!();
    }

    fn on_map(&mut self, phase: Phase, pass: usize, seg: &Segmenter) -> Result<(), String> {
        let kind = match phase {
            Phase::Normal => "rmap",
            Phase::Auxiliary => "armap",
        };
        println!("Writing region map image");
        let nbytes = seg.region_map_nbytes();
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
                )
                .map_err(|e| e.to_string())?;
                p
            }
            OutFormat::Tiff => {
                let p = self.outdir.join(format!("{}.{kind}.{pass}.tif", self.base));
                io::tiff::write_region_map(&p, &seg.bands.rband, seg.nlines, seg.nsamps, nbytes)
                    .map_err(|e| e.to_string())?;
                p
            }
        };
        let _ = &path;
        println!("{}.{kind}.{pass} contains the final region map image", self.base);
        println!();
        let _ = self.nbands;
        Ok(())
    }
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
    let cfg = SegConfig { threads: cli.threads, ..cfg };
    if cli.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .map_err(|e| e.to_string())?;
    }

    let (img, file_nodata) = io::read_with_nodata(&cli.image).map_err(|e| e.to_string())?;
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
    let mask = build_mask(&img, cli.mask.as_ref(), nodata, cli.nodata_any)?;

    std::fs::create_dir_all(&cli.outdir).map_err(|e| e.to_string())?;
    let mut obs = CLog {
        base: cli.base.clone(),
        outdir: cli.outdir.clone(),
        nbands: img.nbands,
        masked: mask.is_some(),
        geo: img.geo.clone(),
        format: cli.format,
    };

    let r = run(img, &cfg, mask.as_deref(), &mut obs)?;
    println!(
        "Normal segmentation completed in {} passes",
        r.normal_passes
    );
    println!(
        "Auxiliary segmentation complete in {} passes",
        r.aux_passes
    );
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
