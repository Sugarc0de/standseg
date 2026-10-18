//! `standseg temporal`: segment a time series, carrying segments forward and
//! re-segmenting only what changed. See `standseg::temporal` for the algorithm.

use std::path::{Path, PathBuf};

use clap::Parser;

use standseg::config::SegConfig;
use standseg::driver::Observer;
use standseg::image::{GeoRef, Image};
use standseg::io;
use standseg::stage2::Stage2Config;
use standseg::temporal::{change_mask, normalize, Series};

use super::{print_stage2, write_gpkg, OutFormat};

/// Segment a time series of images over one grid. Step 0 is the two-stage run;
/// each later step re-segments only the pixels that changed since the step
/// before, pastes those segments into the previous result, and develops the
/// combined map against that step's second-stage layer.
///
/// Give the steps in time order, either as repeated flags (--image, --stage2,
/// and --classes or --change, each once per step) or as a --manifest table.
#[derive(Parser, Debug)]
#[command(name = "standseg temporal", version)]
pub struct TemporalCli {
    /// Stage-1 image for each step, in time order (repeat the flag)
    #[arg(long = "image", value_name = "IMAGE")]
    images: Vec<PathBuf>,

    /// Second-stage layer for each step, in time order (repeat the flag)
    #[arg(long = "stage2", value_name = "IMAGE")]
    layers: Vec<PathBuf>,

    /// Categorical layer for each step (repeat the flag). A pixel is
    /// re-segmented at step t when its class differs from step t-1.
    #[arg(long = "classes", value_name = "IMAGE", conflicts_with = "changes")]
    classes: Vec<PathBuf>,

    /// Change mask for each step *after the first* (repeat the flag); nonzero
    /// pixels are re-segmented. The alternative to --classes.
    #[arg(long = "change", value_name = "IMAGE")]
    changes: Vec<PathBuf>,

    /// CSV with a header row and columns `label,image,stage2` plus `classes` or
    /// `change` (empty on the first row). Relative paths are taken from the
    /// manifest's directory. Replaces the per-step flags.
    #[arg(long, value_name = "CSV")]
    manifest: Option<PathBuf>,

    /// Names for the steps, used in output filenames (default 0, 1, 2, ...)
    #[arg(long, value_delimiter = ',')]
    labels: Vec<String>,

    /// Rescale every image and second-stage layer to 1..255 per band, nodata
    /// (as the file declares it) to 0, before use. This is how the paper's
    /// inputs were prepared; without it the files are used as they are.
    #[arg(long, default_value_t = false)]
    normalize: bool,

    /// Segmentation tolerances, comma separated
    #[arg(short = 't', value_delimiter = ',', required = true)]
    tols: Vec<f32>,

    /// Merge coefficient, 0 < cm <= 1
    #[arg(short = 'm', default_value_t = 1.0)]
    cm: f32,

    /// Nabsmin,Nnormin,Nviable,Nmax,Nabsmax for stage 1
    #[arg(short = 'n', value_delimiter = ',')]
    n: Vec<u32>,

    /// Consider 8-way neighbours in stage 1
    #[arg(short = '8', default_value_t = false)]
    eight: bool,

    /// Nmin,Nmax for stage 2
    #[arg(long = "n2", value_delimiter = ',', value_name = "NMIN,NMAX", required = true)]
    n2: Vec<u32>,

    /// Basename for output files: <base>.<label>.armap.<passes> and friends
    #[arg(short = 'o')]
    base: String,

    /// Directory for output files
    #[arg(long, default_value = ".")]
    outdir: PathBuf,

    /// Output format for the final map of each step
    #[arg(long, value_enum, default_value_t = OutFormat::Envi)]
    format: OutFormat,

    /// Also write each step's stage-1 map, change mask and combined map (ENVI)
    #[arg(long, default_value_t = false)]
    intermediates: bool,

    /// Print stage 2's per-pass counters for every step
    #[arg(long, default_value_t = false)]
    verbose: bool,

    /// Worker threads for stage 1's nearest-neighbour sweep (0 = one per core)
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

enum ChangeSource {
    Classes(Vec<PathBuf>),
    Masks(Vec<PathBuf>),
}

struct Plan {
    labels: Vec<String>,
    images: Vec<PathBuf>,
    layers: Vec<PathBuf>,
    change: ChangeSource,
}

/// Stage 1 runs silently here: forty steps of the C's pass-by-pass commentary
/// would bury the one line per step that matters.
struct Quiet;
impl Observer for Quiet {}

fn read_manifest(path: &Path) -> Result<Plan, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("can't read {}: {e}", path.display()))?;
    let dir = path.parent().unwrap_or(Path::new("."));
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let header: Vec<String> = lines
        .next()
        .ok_or("manifest is empty")?
        .split(',')
        .map(|s| s.trim().to_ascii_lowercase())
        .collect();
    let col = |name: &str| header.iter().position(|h| h == name);
    let (ci, cs) = (
        col("image").ok_or("manifest needs an `image` column")?,
        col("stage2").ok_or("manifest needs a `stage2` column")?,
    );
    let cl = col("label");
    let (cc, is_classes) = match (col("classes"), col("change")) {
        (Some(c), None) => (Some(c), true),
        (None, Some(c)) => (Some(c), false),
        (None, None) => (None, true),
        (Some(_), Some(_)) => return Err("manifest has both `classes` and `change`".into()),
    };
    let resolve = |s: &str| {
        let p = PathBuf::from(s);
        if p.is_absolute() {
            p
        } else {
            dir.join(p)
        }
    };
    let mut plan = Plan {
        labels: vec![],
        images: vec![],
        layers: vec![],
        change: ChangeSource::Classes(vec![]),
    };
    let mut chg: Vec<PathBuf> = vec![];
    for (row, line) in lines.enumerate() {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if f.len() != header.len() {
            return Err(format!(
                "manifest row {} has {} fields, header has {}",
                row + 1,
                f.len(),
                header.len()
            ));
        }
        plan.labels
            .push(cl.map_or_else(|| row.to_string(), |c| f[c].to_string()));
        plan.images.push(resolve(f[ci]));
        plan.layers.push(resolve(f[cs]));
        if let Some(c) = cc {
            if is_classes || row > 0 {
                if f[c].is_empty() {
                    return Err(format!("manifest row {} has no {}", row + 1, header[c]));
                }
                chg.push(resolve(f[c]));
            }
        }
    }
    plan.change = if is_classes {
        ChangeSource::Classes(chg)
    } else {
        ChangeSource::Masks(chg)
    };
    Ok(plan)
}

fn plan_from(cli: &TemporalCli) -> Result<Plan, String> {
    let mut plan = match &cli.manifest {
        Some(m) => {
            if !cli.images.is_empty()
                || !cli.layers.is_empty()
                || !cli.classes.is_empty()
                || !cli.changes.is_empty()
            {
                return Err("--manifest replaces --image/--stage2/--classes/--change".into());
            }
            read_manifest(m)?
        }
        None => Plan {
            labels: (0..cli.images.len()).map(|i| i.to_string()).collect(),
            images: cli.images.clone(),
            layers: cli.layers.clone(),
            change: if cli.changes.is_empty() {
                ChangeSource::Classes(cli.classes.clone())
            } else {
                ChangeSource::Masks(cli.changes.clone())
            },
        },
    };
    let t = plan.images.len();
    if t == 0 {
        return Err("no steps: give --image (and --stage2) once per step, or --manifest".into());
    }
    if plan.layers.len() != t {
        return Err(format!(
            "{t} images but {} second-stage layers; give one per step",
            plan.layers.len()
        ));
    }
    match &plan.change {
        ChangeSource::Classes(c) if t > 1 && c.len() != t => {
            return Err(format!(
                "{t} steps but {} class layers; give one per step (or --change, one per step \
                 after the first)",
                c.len()
            ))
        }
        ChangeSource::Masks(c) if c.len() != t - 1 => {
            return Err(format!(
                "{t} steps need {} change masks (one per step after the first), got {}",
                t - 1,
                c.len()
            ))
        }
        _ => {}
    }
    if !cli.labels.is_empty() {
        if cli.labels.len() != t {
            return Err(format!("{t} steps but {} labels", cli.labels.len()));
        }
        plan.labels = cli.labels.clone();
    }
    let mut seen = std::collections::HashSet::new();
    for l in &plan.labels {
        if l.is_empty() || l.contains(['/', '\\']) || !seen.insert(l) {
            return Err(format!("step label {l:?} is empty, repeated, or contains a slash"));
        }
    }
    Ok(plan)
}

/// Read an image for one step, rescaled if asked. Returns the image and the
/// georeferencing to stamp on outputs.
fn load(path: &Path, rescale: bool) -> Result<Image, String> {
    let (img, nd) = io::read_with_nodata(path).map_err(|e| e.to_string())?;
    if !rescale {
        return Ok(img);
    }
    normalize(&img, nd).map_err(|e| format!("{}: {e}", path.display()))
}

fn write_envi(
    dir: &Path,
    name: &str,
    band: &[u32],
    nl: usize,
    ns: usize,
    nbytes: usize,
    geo: &GeoRef,
    masked: bool,
    prov: &io::Provenance,
) -> Result<(), String> {
    io::envi::write_region_map(&dir.join(name), band, nl, ns, nbytes, geo, masked, prov)
        .map_err(|e| e.to_string())
}

pub fn run(args: impl IntoIterator<Item = String>) -> Result<(), String> {
    let cli = TemporalCli::parse_from(args);
    let plan = plan_from(&cli)?;
    let s2 = match cli.n2.as_slice() {
        [nmin, nmax] => Stage2Config {
            nmin: *nmin,
            nmax: *nmax,
        },
        _ => return Err("stage-2 size rules are --n2 Nmin,Nmax (two values)".into()),
    };
    if cli.cm <= 0.0 || cli.cm > 1.0 {
        return Err("merge coefficient must be > 0 and <= 1".into());
    }
    if cli.tols.iter().any(|t| *t < 0.0) {
        return Err("segmentation tolerance must be > 0".into());
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
        ..cfg
    };
    if cli.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .map_err(|e| e.to_string())?;
    }
    std::fs::create_dir_all(&cli.outdir).map_err(|e| e.to_string())?;
    let prov = io::Provenance::from_args(std::env::args());

    println!(
        "{} steps; stage 1 -t {:?} -m {} ; stage 2 Nmin {} Nmax {}{}",
        plan.images.len(),
        cfg.tols,
        cfg.cm,
        s2.nmin,
        s2.nmax,
        if cli.normalize { " ; inputs rescaled to 1..255" } else { "" }
    );

    let mut series: Option<Series> = None;
    let mut prev_classes: Option<Image> = None;
    let mut geo = GeoRef::default();
    let (mut nl, mut ns) = (0, 0);
    for (t, label) in plan.labels.iter().enumerate() {
        let image = load(&plan.images[t], cli.normalize)?;
        if image.data.is_float() {
            return Err(format!(
                "{}: stage 1 segments 8- and 16-bit integer imagery only (try --normalize)",
                plan.images[t].display()
            ));
        }
        let layer = load(&plan.layers[t], cli.normalize)?;
        if t == 0 {
            (nl, ns) = (image.nlines, image.nsamps);
            geo = image.geo.clone();
        }

        let change: Option<Vec<u8>> = match &plan.change {
            ChangeSource::Classes(c) if !c.is_empty() => {
                let cur = io::read(&c[t]).map_err(|e| e.to_string())?;
                let m = match &prev_classes {
                    None => None,
                    Some(p) => Some(change_mask(p, &cur).map_err(|e| format!("step {label}: {e}"))?),
                };
                prev_classes = Some(cur);
                m
            }
            ChangeSource::Classes(_) => None,
            ChangeSource::Masks(c) => match t {
                0 => None,
                _ => {
                    let m = io::read(&c[t - 1]).map_err(|e| e.to_string())?;
                    if m.nbands != 1 {
                        return Err(format!("change mask {} has {} bands", c[t - 1].display(), m.nbands));
                    }
                    Some(m.to_mask())
                }
            },
        };

        let s = series.get_or_insert_with(|| Series::new(&cfg, s2.clone(), nl, ns));
        let r = s
            .step(image, change.as_deref(), &layer, &mut Quiet)
            .map_err(|e| format!("step {label}: {e}"))?;

        let stem = format!("{}.{label}", cli.base);
        if cli.intermediates {
            write_envi(&cli.outdir, &format!("{stem}.rmap.{}", r.stage1_passes), &r.stage1,
                       nl, ns, r.stage1_nbytes, &geo, t > 0, &prov)?;
            if let Some(m) = &change {
                let band: Vec<u32> = m.iter().map(|&v| u32::from(v)).collect();
                write_envi(&cli.outdir, &format!("{stem}.change"), &band, nl, ns, 1, &geo, false, &prov)?;
            }
            if let Some(c) = &r.combined {
                write_envi(&cli.outdir, &format!("{stem}.combined"), c, nl, ns, r.nbytes, &geo, true, &prov)?;
            }
        }
        let name = format!("{stem}.armap.{}", r.stage2.passes);
        match cli.format {
            OutFormat::Envi => write_envi(&cli.outdir, &name, &r.map, nl, ns, r.nbytes, &geo, true, &prov)?,
            OutFormat::Tiff => io::tiff::write_region_map(
                &cli.outdir.join(format!("{name}.tif")), &r.map, nl, ns, r.nbytes, &prov,
            )
            .map_err(|e| e.to_string())?,
            OutFormat::Gpkg => {
                write_gpkg(&cli.outdir.join(format!("{name}.gpkg")), &r.map, nl, ns, true, &geo,
                           &stem, &prov)?;
            }
        }
        println!(
            "step {label}: {} px re-segmented, stage 1 {} passes, stage 2 {} passes, {} segments -> {name}",
            r.changed, r.stage1_passes, r.stage2.passes, r.stage2.nreg
        );
        if cli.verbose {
            print_stage2(&r.stage2);
        }
    }
    Ok(())
}
