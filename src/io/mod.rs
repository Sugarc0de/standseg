//! Image input and output.
//!
//! Formats: ENVI (raw + `.hdr`) and IPW natively; TIFF and PNG to follow.
//! No libgdal -- the original's GDAL dependency covered exactly the two formats
//! we can write by hand. See PLAN.md section 9.

pub mod envi;
pub mod gpkg;
pub mod ipw;
pub mod png;
pub mod tiff;

use std::fmt;
use std::path::Path;

use crate::image::Image;

#[derive(Debug)]
pub struct IoError(String);

impl IoError {
    pub fn new(msg: impl Into<String>) -> Self {
        IoError(msg.into())
    }
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for IoError {}

pub type Result<T> = std::result::Result<T, IoError>;

/// What produced an output file.
///
/// Modelled on IPW's `history` record -- `segment -t 10 -m .1 -n ...` -- which
/// is how the command that made the golden fixtures was recovered eleven years
/// later. Our ENVI output had no equivalent, which was a regression against
/// 1992. Deliberately deterministic: no timestamp, so running the same command
/// twice produces identical files.
#[derive(Debug, Clone, Default)]
pub struct Provenance {
    /// The command line, reassembled and quoted so it can be pasted back.
    pub command: String,
    /// Program name and version.
    pub software: String,
}

impl Provenance {
    /// Build from an argument iterator, e.g. `std::env::args()`.
    pub fn from_args<I: IntoIterator<Item = String>>(args: I) -> Self {
        let command = args
            .into_iter()
            .map(|a| quote_arg(&a))
            .collect::<Vec<_>>()
            .join(" ");
        Self {
            command,
            software: format!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION")),
        }
    }
}

/// Make one argument safe to embed in a header and to paste back into a shell.
///
/// `{` and `}` would terminate an ENVI block early, and a newline would end the
/// record, so both are replaced rather than escaped -- a mangled history line is
/// better than an unparseable header.
fn quote_arg(a: &str) -> String {
    let clean: String = a
        .chars()
        .map(|c| match c {
            '\n' | '\r' => ' ',
            '{' => '(',
            '}' => ')',
            c => c,
        })
        .collect();
    if clean.is_empty() || clean.contains(|c: char| c.is_whitespace() || c == '\'' || c == '"') {
        format!("'{}'", clean.replace('\'', "'\\''"))
    } else {
        clean
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Envi,
    Ipw,
    Tiff,
    Png,
}

/// Identify a file by content first, extension second.
///
/// Content wins because the fixtures are inconsistent about extensions: the
/// Case 1 ENVI input is named `temp_byte_bip` with no extension at all.
pub fn detect(path: &Path) -> Result<Format> {
    let mut head = [0u8; 16];
    {
        use std::io::Read;
        let mut f = std::fs::File::open(path)
            .map_err(|e| IoError::new(format!("can't open {}: {e}", path.display())))?;
        let n = f.read(&mut head).unwrap_or(0);
        let head = &head[..n];
        if ipw::sniff(head) {
            return Ok(Format::Ipw);
        }
        if png::sniff(head) {
            return Ok(Format::Png);
        }
        if tiff::sniff(head) {
            return Ok(Format::Tiff);
        }
    }
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        match ext.to_lowercase().as_str() {
            "ipw" => return Ok(Format::Ipw),
            "tif" | "tiff" => return Ok(Format::Tiff),
            "png" => return Ok(Format::Png),
            _ => {}
        }
    }
    // Anything with a readable ENVI sidecar is ENVI.
    let hdr = envi::header_path(path);
    if hdr.exists() {
        return Ok(Format::Envi);
    }
    Err(IoError::new(format!(
        "can't identify the format of {} -- not IPW, TIFF or PNG, and no ENVI .hdr \
         sidecar found (looked for {})",
        path.display(),
        hdr.display()
    )))
}

/// Read a region map produced by an earlier run (or by another program).
///
/// ENVI only for now: it is the container both the golden and the stage-2
/// fixtures use, and the one our own writer emits by default.
pub fn read_region_map(path: &Path) -> Result<envi::RegionMapImage> {
    match detect(path)? {
        Format::Envi => envi::read_region_map(path),
        other => Err(IoError::new(format!(
            "{}: region maps can only be read from ENVI, not {other:?}",
            path.display()
        ))),
    }
}

pub fn read(path: &Path) -> Result<Image> {
    Ok(read_with_nodata(path)?.0)
}

/// Read an image, along with any nodata value the format itself declares.
///
/// ENVI carries it as `data ignore value`, GeoTIFF as the `GDAL_NODATA` tag.
/// A `--nodata` on the command line overrides whatever comes back here.
pub fn read_with_nodata(path: &Path) -> Result<(Image, Option<f64>)> {
    match detect(path)? {
        Format::Envi => {
            let hdr = envi::read_header(&envi::header_path(path))?;
            let nd = hdr.data_ignore_value;
            Ok((envi::read(path)?, nd))
        }
        Format::Ipw => Ok((ipw::read(path)?, None)),
        Format::Tiff => {
            let r = tiff::read(path)?;
            Ok((r.image, r.nodata))
        }
        Format::Png => Ok((png::read(path)?, None)),
    }
}
