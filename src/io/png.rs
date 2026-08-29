//! PNG.
//!
//! 8- and 16-bit samples. Channel count becomes the band count: grayscale is 1
//! band, RGB is 3, RGBA is 4 -- note that an alpha channel is read as an
//! ordinary band and will take part in the spectral distance, which is almost
//! never what you want. Use `--nodata` or a mask for transparency instead.
//!
//! Sub-byte depths (1, 2, 4) are rejected rather than expanded: they are palette
//! or bilevel images, where the stored number is an index, not a radiance.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::image::{GeoRef, Image, Samples};
use crate::io::{IoError, Result};

pub fn read(path: &Path) -> Result<Image> {
    let f = File::open(path)
        .map_err(|e| IoError::new(format!("can't open {}: {e}", path.display())))?;
    let dec = png::Decoder::new(BufReader::new(f));
    let mut reader = dec
        .read_info()
        .map_err(|e| IoError::new(format!("{}: not a readable PNG: {e}", path.display())))?;

    let info = reader.info();
    let sample_bytes = match info.bit_depth {
        png::BitDepth::Eight => 1usize,
        png::BitDepth::Sixteen => 2usize,
        other => {
            return Err(IoError::new(format!(
                "{}: PNG is {:?} bits per sample; this program segments 8- and \
                 16-bit imagery only",
                path.display(),
                other
            )))
        }
    };
    let (w, h) = (info.width as usize, info.height as usize);

    let mut buf = vec![0u8; reader.output_buffer_size().unwrap_or(0)];
    let out = reader
        .next_frame(&mut buf)
        .map_err(|e| IoError::new(format!("{}: can't decode: {e}", path.display())))?;
    let nbands = out.color_type.samples();
    buf.truncate(out.buffer_size());

    if buf.len() != w * h * nbands * sample_bytes {
        return Err(IoError::new(format!(
            "{}: decoded {} bytes, expected {}x{}x{} at {} bytes/sample",
            path.display(),
            buf.len(),
            w,
            h,
            nbands,
            sample_bytes
        )));
    }

    // PNG rows are already interleaved by pixel, i.e. BIP. 16-bit samples are
    // big-endian on the wire and in the decoder's output buffer.
    let data = if sample_bytes == 1 {
        Samples::U8(buf)
    } else {
        Samples::U16(
            buf.chunks_exact(2)
                .map(|c| u16::from_be_bytes([c[0], c[1]]))
                .collect(),
        )
    };
    let mut image = Image::from_samples(h, w, nbands, data);
    image.geo = GeoRef {
        description: path.file_name().map(|s| s.to_string_lossy().to_string()),
        ..Default::default()
    };
    Ok(image)
}

pub fn sniff(head: &[u8]) -> bool {
    head.starts_with(b"\x89PNG\r\n\x1a\n")
}
