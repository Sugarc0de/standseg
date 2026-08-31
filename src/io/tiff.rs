//! TIFF / GeoTIFF.
//!
//! Same sample widths as every other reader here: 8- and 16-bit integers, plus
//! 32-bit float for stage-2 layers (height, biomass, age, z-scores -- the
//! structural imagery the second stage segments against, which is routinely
//! float). Bands map to TIFF samples-per-pixel, so an RGB TIFF is a 3-band image
//! and a 6-band satellite stack is a 6-sample TIFF.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use tiff::decoder::{Decoder, DecodingResult, Limits};
use tiff::tags::Tag;

use crate::image::{GeoRef, Image, Samples};
use crate::io::{IoError, Provenance, Result};

/// GDAL writes its nodata value into this private tag, as an ASCII string.
const GDAL_NODATA: u16 = 42113;

pub struct TiffRead {
    pub image: Image,
    /// Nodata declared by the file, if any.
    pub nodata: Option<f64>,
}

pub fn read(path: &Path) -> Result<TiffRead> {
    let f = File::open(path)
        .map_err(|e| IoError::new(format!("can't open {}: {e}", path.display())))?;
    // The crate's default 256 MB decode budget is smaller than the imagery this
    // program exists for: a 5000 x 5000 x 6 uint8 stack is 150 MB and a 16-bit
    // one is 300 MB. Raise it to something that admits the scenes we advertise
    // (15000^2 x 6 at 16 bits) rather than lifting the limit entirely.
    const MAX_DECODE: usize = 3 << 30;
    let mut dec = Decoder::new(BufReader::new(f))
        .map_err(|e| IoError::new(format!("{}: not a readable TIFF: {e}", path.display())))?
        .with_limits({
            let mut l = Limits::default();
            l.decoding_buffer_size = MAX_DECODE;
            l
        });

    let (w, h) = dec
        .dimensions()
        .map_err(|e| IoError::new(format!("{}: no dimensions: {e}", path.display())))?;

    // `read_image` reads *one plane*. For PlanarConfiguration = 2 -- separate
    // planes, which is how GDAL writes some multi-band stacks -- that is the
    // first band, and the result then divides evenly into a "1-band image":
    // silently the wrong answer on a 6-band file. Read into a buffer we can
    // measure against `complete_len` instead, and refuse anything short.
    let mut img = DecodingResult::U8(Vec::new());
    let layout = dec
        .read_image_to_buffer(&mut img)
        .map_err(|e| IoError::new(format!("{}: can't decode: {e}", path.display())))?;
    let got = img.as_buffer(0).as_bytes().len();
    if got < layout.complete_len {
        return Err(IoError::new(format!(
            "{}: only {got} of {} bytes were decoded ({} planes); the file is \
             larger than this program's {} MB decode budget",
            path.display(),
            layout.complete_len,
            layout.planes,
            MAX_DECODE >> 20
        )));
    }
    let nplanes = layout.planes;

    let data = match img {
        DecodingResult::U8(v) => Samples::U8(v),
        DecodingResult::U16(v) => Samples::U16(v),
        DecodingResult::I16(v) => Samples::I16(v),
        DecodingResult::F32(v) => Samples::F32(v),
        other => {
            return Err(IoError::new(format!(
                "{}: TIFF samples are {}; this program reads 8- and 16-bit \
                 integer imagery, and 32-bit float for --stage2 layers",
                path.display(),
                sample_kind(&other)
            )))
        }
    };

    let npix = w as usize * h as usize;
    if npix == 0 || data.len() % npix != 0 {
        return Err(IoError::new(format!(
            "{}: {} samples does not divide {}x{} pixels evenly",
            path.display(),
            data.len(),
            w,
            h
        )));
    }
    let nbands = data.len() / npix;

    // Chunky TIFFs are already BIP. Planar ones arrive plane by plane -- BSQ --
    // and everything downstream indexes pixels, so transpose here rather than
    // teaching the segmenter a second layout.
    let data = if nplanes > 1 {
        if nplanes != nbands {
            return Err(IoError::new(format!(
                "{}: {nplanes} planes but {nbands} bands; cannot deinterleave",
                path.display()
            )));
        }
        bsq_to_bip(data, npix, nbands)
    } else {
        data
    };

    let mut image = Image::from_samples(h as usize, w as usize, nbands, data);

    let nodata = dec
        .get_tag_ascii_string(Tag::Unknown(GDAL_NODATA))
        .ok()
        .and_then(|s| s.trim().parse::<f64>().ok());

    image.geo = GeoRef {
        description: path.file_name().map(|s| s.to_string_lossy().to_string()),
        ..Default::default()
    };

    Ok(TiffRead { image, nodata })
}

/// Band-sequential to band-interleaved-by-pixel.
fn bsq_to_bip(data: Samples, npix: usize, nbands: usize) -> Samples {
    fn t<T: Copy + Default>(v: &[T], npix: usize, nbands: usize) -> Vec<T> {
        let mut out = vec![T::default(); v.len()];
        for b in 0..nbands {
            let plane = &v[b * npix..(b + 1) * npix];
            for (p, &s) in plane.iter().enumerate() {
                out[p * nbands + b] = s;
            }
        }
        out
    }
    match data {
        Samples::U8(v) => Samples::U8(t(&v, npix, nbands)),
        Samples::U16(v) => Samples::U16(t(&v, npix, nbands)),
        Samples::I16(v) => Samples::I16(t(&v, npix, nbands)),
        Samples::F32(v) => Samples::F32(t(&v, npix, nbands)),
    }
}

fn sample_kind(r: &DecodingResult) -> &'static str {
    match r {
        DecodingResult::U8(_) => "8-bit unsigned",
        DecodingResult::U16(_) => "16-bit unsigned",
        DecodingResult::U32(_) => "32-bit unsigned",
        DecodingResult::U64(_) => "64-bit unsigned",
        DecodingResult::I8(_) => "8-bit signed",
        DecodingResult::I16(_) => "16-bit signed",
        DecodingResult::I32(_) => "32-bit signed",
        DecodingResult::I64(_) => "64-bit signed",
        DecodingResult::F16(_) => "16-bit float",
        DecodingResult::F32(_) => "32-bit float",
        DecodingResult::F64(_) => "64-bit float",
    }
}

pub fn sniff(head: &[u8]) -> bool {
    head.starts_with(b"II\x2a\x00") || head.starts_with(b"MM\x00\x2a")
}

/// Write a single-band region map as a TIFF.
///
/// The command that produced the file goes into `ImageDescription` and the
/// program into `Software`, which is the TIFF equivalent of IPW's `history`
/// record.
pub fn write_region_map(
    path: &Path,
    rband: &[u32],
    nlines: usize,
    nsamps: usize,
    nbytes: usize,
    prov: &Provenance,
) -> Result<()> {
    use tiff::encoder::colortype::ColorType;
    use tiff::encoder::{colortype, TiffEncoder};
    use tiff::tags::Tag as T;

    let f = File::create(path)
        .map_err(|e| IoError::new(format!("can't create {}: {e}", path.display())))?;
    let mut enc =
        TiffEncoder::new(f).map_err(|e| IoError::new(format!("{}: {e}", path.display())))?;

    let (w, h) = (nsamps as u32, nlines as u32);
    let n = nlines * nsamps;

    /// One `write_image` with the provenance tags attached first.
    fn tagged<W, K, C>(
        enc: &mut TiffEncoder<W, K>,
        w: u32,
        h: u32,
        data: &[C::Inner],
        prov: &Provenance,
    ) -> tiff::TiffResult<()>
    where
        W: std::io::Write + std::io::Seek,
        K: tiff::encoder::TiffKind,
        C: ColorType,
        [C::Inner]: tiff::encoder::TiffValue,
    {
        let mut img = enc.new_image::<C>(w, h)?;
        if !prov.software.is_empty() {
            img.encoder()
                .write_tag(T::Software, prov.software.as_str())?;
        }
        if !prov.command.is_empty() {
            img.encoder()
                .write_tag(T::ImageDescription, prov.command.as_str())?;
        }
        img.write_data(data)
    }

    let res = match nbytes {
        1 => {
            let v: Vec<u8> = rband[..n].iter().map(|&x| x as u8).collect();
            tagged::<_, _, colortype::Gray8>(&mut enc, w, h, &v, prov)
        }
        2 => {
            let v: Vec<u16> = rband[..n].iter().map(|&x| x as u16).collect();
            tagged::<_, _, colortype::Gray16>(&mut enc, w, h, &v, prov)
        }
        4 => tagged::<_, _, colortype::Gray32>(&mut enc, w, h, &rband[..n], prov),
        other => return Err(IoError::new(format!("bad region map pixel size {other}"))),
    };
    res.map_err(|e| IoError::new(format!("can't write {}: {e}", path.display())))
}
