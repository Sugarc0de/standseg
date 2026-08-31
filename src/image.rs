//! In-memory image, band-interleaved-by-pixel.
//!
//! The original stores this as a dope vector (`uchar_t **`: an array of row
//! pointers into one contiguous block). That is a C ergonomics hack, not a space
//! saving -- it costs 8 bytes per row and an indirection per access. A flat
//! buffer with computed offsets is smaller and faster. See PLAN.md section 4, trick 1.
//!
//! # Sample width
//!
//! The 1992 program was uint8-only, because that is what an 8-bit TM scene was.
//! Landsat 8/9 and Sentinel-2 are 12-bit data delivered in 16-bit containers, so
//! forcing them down to 8 bits throws away radiometry *before* segmentation and
//! changes the answer. We therefore carry three sample widths -- `u8`, `u16` and
//! `i16` (the container Landsat Collection 2 surface reflectance ships in).
//!
//! Nothing downstream of Phase 0 sees a pixel: the region centroids are f32 and
//! the image is freed as soon as the region list exists. So widening the input
//! is confined to this file, `pixel.rs` and `RegionList::from_pixel`, and the
//! `u8` path keeps its exact arithmetic -- which the golden fixtures pin.
//!
//! # Why there are two traits
//!
//! Stage 2 (Ye et al. 2025, `stage2.rs`) reads a *second* image that is not
//! reflectance -- height, biomass, age, a z-score -- and those ship as 32-bit
//! float. Stage 1 cannot: `pix_dist2` is integer arithmetic and the normality
//! band is compared against integer DN limits, both pinned by the golden
//! fixtures. So `Sample` carries what both stages need and `IntSample` carries
//! what only stage 1 needs. A float image therefore *cannot* reach stage 1 --
//! not by a runtime check that might be forgotten, but because `f32` does not
//! implement the trait stage 1 is generic over.

/// A pixel sample type. Everything both stages can read.
pub trait Sample: Copy + PartialEq + Send + Sync + 'static {
    /// Name used in diagnostics.
    const KIND: &'static str;
    /// Whether region means must be accumulated in the oracle's own order.
    ///
    /// False for the integer widths: every partial sum is an exact integer well
    /// below 2^53, so order cannot change the result and one `i64` sum will do.
    /// True for `f32`, where it demonstrably can -- see `stage2::mean_f32`.
    const ORDERED_MEAN: bool;
    fn to_f64(self) -> f64;
    /// Exact comparison against zero, which is what the oracle's
    /// `np.all(b_images == 0, axis=0)` does. NaN is not zero, in both.
    fn is_zero(self) -> bool;
}

/// The integer sample widths -- the only ones stage 1 accepts.
pub trait IntSample: Sample {
    const MIN_VALUE: i64;
    const MAX_VALUE: i64;
    fn to_i64(self) -> i64;
    fn to_f32(self) -> f32;
}

impl Sample for u8 {
    const KIND: &'static str = "8-bit unsigned";
    const ORDERED_MEAN: bool = false;
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0
    }
}

impl IntSample for u8 {
    const MIN_VALUE: i64 = 0;
    const MAX_VALUE: i64 = u8::MAX as i64;
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }
}

impl Sample for u16 {
    const KIND: &'static str = "16-bit unsigned";
    const ORDERED_MEAN: bool = false;
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0
    }
}

impl IntSample for u16 {
    const MIN_VALUE: i64 = 0;
    const MAX_VALUE: i64 = u16::MAX as i64;
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }
}

impl Sample for i16 {
    const KIND: &'static str = "16-bit signed";
    const ORDERED_MEAN: bool = false;
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0
    }
}

impl IntSample for i16 {
    const MIN_VALUE: i64 = i16::MIN as i64;
    const MAX_VALUE: i64 = i16::MAX as i64;
    #[inline]
    fn to_i64(self) -> i64 {
        self as i64
    }
    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }
}

/// 32-bit float -- stage 2 only. Deliberately no `IntSample`: see the module
/// docs. Note there is no `MIN_VALUE`/`MAX_VALUE` to state, which is the point.
impl Sample for f32 {
    const KIND: &'static str = "32-bit float";
    const ORDERED_MEAN: bool = true;
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0.0
    }
}

/// Pixel storage, one variant per supported sample width. BIP order throughout:
/// `nlines * nsamps * nbands` samples, bands adjacent within a pixel.
#[derive(Debug, Clone, PartialEq)]
pub enum Samples {
    U8(Vec<u8>),
    U16(Vec<u16>),
    I16(Vec<i16>),
    /// Stage-2 imagery only; stage 1 refuses it at the door.
    F32(Vec<f32>),
}

impl Samples {
    pub fn len(&self) -> usize {
        match self {
            Samples::U8(v) => v.len(),
            Samples::U16(v) => v.len(),
            Samples::I16(v) => v.len(),
            Samples::F32(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Bytes of memory per sample -- what the memory report needs.
    pub fn bytes_per_sample(&self) -> usize {
        match self {
            Samples::U8(_) => 1,
            Samples::U16(_) | Samples::I16(_) => 2,
            Samples::F32(_) => 4,
        }
    }

    pub fn kind(&self) -> &'static str {
        match self {
            Samples::U8(_) => <u8 as Sample>::KIND,
            Samples::U16(_) => <u16 as Sample>::KIND,
            Samples::I16(_) => <i16 as Sample>::KIND,
            Samples::F32(_) => <f32 as Sample>::KIND,
        }
    }

    pub fn as_u8(&self) -> Option<&[u8]> {
        match self {
            Samples::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_u16(&self) -> Option<&[u16]> {
        match self {
            Samples::U16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_i16(&self) -> Option<&[i16]> {
        match self {
            Samples::I16(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_f32(&self) -> Option<&[f32]> {
        match self {
            Samples::F32(v) => Some(v),
            _ => None,
        }
    }

    /// True for sample types stage 1 cannot segment.
    pub fn is_float(&self) -> bool {
        matches!(self, Samples::F32(_))
    }

    /// Inclusive range an integer value must sit in to be a valid sample of this
    /// type -- used to validate `--nodata` against the input.
    ///
    /// For `F32` this is the range over which integers are *exactly*
    /// representable (±2^24); beyond it two adjacent integers share a float, so
    /// an equality test against one of them would be meaningless. Only the
    /// stage-1 path calls this, and that path never sees a float image.
    /// Smallest and largest sample actually present. Diagnostics only -- the
    /// algorithm never consults it. It exists so the CLI can tell a user that
    /// their tolerance is in 8-bit units and their image is not.
    pub fn observed_range(&self) -> Option<(f64, f64)> {
        fn mm(mut it: impl Iterator<Item = f64>) -> Option<(f64, f64)> {
            let first = it.next()?;
            let (mut lo, mut hi) = (first, first);
            for v in it {
                if v < lo {
                    lo = v;
                } else if v > hi {
                    hi = v;
                }
            }
            Some((lo, hi))
        }
        match self {
            Samples::U8(v) => mm(v.iter().map(|&s| s as f64)),
            Samples::U16(v) => mm(v.iter().map(|&s| s as f64)),
            Samples::I16(v) => mm(v.iter().map(|&s| s as f64)),
            Samples::F32(v) => mm(v.iter().map(|&s| s as f64)),
        }
    }

    pub fn value_range(&self) -> (i64, i64) {
        match self {
            Samples::U8(_) => (<u8 as IntSample>::MIN_VALUE, <u8 as IntSample>::MAX_VALUE),
            Samples::U16(_) => (<u16 as IntSample>::MIN_VALUE, <u16 as IntSample>::MAX_VALUE),
            Samples::I16(_) => (<i16 as IntSample>::MIN_VALUE, <i16 as IntSample>::MAX_VALUE),
            Samples::F32(_) => (-(1 << 24), 1 << 24),
        }
    }
}

/// Georeferencing carried through from input to output where the format has it.
#[derive(Debug, Clone, Default)]
pub struct GeoRef {
    pub map_info: Option<String>,
    pub coord_sys: Option<String>,
    pub band_names: Vec<String>,
    pub description: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Image {
    pub nlines: usize,
    pub nsamps: usize,
    pub nbands: usize,
    pub data: Samples,
    pub geo: GeoRef,
}

/// A typed, borrowed view of the pixel buffer. Phase 0 is generic over this so
/// each sample width gets its own specialised code, rather than paying a match
/// per pixel access.
pub struct Raster<'a, T> {
    pub nlines: usize,
    pub nsamps: usize,
    pub nbands: usize,
    pub data: &'a [T],
}

impl<'a, T: Sample> Raster<'a, T> {
    /// The `nbands` values at (line, samp). Mirrors the C's `pixvec` macro.
    #[inline]
    pub fn pixel(&self, line: usize, samp: usize) -> &[T] {
        let off = (line * self.nsamps + samp) * self.nbands;
        &self.data[off..off + self.nbands]
    }

    #[inline]
    pub fn pixel_at(&self, idx: usize) -> &[T] {
        let off = idx * self.nbands;
        &self.data[off..off + self.nbands]
    }

    #[inline]
    pub fn npixels(&self) -> usize {
        self.nlines * self.nsamps
    }
}

impl Image {
    /// An 8-bit image of zeros. 8-bit stays the default because it is what the
    /// fixtures are and what most callers construct.
    pub fn new(nlines: usize, nsamps: usize, nbands: usize) -> Self {
        Self::from_samples(
            nlines,
            nsamps,
            nbands,
            Samples::U8(vec![0u8; nlines * nsamps * nbands]),
        )
    }

    pub fn from_samples(nlines: usize, nsamps: usize, nbands: usize, data: Samples) -> Self {
        debug_assert_eq!(data.len(), nlines * nsamps * nbands);
        Self {
            nlines,
            nsamps,
            nbands,
            data,
            geo: GeoRef::default(),
        }
    }

    #[inline]
    pub fn npixels(&self) -> usize {
        self.nlines * self.nsamps
    }

    /// Bytes the pixel buffer occupies.
    pub fn nbytes(&self) -> usize {
        self.data.len() * self.data.bytes_per_sample()
    }

    /// Run `f` against whichever typed view the image actually holds.
    ///
    /// This is the single dispatch point: everything below it is monomorphic.
    pub fn with_raster<R>(&self, f: impl FnOnce(RasterRef<'_>) -> R) -> R {
        let (nl, ns, nb) = (self.nlines, self.nsamps, self.nbands);
        match &self.data {
            Samples::U8(v) => f(RasterRef::U8(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
            Samples::U16(v) => f(RasterRef::U16(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
            Samples::I16(v) => f(RasterRef::I16(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
            Samples::F32(v) => f(RasterRef::F32(Raster {
                nlines: nl,
                nsamps: ns,
                nbands: nb,
                data: v,
            })),
        }
    }

    /// Every pixel whose bands match `nd` is marked 0 in `mask`.
    ///
    /// `any` picks the "nodata if any band matches" reading; the default is that
    /// all bands must match. A value outside the sample type's range simply
    /// matches nothing.
    pub fn apply_nodata(&self, nd: i64, any: bool, mask: &mut [u8]) {
        fn scan<T: Sample>(r: &Raster<'_, T>, nd: f64, any: bool, mask: &mut [u8]) {
            for (p, m) in mask.iter_mut().enumerate().take(r.npixels()) {
                let px = r.pixel_at(p);
                let hit = if any {
                    px.iter().any(|s| s.to_f64() == nd)
                } else {
                    px.iter().all(|s| s.to_f64() == nd)
                };
                if hit {
                    *m = 0;
                }
            }
        }
        let nd = nd as f64;
        self.with_raster(|r| match r {
            RasterRef::U8(r) => scan(&r, nd, any, mask),
            RasterRef::U16(r) => scan(&r, nd, any, mask),
            RasterRef::I16(r) => scan(&r, nd, any, mask),
            RasterRef::F32(r) => scan(&r, nd, any, mask),
        })
    }

    /// Flatten a single-band image to a 0/1 mask: any nonzero sample is valid.
    pub fn to_mask(&self) -> Vec<u8> {
        fn flat<T: Sample>(r: &Raster<'_, T>) -> Vec<u8> {
            r.data.iter().map(|s| u8::from(!s.is_zero())).collect()
        }
        self.with_raster(|r| match r {
            RasterRef::U8(r) => flat(&r),
            RasterRef::U16(r) => flat(&r),
            RasterRef::I16(r) => flat(&r),
            RasterRef::F32(r) => flat(&r),
        })
    }
}

/// The three typed views `Image::with_raster` can hand back.
pub enum RasterRef<'a> {
    U8(Raster<'a, u8>),
    U16(Raster<'a, u16>),
    I16(Raster<'a, i16>),
    F32(Raster<'a, f32>),
}
