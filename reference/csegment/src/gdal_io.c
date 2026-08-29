/*
 * gdal_io.c -- minimal, non-GDAL replacement.
 *
 * GDAL is not available on this machine and the original file was the only
 * GDAL consumer in the program.  This provides the same three entry points
 * with behaviour identical to the original for the raw-ENVI case that the
 * golden fixtures use:
 *
 *   GDAL_read_image     -- raw uint8 BIP file -> the same dope-vectored
 *                          uchar_t ** the original built.
 *   GDAL_write_image    -- same nbits/nbytes ladder as the original, raw
 *                          little-endian pixels row-major + an ENVI .hdr.
 *   GDAL_process_headers-- fills nlines/nsamps/nbands from the .hdr sidecar
 *                          and calls the read.
 *
 * Projection / geotransform are kept as harmless no-ops so that the rest of
 * the program (which only prints them) is untouched.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "segment.h"

/* ENVI "data type" codes we care about */
#define ENVI_DT_BYTE    1
#define ENVI_DT_UINT16  12
#define ENVI_DT_UINT32  13

typedef struct {
    int  samples;
    int  lines;
    int  bands;
    int  data_type;
    char interleave[8];         /* "bip" / "bsq" / "bil" */
    int  header_offset;
    int  byte_order;
} envi_hdr;


/*
 * Locate the ENVI header sidecar for an image file.  GDAL's ENVI driver
 * tries "<file>.hdr" first and then "<base>.hdr".  Returns a malloc'd path
 * or NULL.
 */
static char *
find_hdr(const char *fname)
{
    char           *path;
    FILE           *fp;
    char           *dot;
    size_t          n = strlen(fname);

    path = (char *) malloc(n + 8);
    if (path == NULL)
        error("out of memory");

    (void) sprintf(path, "%s.hdr", fname);
    if ((fp = fopen(path, "r")) != NULL) {
        (void) fclose(fp);
        return (path);
    }

    (void) strcpy(path, fname);
    dot = strrchr(path, '.');
    if (dot != NULL && strchr(dot, '/') == NULL) {
        (void) strcpy(dot, ".hdr");
        if ((fp = fopen(path, "r")) != NULL) {
            (void) fclose(fp);
            return (path);
        }
    }

    free(path);
    return (NULL);
}


static void
lower_str(char *s)
{
    for (; *s != '\0'; s++)
        *s = (char) tolower((unsigned char) *s);
}


/*
 * Parse the handful of ENVI header fields the program needs.  Records are
 * "key = value"; brace-delimited values may span lines but none of the keys
 * we read use them, so we simply skip over any line we do not recognise.
 */
static void
read_envi_hdr(const char *fname, envi_hdr * h)
{
    FILE           *fp;
    char           *hdrpath;
    char            line[4096];

    h->samples = h->lines = h->bands = 0;
    h->data_type = 0;
    h->header_offset = 0;
    h->byte_order = 0;
    (void) strcpy(h->interleave, "bip");

    hdrpath = find_hdr(fname);
    if (hdrpath == NULL)
        error("Can't find ENVI header (.hdr) for image file \"%s\"", fname);

    fp = fopen(hdrpath, "r");
    if (fp == NULL)
        error("Can't open ENVI header \"%s\"", hdrpath);

    while (fgets(line, sizeof(line), fp) != NULL) {
        char           *eq;
        char           *key;
        char           *val;
        char           *p;

        eq = strchr(line, '=');
        if (eq == NULL)
            continue;
        *eq = '\0';
        key = line;
        val = eq + 1;

        /* trim + fold the key */
        while (*key == ' ' || *key == '\t')
            key++;
        p = key + strlen(key);
        while (p > key && (p[-1] == ' ' || p[-1] == '\t'))
            *--p = '\0';
        lower_str(key);

        while (*val == ' ' || *val == '\t')
            val++;
        p = val + strlen(val);
        while (p > val && (p[-1] == '\n' || p[-1] == '\r' ||
                           p[-1] == ' ' || p[-1] == '\t'))
            *--p = '\0';

        if (strcmp(key, "samples") == 0)
            h->samples = atoi(val);
        else if (strcmp(key, "lines") == 0)
            h->lines = atoi(val);
        else if (strcmp(key, "bands") == 0)
            h->bands = atoi(val);
        else if (strcmp(key, "data type") == 0)
            h->data_type = atoi(val);
        else if (strcmp(key, "header offset") == 0)
            h->header_offset = atoi(val);
        else if (strcmp(key, "byte order") == 0)
            h->byte_order = atoi(val);
        else if (strcmp(key, "interleave") == 0) {
            lower_str(val);
            (void) strncpy(h->interleave, val, sizeof(h->interleave) - 1);
            h->interleave[sizeof(h->interleave) - 1] = '\0';
        }
    }
    (void) fclose(fp);
    free(hdrpath);

    if (h->samples <= 0 || h->lines <= 0 || h->bands <= 0)
        error("Bad ENVI header for \"%s\": samples/lines/bands missing", fname);
}


/*
 * Read a raw uint8 image into the dope-vectored, band-interleaved-by-pixel
 * layout the rest of the program expects.  The allocation is byte-for-byte
 * the one the original made.
 */
uchar_t **
GDAL_read_image(const char *fname, int nlines, int nsamps, int nbands,
                const char *interleave, int header_offset)
{
    uchar_t       **image;
    int             image_size;
    int             line_size;
    int             stored_image_size;
    uchar_t        *row_p;
    uchar_t        *raw;
    FILE           *fp;
    int             line;

    line_size = nbands * nsamps;
    image_size = line_size * nlines;
    stored_image_size = image_size + nlines * sizeof(uchar_t *);

    image = (uchar_t **) LINT_CAST(ecalloc(stored_image_size, 1));
    if (image == NULL)
        error("Can't allocate space for image");

    /* Initialize image scanline pointers */
    row_p = (uchar_t *) LINT_CAST(image + nlines);
    for (line = 0; line < nlines; line++) {
        image[line] = row_p;
        row_p += line_size;
    }

    printf("Trying to read in image\n");

    fp = fopen(fname, "rb");
    if (fp == NULL)
        error("Can't open input image file \"%s\"", fname);
    if (header_offset > 0 && fseek(fp, (long) header_offset, SEEK_SET) != 0)
        error("Can't seek past header offset in \"%s\"", fname);

    if (strcmp(interleave, "bip") == 0) {
        /* On-disk layout already matches the in-memory layout. */
        if (fread(image[0], 1, (size_t) image_size, fp) != (size_t) image_size)
            error("Short read on input image file \"%s\"", fname);
    } else {
        raw = (uchar_t *) malloc((size_t) image_size);
        if (raw == NULL)
            error("Can't allocate space for image");
        if (fread(raw, 1, (size_t) image_size, fp) != (size_t) image_size)
            error("Short read on input image file \"%s\"", fname);

        if (strcmp(interleave, "bsq") == 0) {
            int             b, l, s;
            for (b = 0; b < nbands; b++)
                for (l = 0; l < nlines; l++)
                    for (s = 0; s < nsamps; s++)
                        image[l][s * nbands + b] =
                            raw[((b * nlines) + l) * nsamps + s];
        } else if (strcmp(interleave, "bil") == 0) {
            int             b, l, s;
            for (l = 0; l < nlines; l++)
                for (b = 0; b < nbands; b++)
                    for (s = 0; s < nsamps; s++)
                        image[l][s * nbands + b] =
                            raw[(l * nbands + b) * nsamps + s];
        } else {
            error("Unsupported ENVI interleave \"%s\"", interleave);
        }
        free(raw);
    }
    (void) fclose(fp);

    printf("Read in data\n");

    return (image);
}


/*
 * Write the region band.  The nbits ladder is exactly the original's; the
 * pixels are the low nbytes of each REGION_ID, little-endian, row-major --
 * which is what GDAL's ENVI driver produced.
 */
void
GDAL_write_image(Seg_proc Spr, char *fname)
{
    int             nbits, nbytes;
    long            nregions;
    int             envi_dt;
    FILE           *fp;
    char           *hdrname;
    char           *dot;
    int             l;
    int             c;
    unsigned char  *rowbuf;

    /* Figure out minimum datatype to use */
    for (nbits = 0, nregions = Spr->nreg; nregions; nbits++, nregions >>= 1);
    if (nbits <= 8) {
        envi_dt = ENVI_DT_BYTE;
        nbytes = 1;
    } else if (nbits <= 16) {
        envi_dt = ENVI_DT_UINT16;
        nbytes = 2;
    } else if (nbits <= 32) {
        envi_dt = ENVI_DT_UINT32;
        nbytes = 4;
    } else {
        error("Cannot determine datatype\n");
        return;                 /* NOTREACHED */
    }

    fp = fopen(fname, "wb");
    if (fp == NULL)
        error("Error opening output file %s\n", fname);

    rowbuf = (unsigned char *) malloc((size_t) Spr->nsamps * nbytes);
    if (rowbuf == NULL)
        error("Can't allocate output row buffer");

    for (l = 0; l < Spr->nlines; l++) {
        for (c = 0; c < Spr->nsamps; c++) {
            REGION_ID       v = Spr->rband[l][c];
            int             k;
            for (k = 0; k < nbytes; k++)
                rowbuf[c * nbytes + k] = (unsigned char) ((v >> (8 * k)) & 0xff);
        }
        if (fwrite(rowbuf, (size_t) nbytes, (size_t) Spr->nsamps, fp)
            != (size_t) Spr->nsamps)
            error("Error writing output file %s (line %d)\n", fname, l);
    }
    free(rowbuf);
    if (fclose(fp) != 0)
        error("Error closing output file %s\n", fname);

    /*
     * Matching ENVI header.  GDAL's ENVI driver replaces the file's
     * extension with ".hdr"; do the same.
     */
    hdrname = (char *) malloc(strlen(fname) + 8);
    if (hdrname == NULL)
        error("out of memory");
    (void) strcpy(hdrname, fname);
    dot = strrchr(hdrname, '.');
    if (dot != NULL && strchr(dot, '/') == NULL)
        (void) strcpy(dot, ".hdr");
    else
        (void) strcat(hdrname, ".hdr");

    fp = fopen(hdrname, "w");
    if (fp == NULL)
        error("Error opening output header %s\n", hdrname);

    {
        const char     *base = strrchr(fname, '/');
        base = (base == NULL) ? fname : base + 1;
        (void) fprintf(fp, "ENVI\n");
        (void) fprintf(fp, "description = {\n%s}\n", base);
        (void) fprintf(fp, "samples = %d\n", Spr->nsamps);
        (void) fprintf(fp, "lines   = %d\n", Spr->nlines);
        (void) fprintf(fp, "bands   = 1\n");
        (void) fprintf(fp, "header offset = 0\n");
        (void) fprintf(fp, "file type = ENVI Standard\n");
        (void) fprintf(fp, "data type = %d\n", envi_dt);
        (void) fprintf(fp, "interleave = bsq\n");
        (void) fprintf(fp, "byte order = 0\n");
        (void) fprintf(fp, "band names = {\nBand 1}\n");
    }
    (void) fclose(fp);
    free(hdrname);
}


void
GDAL_process_headers(Seg_proc Spr)
{
    envi_hdr        h;

    read_envi_hdr(Spr->image_fn, &h);

    Spr->nlines = h.lines;
    if (Spr->nlines > MAXSHORT)
        error("Image has too many (%d) lines\n", Spr->nlines);
    Spr->nsamps = h.samples;
    if (Spr->nsamps > MAXSHORT)
        error("Image has too many (%d) samps\n", Spr->nsamps);
    Spr->nbands = h.bands;
    if (Spr->nbands > MAXSHORT)
        error("Image has too many (%d) bands\n", Spr->nbands);

    /* Projection / geotransform: harmless no-ops. */
    Spr->pszProjection = (char *) malloc(1);
    if (Spr->pszProjection == NULL)
        error("out of memory");
    Spr->pszProjection[0] = '\0';
    if (strlen(Spr->pszProjection) == 0)
        warn("Could not get image projection\n");
    printf("Projection is: %s\n", Spr->pszProjection);
    Spr->adfGeoTransform[0] = 0.0;
    Spr->adfGeoTransform[1] = 1.0;
    Spr->adfGeoTransform[2] = 0.0;
    Spr->adfGeoTransform[3] = 0.0;
    Spr->adfGeoTransform[4] = 0.0;
    Spr->adfGeoTransform[5] = 1.0;

    Spr->image = GDAL_read_image(Spr->image_fn, Spr->nlines, Spr->nsamps,
                                 Spr->nbands, h.interleave, h.header_offset);

    if (h.data_type != ENVI_DT_BYTE) {
        error("Image must be Byte datatype, not data type %d\n", h.data_type);
    }

    /*
     * Mask image, if there is one.
     */
    if (strlen(Spr->mask_fn) > 0) {
        envi_hdr        mh;

        read_envi_hdr(Spr->mask_fn, &mh);

        if (mh.lines != Spr->nlines)
            error("Input and mask images have different number of lines");
        if (mh.samples != Spr->nsamps)
            error("Input and mask images have different number of samples");
        if (mh.bands != 1)
            error("Input mask image must have 1 band");
        if (mh.data_type != ENVI_DT_BYTE)
            error("The mask image is not 1 byte per pixel");

        Spr->imask = GDAL_read_image(Spr->mask_fn, mh.lines, mh.samples,
                                     mh.bands, mh.interleave,
                                     mh.header_offset);
    }
}
