#!/usr/bin/env python3
"""Regenerate docs/img/stands-landsat.png -- the README's worked example.

Draws the stand boundaries standseg produced over the Landsat composite it
produced them from. Run it from the repository root, after the worked example
in the README:

    standseg -t 20 -m .2 -n 50,100,200 --format gpkg \\
        -o stands --outdir out tests/golden/misc/temp_byte_bip
    python3 docs/make_figure.py out/stands.armap.69

Needs numpy and pillow. Nothing in the program itself needs Python; this is a
documentation tool.
"""
import sys

import numpy as np
from PIL import Image

SCENE = "tests/golden/misc/temp_byte_bip"
OUT = "docs/img/stands-landsat.png"
N = 250
ZOOM = 8
LINE = 3
BOUNDARY = (255, 40, 40)


def stretch(band, lo=2, hi=98):
    """Percentile stretch to 0-255, the usual way a composite is made legible."""
    a, b = np.percentile(band, [lo, hi])
    return np.clip((band.astype(np.float32) - a) * 255.0 / max(b - a, 1e-6), 0, 255)


def read_region_map(path):
    """An ENVI region map is 1, 2 or 4 bytes per pixel; the header says which."""
    hdr = {}
    for line in open(f"{path}.hdr"):
        if "=" in line:
            k, v = line.split("=", 1)
            hdr[k.strip()] = v.strip()
    dtype = {"1": np.uint8, "12": np.uint16, "13": np.uint32}[hdr["data type"]]
    return np.fromfile(path, dtype=dtype).reshape(int(hdr["lines"]), int(hdr["samples"]))


def main(rmap_path):
    # Bands are Red, NIR, SWIR1, SWIR2 (Landsat 8 OLI 4, 5, 6, 7), interleaved
    # by pixel.
    scene = np.fromfile(SCENE, dtype=np.uint8).reshape(N, N, 4)
    # SWIR1 / NIR / Red: the standard forestry composite, vegetation in green.
    rgb = np.dstack([stretch(scene[:, :, i]) for i in (2, 1, 0)]).astype(np.uint8)

    img = np.asarray(Image.fromarray(rgb).resize((N * ZOOM, N * ZOOM), Image.NEAREST)).copy()

    regions = read_region_map(rmap_path)

    # A boundary is where the region id changes. Drawn on the shared edge, so
    # neighbouring stands share one line rather than each drawing their own.
    vert = regions[:, 1:] != regions[:, :-1]
    horiz = regions[1:, :] != regions[:-1, :]
    for row, col in zip(*np.nonzero(vert)):
        x = (col + 1) * ZOOM
        img[row * ZOOM : (row + 1) * ZOOM, x - LINE // 2 : x + LINE - LINE // 2] = BOUNDARY
    for row, col in zip(*np.nonzero(horiz)):
        y = (row + 1) * ZOOM
        img[y - LINE // 2 : y + LINE - LINE // 2, col * ZOOM : (col + 1) * ZOOM] = BOUNDARY

    Image.fromarray(img).save(OUT)
    print(f"{OUT}: {img.shape[1]} x {img.shape[0]}, {len(np.unique(regions))} stands")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "out/stands.armap.69")
