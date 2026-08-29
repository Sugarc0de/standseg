#!/usr/bin/env python3
"""Build a large synthetic scene for latency testing by tiling a real one.

Tiling keeps realistic local spectral statistics (so region counts and pass
counts stay plausible) without needing a real 15000^2 acquisition. Per-tile
offsets stop every tile from being pixel-identical.
"""
import sys, numpy as np

src, out, size, nbands = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

raw = open(src, 'rb').read()
off = raw.rindex(b'\x0c') + 2 if b'\x0c' in raw else 0
a = np.frombuffer(raw[off:off + 250*250*8], dtype=np.uint8).reshape(250, 250, 8)[:, :, :nbands]

reps = (size + 249) // 250
big = np.empty((size, size, nbands), dtype=np.uint8)
rng = np.random.default_rng(0)
for ty in range(reps):
    for tx in range(reps):
        y0, x0 = ty*250, tx*250
        y1, x1 = min(y0+250, size), min(x0+250, size)
        # Small per-tile offset so tiles are not byte-identical.
        d = int(rng.integers(-6, 7))
        big[y0:y1, x0:x1] = np.clip(a[:y1-y0, :x1-x0].astype(np.int16) + d, 0, 255).astype(np.uint8)

big.tofile(out)
open(out + '.hdr', 'w').write(
    f"ENVI\ndescription = {{\n{out}}}\nsamples = {size}\nlines   = {size}\n"
    f"bands   = {nbands}\nheader offset = 0\nfile type = ENVI Standard\n"
    f"data type = 1\ninterleave = bip\nbyte order = 0\n")
print(f"wrote {out}: {size}x{size}x{nbands} = {big.nbytes/1e9:.2f} GB")
