"""Drive the stage-2 oracle (segment.py wind_up) on in-memory arrays.

Nothing here changes the algorithm: it calls wind_up() exactly as main.py does,
but skips the file IO so fixtures can be produced from numpy crops.
"""
import io, sys, contextlib, random
import numpy as np
import region as region_mod
import segment as segment_mod


def run(region_map, image, min_region_size, max_region_size, seed=None, quiet=True):
    """Return (new_region_map, num_passes, log_text, stats)."""
    region_map = np.array(region_map, dtype=np.uint32, copy=True)
    image = np.array(image, copy=True)
    if image.ndim == 2:
        image = image[None, :, :]
    if seed is not None:
        random.seed(seed)
    region_mod.STATS["ties"] = 0
    region_mod.STATS["cmp"] = 0
    region_mod.STATS["passes"] = []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf if quiet else sys.stdout):
        out_map, num_passes = segment_mod.wind_up(
            region_map, image, min_region_size, max_region_size
        )
    return out_map, num_passes, buf.getvalue(), dict(region_mod.STATS)
