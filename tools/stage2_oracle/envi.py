"""Minimal ENVI writer, matching the header conventions this repo emits."""
import numpy as np, os

DT = {np.dtype('uint8'): 1, np.dtype('int16'): 2, np.dtype('uint16'): 12,
      np.dtype('uint32'): 13}


def write(path, arr, map_info=None, coord_sys=None, ignore0=False, desc=None,
          band_names=None):
    """arr: (rows, cols) or (bands, rows, cols). Written BSQ, little-endian."""
    if arr.ndim == 2:
        arr = arr[None, :, :]
    nb, nl, ns = arr.shape
    arr = np.ascontiguousarray(arr)
    with open(path, 'wb') as f:
        f.write(arr.tobytes(order='C'))          # BSQ: band, row, col
    name = os.path.basename(path)
    h = ["ENVI", "description = {", f"{desc or name}}}",
         f"samples = {ns}", f"lines   = {nl}", f"bands   = {nb}",
         "header offset = 0", "file type = ENVI Standard",
         f"data type = {DT[arr.dtype]}", "interleave = bsq", "byte order = 0"]
    if ignore0:
        h.append("data ignore value = 0")
    if map_info:
        h.append("map info = {%s}" % map_info)
    if coord_sys:
        h.append("coordinate system string = {%s}" % coord_sys)
    names = band_names or [f"Band {i+1}" for i in range(nb)]
    h.append("band names = {\n" + ",\n".join(names) + "}")
    with open(path + ".hdr", 'w') as f:
        f.write("\n".join(h) + "\n")


def crop_map_info(map_info, r0, c0):
    """Shift an ENVI `map info` tie point to the origin of a crop."""
    p = [x.strip() for x in map_info.split(',')]
    # proj, ref_x(1-based), ref_y(1-based), east, north, xsize, ysize, datum...
    east, north = float(p[3]), float(p[4])
    xs, ys = float(p[5]), float(p[6])
    p[3] = f"{east + c0 * xs:.4f}"
    p[4] = f"{north - r0 * ys:.4f}"
    return ", ".join(p[:7]) + "," + ",".join(p[7:])
