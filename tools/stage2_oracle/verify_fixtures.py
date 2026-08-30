"""Independent invariant checks on tests/stage2 -- no shared code with the generator."""
import json, os, sys
import numpy as np

ROOT = "/Users/elaineye/mac2025/fast_segment/tests/stage2"
DT = {1: np.uint8, 2: np.int16, 12: np.uint16, 13: np.uint32}


def read_envi(path):
    h = {}
    for line in open(path + ".hdr"):
        if "=" in line:
            k, v = line.split("=", 1)
            h[k.strip()] = v.strip().strip("{}")
    nl, ns, nb, dt = int(h["lines"]), int(h["samples"]), int(h["bands"]), int(h["data type"])
    assert h["interleave"] == "bsq" and h["byte order"] == "0", path
    a = np.fromfile(path, dtype=DT[dt])
    assert a.size == nl * ns * nb, f"{path}: {a.size} != {nl*ns*nb}"
    return a.reshape(nb, nl, ns), dt


fail = 0
for c in json.load(open(f"{ROOT}/cases.json")):
    d = f"{ROOT}/{c['name']}"
    rm, rdt = read_envi(f"{d}/input/rmap")
    lay, ldt = read_envi(f"{d}/input/layer")
    ex, edt = read_envi(f"{d}/{c['expected']}")
    rm, ex = rm[0], ex[0]
    errs = []

    if rm.shape != tuple(c["shape"]) or ex.shape != rm.shape or lay.shape[1:] != rm.shape:
        errs.append(f"shape mismatch {rm.shape} {ex.shape} {lay.shape}")
    if lay.shape[0] != c["nbands"]:
        errs.append(f"nbands {lay.shape[0]} != {c['nbands']}")
    if np.dtype(DT[rdt]).name != c["region_map_dtype"] or rdt != edt:
        errs.append(f"dtype rmap={rdt} expected={edt} json={c['region_map_dtype']}")

    # 1. stage 2 never splits: pixels sharing an input region share an output region.
    order = np.argsort(rm.ravel(), kind="stable")
    ri, ei = rm.ravel()[order], ex.ravel()[order]
    bnd = np.flatnonzero(np.diff(ri)) + 1
    for grp in np.split(ei, bnd):
        if len(np.unique(grp)) != 1:
            errs.append("an input region was split across output regions")
            break

    # 2. output ids are input ids or 0
    extra = set(np.unique(ex)) - set(np.unique(rm)) - {0}
    if extra:
        errs.append(f"output invents region ids {sorted(extra)[:5]}")

    # 3. masking only grows
    if not np.all(ex[rm == 0] == 0):
        errs.append("a masked input pixel became unmasked")

    # 4. every output region is <= Nmax, or is an untouched input region
    ids, cnt = np.unique(ex[ex > 0], return_counts=True)
    isz = dict(zip(*np.unique(rm[rm > 0], return_counts=True)))
    for i, n in zip(ids, cnt):
        if n > c["max_region_size"] and isz.get(i) != n:
            errs.append(f"region {i} has {n} px > Nmax {c['max_region_size']} after merging")
            break

    # 5. counts agree with the manifest
    if len(ids) != c["regions_out"]:
        errs.append(f"regions_out {len(ids)} != {c['regions_out']}")
    if len(isz) != c["regions_in"]:
        errs.append(f"regions_in {len(isz)} != {c['regions_in']}")

    status = "OK  " if not errs else "FAIL"
    fail += bool(errs)
    print(f"{status} {c['name']:16s} {rm.shape[0]}x{rm.shape[1]}x{c['nbands']} "
          f"rmap={np.dtype(DT[rdt]).name:6s} layer={np.dtype(DT[ldt]).name} "
          f"maskin={c['masked_frac_in']:.3f} maskout={c['masked_frac_out']:.3f}")
    for e in errs:
        print("      !!", e)
print("\nFIXTURES OK" if not fail else f"\n{fail} CASES FAILED")
sys.exit(1 if fail else 0)
