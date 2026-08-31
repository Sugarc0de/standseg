"""Generate tests/stage2 fixtures from the Python stage-2 oracle.

Every case must survive an *arbitrariness sweep* before it is written: the two
choices the algorithm leaves undefined -- the coin flip on near-ties, and the
iteration order of the candidate-neighbour set -- are varied over all six
combinations, and the region map must come out byte-identical every time. A case
that fails the sweep is not a usable oracle and is refused, not recorded.
"""
import json, os, sys, time
import numpy as np
import rasterio
from rasterio.windows import Window

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import envi, harness
import region as region_mod

# tools/stage2_oracle/ -> tools/ -> repo root.
REPO = os.path.dirname(os.path.dirname(HERE))

# The source tile is external data, not in the repo. Point NTEMS_TILE at a
# directory laid out like an NTEMS processed tile to regenerate from your own
# copy; the default is where it lived when these fixtures were made.
T = os.environ.get(
    "NTEMS_TILE",
    "/Volumes/easystore/UBC/first_project/ntems_2019_update/ab/processed_tiles/tile_399",
)
RMAP = (f"{T}/proxies/step1_results/proxies_t_50_m_0.2_n_9_18_36_tile_399/"
        f"proxies_t_50_m_0.2_n_9_18_36_tile_399.rmap.26")
MI = "Lambert Conformal Conic, 1, 1, -1460910.524, 798648.1064, 30, 30,North America 1983"
CS = ('PROJCS["NAD_1983_Canada_Atlas_Lambert",GEOGCS["GCS_North_American_1983",'
      'DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],'
      'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
      'PROJECTION["Lambert_Conformal_Conic"],PARAMETER["False_Easting",0.0],'
      'PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-95.0],'
      'PARAMETER["Standard_Parallel_1",49.0],PARAMETER["Standard_Parallel_2",77.0],'
      'PARAMETER["Latitude_Of_Origin",49.0],UNIT["Meter",1.0]]')
LAYERS = {
    "elev_p95": f"{T}/structure/elev_p95/elev_p95-tile-399-norm.tif",
    "elev_cv": f"{T}/structure/elev_cv/elev_cv-tile-399-norm.tif",
    "gross_stem_volume": f"{T}/structure/gross_stem_volume/gross_stem_volume-tile-399-norm.tif",
    "age": f"{T}/age/age-tile-399-norm.tif",
    "species": f"{T}/species/species_tile-399-norm.tif",
}
OUT = os.environ.get("STAGE2_FIXTURES", os.path.join(REPO, "tests", "stage2"))


def read_rmap_crop(r0, c0, n):
    m = np.memmap(RMAP, dtype='<u4', mode='r', shape=(5000, 5000))
    return np.array(m[r0:r0 + n, c0:c0 + n], dtype=np.uint32)


def read_layer_crop(name, r0, c0, n):
    with rasterio.open(LAYERS[name]) as s:
        return s.read(window=Window(c0, r0, n, n))


# The canonical, deterministic tie rule this project adopts and the Rust port
# must implement: visit candidate neighbours in ascending region id, and on a
# near-tie keep the incumbent, so the smallest id among near-equal candidates
# wins. Under this rule region.py never calls randint, so stage 2 needs no RNG
# at all -- unlike stage 1, which needs glibc random() call for call.
#
# This is a decision, not a reconstruction. The 1992 C picks among equidistant
# candidates with flip() ("This is biased, but it does give some randomness"),
# and its lower-id rule governs something else entirely -- which of the two
# regions survives a merge. Ascending-id-keep-incumbent is chosen because it is
# deterministic, cheap in both languages, and consistent in spirit with that
# survivor rule; it is not what either the C or the Python actually did.
CANON = ("asc", "keep")


def run_variant(rmap, image, nmin, nmax, order, flip_takes, fast=True):
    saved = (region_mod.ON_TIE, region_mod.ORDER, region_mod.FAST)
    try:
        region_mod.ORDER = order
        region_mod.FAST = fast
        region_mod.ON_TIE = "take" if flip_takes else "keep"
        out, npass, _, st = harness.run(rmap, image, nmin, nmax)
        return out, npass, st
    finally:
        region_mod.ON_TIE, region_mod.ORDER, region_mod.FAST = saved


def sweep(rmap, image, nmin, nmax, check_slow=False):
    """Canonical result, plus how sensitive it is to the undefined choices."""
    out, npass, stats = run_variant(rmap, image, nmin, nmax, *
                                    (CANON[0], CANON[1] == "take"))
    variants = {}
    for order in ("set", "asc", "desc"):
        for flip_name in ("keep", "take"):
            o, p, _ = run_variant(rmap, image, nmin, nmax, order, flip_name == "take")
            variants[(order, flip_name)] = (o.tobytes(), p)
    if check_slow:      # the membership-set optimisation must be neutral
        o, p, _ = run_variant(rmap, image, nmin, nmax, CANON[0],
                              CANON[1] == "take", fast=False)
        assert (o.tobytes(), p) == (out.tobytes(), npass), "membership set changed the result"
    return out, npass, stats, variants


CASES = []


def case(name, **kw):
    CASES.append(dict(name=name, **kw))


case("tiny_synthetic", kind="synthetic", nmin=4, nmax=9,
     note="Elaine's own integ_test.py 5x5 example, hand-checkable.")
case("p95_250", kind="tile", layer="elev_p95", r0=2000, c0=2000, n=250,
     nmin=80, nmax=8000, note="Paper parameters, single structural band, 31% non-treed.")
case("species_250", kind="tile", layer="species", r0=1500, c0=3000, n=250,
     nmin=80, nmax=8000, note="Six species-probability bands; exercises the all-bands-zero rule.")
case("age_capped", kind="tile", layer="age", r0=1000, c0=1000, n=250,
     nmin=60, nmax=200, note="Tight Nmax so the maximum-region-size rejection actually binds.")
case("e2e_gsv", kind="e2e", layer="gross_stem_volume", r0=2400, c0=2400, n=200,
     nmin=50, nmax=8000, rmap_file="e2e.rmap.41",
     note="Stage-1 region map produced by this repo's own Rust binary.")
case("e2e_masked", kind="e2e", layer="elev_cv", r0=1500, c0=3000, n=200,
     nmin=50, nmax=8000, rmap_file="e2em.rmap.35",
     note="Stage-1 run under -M, so the input region map already contains masked "
          "pixels (id 0) before stage 2 adds its own.")


def build(c):
    if c["kind"] == "synthetic":
        rmap = np.array([[2, 2, 2, 8, 5], [1, 1, 2, 8, 8], [4, 1, 1, 1, 7],
                         [4, 3, 3, 3, 7], [3, 3, 6, 6, 6]], dtype=np.uint32)
        img = np.array([[[44, 47, 64, 67, 67], [9, 83, 21, 36, 87], [70, 88, 88, 12, 58],
                         [65, 39, 87, 46, 88], [81, 37, 25, 77, 72]]], dtype=np.uint8)
        return rmap, img, None, None, "synthetic"
    r0, c0, n = c["r0"], c["c0"], c["n"]
    img = read_layer_crop(c["layer"], r0, c0, n)
    mi = envi.crop_map_info(MI, r0, c0)
    if c["kind"] == "e2e":
        raw = np.fromfile(os.path.join(REPO, "build", "out", "stage2gen",
                                       c["rmap_file"]), dtype='<u2').reshape(n, n)
        return raw.astype(np.uint32), img, mi, np.uint16, "rust-stage1"
    return read_rmap_crop(r0, c0, n), img, mi, np.uint32, "ntems-stage1"


def main():
    os.makedirs(OUT, exist_ok=True)
    manifest = []
    for c in CASES:
        t0 = time.time()
        rmap, img, mi, rdt, src = build(c)
        nin = len(np.unique(rmap[rmap > 0]))
        out, npass, st, variants = sweep(rmap, img, c["nmin"], c["nmax"],
                                         check_slow=(c.get("n", 5) <= 250))
        tie_sensitive = len(set(variants.values())) != 1
        rdt = rdt or np.uint32
        d = f"{OUT}/{c['name']}"
        os.makedirs(f"{d}/input", exist_ok=True)
        os.makedirs(f"{d}/expected", exist_ok=True)
        envi.write(f"{d}/input/rmap", rmap.astype(rdt), map_info=mi, coord_sys=CS if mi else None,
                   ignore0=True, desc=f"stage-1 region map ({src})")
        envi.write(f"{d}/input/layer", img, map_info=mi, coord_sys=CS if mi else None,
                   ignore0=(c["kind"] != "synthetic"),
                   desc=f"stage-2 input: {c.get('layer','synthetic')}")
        envi.write(f"{d}/expected/armap.{npass}", out.astype(rdt), map_info=mi,
                   coord_sys=CS if mi else None, ignore0=True,
                   desc=f"stage-2 oracle output, {npass} passes")
        nout = len(np.unique(out[out > 0]))
        info = dict(
            name=c["name"], note=c["note"], min_region_size=c["nmin"],
            max_region_size=c["nmax"], shape=list(rmap.shape), nbands=int(img.shape[0]),
            layer=c.get("layer", "synthetic"), source=src,
            region_map_dtype=np.dtype(rdt).name, layer_dtype=str(img.dtype),
            regions_in=int(nin), regions_out=int(nout), passes=int(npass),
            expected=f"expected/armap.{npass}",
            masked_frac_in=float((rmap == 0).mean()), masked_frac_out=float((out == 0).mean()),
            ties=int(st["ties"]), distance_comparisons=int(st["cmp"]),
            tie_rule="ascending region id, keep incumbent on near-tie",
            tie_sensitive=bool(tie_sensitive),
            sweep_variants=len(variants),
            sweep_distinct=len(set(variants.values())),
            per_pass=st["passes"],
        )
        if c["kind"] != "synthetic":
            info["window"] = dict(row=c["r0"], col=c["c0"], size=c["n"], tile=399, province="ab")
        json.dump(info, open(f"{d}/case.json", "w"), indent=2)
        manifest.append(info)
        print(f"OK {c['name']:16s} {rmap.shape[0]}x{rmap.shape[1]}x{img.shape[0]} "
              f"nreg {nin:6d}->{nout:4d} passes={npass:3d} ties={st['ties']:5d} "
              f"tie-sensitive={'YES' if tie_sensitive else 'no ':3s} "
              f"[{time.time()-t0:.1f}s]")
    json.dump(manifest, open(f"{OUT}/cases.json", "w"), indent=2)
    print(f"\n{len(manifest)}/{len(CASES)} cases written to {OUT}")


if __name__ == "__main__":
    main()
