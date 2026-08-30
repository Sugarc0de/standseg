# This file implement the auxiliary steps for the Woodcock's segemnt algorithm. The file reads in the region map produced by the general pass.
# It also get the new input for step 2, three additional region sizes min_region_size, nviable, max_region_size, then perform the pass and output the new region map.
# No need to explicitly give it the mask because DN = 0 in the input region map represents the mask.

from collections import defaultdict
import math
import gdal_io
import numpy as np
from region import Region
from numba import jit


def crop_center_multiband(img, cropx, cropy):
    bands, y, x = img.shape
    startx = x // 2 - (cropx // 2)
    starty = y // 2 - (cropy // 2)
    return img[:, starty : starty + cropy, startx : startx + cropx]


def crop_center(img, cropx, cropy):
    y, x = img.shape
    startx = x // 2 - (cropx // 2)
    starty = y // 2 - (cropy // 2)
    return img[starty : starty + cropy, startx : startx + cropx]


def segment(
    min_region_size,
    max_region_size,
    output_filename,
    region_map_filename,
    image_filename,
    test=False,
):
    # Read in the input region map
    region_map, profile = gdal_io.read_region_map(region_map_filename)
    # Get the new input for step 2
    image = gdal_io.read_image(image_filename)

    # Crop the arrays
    # crop_size_x = 100
    # crop_size_y = 100
    # region_map = crop_center(region_map_temp, crop_size_x, crop_size_y)
    # image = crop_center_multiband(image_temp, crop_size_x, crop_size_y)
    # profile.update({"height": crop_size_y, "width": crop_size_x})

    print("new region map shape: {}".format(region_map.shape))
    print("image size is {}".format(image.shape))

    # Perform the pass
    new_region_map, num_passes = wind_up(
        region_map, image, min_region_size, max_region_size
    )
    if test:
        output_filename = output_filename
    else:
        output_filename = "{}.armap.{}".format(output_filename, num_passes)
    # Write out the new region map
    gdal_io.write_region_map(new_region_map, output_filename, profile)


def wind_up(region_map, image, min_region_size, max_region_size):
    old_num_regions = 0

    regions_temp = defaultdict(list)
    for r in range(len(region_map)):
        for c in range(len(region_map[0])):
            region_id = region_map[r][c]
            regions_temp[region_id].append((r, c))

    regions = {}
    for region_id, coords in regions_temp.items():
        if region_id == 0:
            continue
        regions[region_id] = Region(region_id, coords)

    print("About to make regions centroids")
    make_regions_centroids_and_remove_zero_centroid_region(image, regions, region_map)
    adjacency_info = create_adjacency_info_for_pixels(region_map)

    num_regions = len(regions)
    num_passes = 1
    while num_regions != old_num_regions:
        (region_map, regions, adjacency_info) = segment_per_pass(
            region_map, regions, adjacency_info, min_region_size, max_region_size
        )
        old_num_regions = num_regions
        log_num_passes(old_num_regions, regions, num_passes)
        num_passes += 1
        num_regions = len(regions)

    print("Pass {} completed without any merge.".format(num_passes))
    return region_map, num_passes


def segment_per_pass(
    region_map, regions, adjacency_info, min_region_size, max_region_size
):
    """
    Perform one pass of the nearest neighbor merging algorithm.
    """
    merged_regions = set()
    deleted_regions = set()
    total_merged = 0
    global_min_dist = math.inf
    from region import STATS
    C = {"considered": 0, "no_cand": 0, "busy": 0, "inf": 0,
         "over_max": 0, "not_mutual": 0, "merged": 0}

    # Make sure all the nearest neighbor distances are max float
    for region in regions.values():
        region.update_nearest_region_dist(math.inf)

    # One pass to update all the nearest neighbor regions and distances
    for region in regions.values():
        if region.size() >= min_region_size:
            continue
        assert sum(region.centroids) != 0
        nearest_region_id, nearest_region_dist = region.find_nearest_region(
            region_map, adjacency_info, regions
        )
        if nearest_region_id == 0:
            continue
        region.update_nearest_region_id_and_dist(nearest_region_id, nearest_region_dist)

        if nearest_region_dist < regions[nearest_region_id].nearest_region_dist:
            regions[nearest_region_id].update_nearest_region_id_and_dist(
                region.id, nearest_region_dist
            )

    for region in regions.values():
        if region.size() >= min_region_size:
            continue
        C["considered"] += 1
        if region.id in merged_regions or region.id in deleted_regions:
            C["busy"] += 1
            continue

        if len(merged_regions) % 200 == 0:
            print(
                "Merged {} regions out of {} total regions".format(
                    len(merged_regions), len(regions)
                )
            )
        nearest_region_id = region.nearest_region_id
        if nearest_region_id == 0:  # mask out region
            C["no_cand"] += 1
            continue
        if nearest_region_id in merged_regions or nearest_region_id in deleted_regions:
            C["busy"] += 1
            continue
        nearest_region_dist = region.nearest_region_dist
        if math.isclose(nearest_region_dist, math.inf):
            C["inf"] += 1
            continue
        nearest_region = regions[nearest_region_id]

        global_min_dist = min(global_min_dist, nearest_region_dist)

        if nearest_region.size() + region.size() > max_region_size:
            C["over_max"] += 1
            continue

        if not math.isclose(nearest_region_dist, nearest_region.nearest_region_dist):
            C["not_mutual"] += 1
            continue
        C["merged"] += 1
        # TODO: merge into the smaller number region (nice to have)
        adjacency_info = region.merge_region(nearest_region, adjacency_info)
        merged_regions.add(region.id)
        deleted_regions.add(nearest_region_id)
        total_merged += 1
    STATS.setdefault("passes", []).append(C)
    clean_up_regions(deleted_regions, regions)
    update_region_map(region_map, regions)
    print(
        f"The minimum nearest neighbor distance on this pass is {math.sqrt(global_min_dist)}"
    )
    return region_map, regions, adjacency_info


def log_num_passes(num_old_regions, new_regions, num_passes):
    num_new_regions = len(new_regions)
    print(
        "Pass %d: %d regions merged into %d regions"
        % (num_passes, num_old_regions, num_new_regions)
    )

    max_region_size = 0
    min_region_size = float("inf")

    for region in new_regions.values():
        min_region_size = min(min_region_size, region.size())
        if region.id != 0:
            max_region_size = max(max_region_size, region.size())

    print(
        "Pass %d: the smallest region size = %d, the largest region size = %d"
        % (num_passes, min_region_size, max_region_size)
    )


def clean_up_regions(deleted_regions, regions):
    for region_id in deleted_regions:
        if region_id in regions:
            del regions[region_id]


def update_region_map(region_map, regions):
    for region in regions.values():
        for coord in region.coords:
            region_map[coord[0]][coord[1]] = region.id


def make_regions_centroids_and_remove_zero_centroid_region(image, regions, region_map):
    majority_invalid_region_ids = []

    for region_id, region in regions.items():
        coords = np.array(region.coords)
        b_images = image[:, coords[:, 0], coords[:, 1]]
        region_centroids = b_images.mean(axis=1)
        # We only count the zero pixels if all the bands are zero
        all_bands_zero = np.all(b_images == 0, axis=0)
        region.num_zero_pixels = np.sum(all_bands_zero)

        if region.num_zero_pixels / len(region.coords) > 0.5:
            majority_invalid_region_ids.append(region_id)
        else:
            region.centroids = region_centroids.tolist()

    print(
        f"About to remove {len(majority_invalid_region_ids)} number of regions with more than half invalid pixels"
    )

    # FIXTURE-GENERATION FIX (see PLAN.md 13.2). With an empty list this built
    # `region_map[()] = 0`, which is the *whole array*, wiping every region id.
    # It never fired on the real runs because the structure/age/species layers
    # always contain some majority-nodata region, but it makes the stage silently
    # a no-op on any input with none. Guarded, not otherwise altered.
    if majority_invalid_region_ids:
        region_map[
            tuple(
                np.array(
                    [
                        coord
                        for region_id in majority_invalid_region_ids
                        for coord in regions[region_id].coords
                    ]
                ).T.tolist()
            )
        ] = 0

    for region_id in majority_invalid_region_ids:
        del regions[region_id]

    print("Finished removing regions with more than half invalid pixels")


# Create_adjacency_info_for_pixels(region_map) computes a 2D array where each element is 4 bits
# representing 4 directions (N, E, S, W). The value is 1 if the neighboring pixel in that direction
# is out of bounds or masked out or if it belongs to a pixel of the same region. Thus, if it is 0,
# the corresponding neighbor pixel is in bounds and belongs to a different region. Input is a 2D
# array of region ids.
@jit(nopython=True)
def create_adjacency_info_for_pixels(region_map):
    num_rows = len(region_map)
    num_cols = len(region_map[0])
    adjacency_info = np.ones((num_rows, num_cols), dtype=np.uint8) * 0b1111
    for row in range(num_rows):
        for col in range(num_cols):
            if region_map[row][col] == 0:
                continue
            if (
                row - 1 >= 0
                and region_map[row - 1][col] != region_map[row][col]
                and region_map[row - 1][col] != 0
            ):  # N
                adjacency_info[row][col] &= 0b1110
            if (  # E
                col + 1 < num_cols
                and region_map[row][col + 1] != region_map[row][col]
                and region_map[row][col + 1] != 0
            ):
                adjacency_info[row][col] &= 0b1101
            if (
                row + 1 < num_rows
                and region_map[row + 1][col] != region_map[row][col]
                and region_map[row + 1][col] != 0
            ):  # S
                adjacency_info[row][col] &= 0b1011
            if (
                col - 1 >= 0
                and region_map[row][col - 1] != region_map[row][col]
                and region_map[row][col - 1] != 0
            ):  # W
                adjacency_info[row][col] &= 0b0111
    print("Finished creating adjacency info for pixels")
    return adjacency_info
