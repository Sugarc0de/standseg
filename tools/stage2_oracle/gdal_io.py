import rasterio as rio
import os
import numpy as np


def read_region_map(region_map_filename):
    print("Reading region map from file: {}".format(region_map_filename))
    with rio.open(region_map_filename) as src:
        region_map = src.read(1)
        profile = src.profile
        print("The profile of the region map is: {}".format(profile))
        rows, columns = region_map.shape
        num_regions = np.max(region_map)
        print(
            "The region map has {} rows and {} columns, and {} regions in total".format(
                rows, columns, num_regions
            )
        )
        print("\n")
    return region_map, profile


def read_image(input_filename):
    print("Reading step 2 image from file: {}".format(input_filename))
    with rio.open(input_filename) as src:
        image = src.read()
        profile = src.profile
        print("The profile of the image is: {}".format(profile))
        bands, rows, columns = image.shape
        print(
            "The image has {} bands, {} rows and {} columns".format(
                bands, rows, columns
            )
        )
    return image


# TODO: convert data type when needed
def write_region_map(new_region_map, output_filename, profile):
    print("Writing new region map to file: {}".format(output_filename))
    with rio.open(output_filename, "w", **profile) as dst:
        dst.write(new_region_map, 1)
    print("Done writing new region map to file: {}".format(output_filename))
