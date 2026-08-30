from typing import List, Tuple
from random import randint
import math

# Instrumentation added for fixture generation. STATS["ties"] counts the tie
# branch in find_nearest_region -- the only place the algorithm consults the RNG
# -- and STATS["passes"] records per-pass merge/rejection counters, which give a
# Rust port the same ability to localise a divergence that myseg.log gives
# stage 1.
STATS = {"ties": 0, "cmp": 0, "passes": []}
# FAST swaps the O(n) list membership test in update_adjacent_regions for an
# O(1) set lookup. Provably equivalent: the test is pure membership, and the
# function only ORs bits, so neither the iteration order of coords nor the
# lookup structure can change its result.
FAST = True
# The two choices the original leaves undefined, both now explicit and both
# defaulting to the deterministic rule this project adopts.
#
# ORDER -- the original iterates a Python set, whose order is an implementation
# detail. "set" reproduces that, "asc"/"desc" sort by region id.
#
# ON_TIE -- the original calls an unseeded randint(0, 1) here, so the phase had
# no reproducible answer at all. "keep"/"take" are deterministic; "random"
# restores the original coin flip, and is the only setting that touches the RNG.
#
# Defaults are the canonical rule (ascending id, keep the incumbent), so
# importing this module gives a deterministic oracle with no RNG in the path.
# gen_fixtures.py varies them to measure how much the choice actually matters.
ORDER = "asc"
ON_TIE = "keep"


def _order(ids):
    if ORDER == "asc":
        return sorted(ids)
    if ORDER == "desc":
        return sorted(ids, reverse=True)
    return ids


def _tie_takes():
    """Should a candidate that ties the running best displace it?"""
    if ON_TIE == "keep":
        return False
    if ON_TIE == "take":
        return True
    return bool(randint(0, 1))


class Region:
    def __init__(self, id: int, coords: List[Tuple[int, int]]):
        self.coords = coords
        self.coord_set = set(coords)
        self.id = id
        self.num_invalid_pixel = 0
        # self.border_coords = self._get_border_coords()
        self.bounding_box = self._initialize_bounding_box()
        # self.adjacent_regions = {}
        self.centroids = []
        self.nearest_region_id = 0
        self.nearest_region_dist = math.inf

    def update_centroids(self, region):
        for i in range(len(self.centroids)):
            self.centroids[i] = (
                (self.size()) * self.centroids[i]
                + ((region.size()) * region.centroids[i])
            ) / (self.size() + region.size())

    def size(self):
        return len(self.coords)

    def extend_coords(self, region):
        # Before extending, update the new bounding box without reiterating through all the old coords
        self.bounding_box = (
            min(self.bounding_box[0], region.bounding_box[0]),
            min(self.bounding_box[1], region.bounding_box[1]),
            max(self.bounding_box[2], region.bounding_box[2]),
            max(self.bounding_box[3], region.bounding_box[3]),
        )
        self.coords.extend(region.coords)
        self.coord_set.update(region.coord_set)

    def merge_region(self, nearest_region, adjacent_info):
        adjacent_info = self.update_adjacent_regions(nearest_region, adjacent_info)
        self.update_centroids(nearest_region)
        self.num_invalid_pixel += nearest_region.num_invalid_pixel
        # update_centroids happens before extend_coords, otherwise the calculation will be wrong
        self.extend_coords(nearest_region)
        return adjacent_info

    def update_adjacent_regions(self, nearest_region, adjacent_info):
        for coord in self.coords:
            r, c = coord
            offset = [(-1, 0), (0, 1), (1, 0), (0, -1)]  # N, E, S, W
            for i in range(len(offset)):
                nr, nc = r + offset[i][0], c + offset[i][1]
                if (
                    nr < 0
                    or nc < 0
                    or nr >= len(adjacent_info)
                    or nc >= len(adjacent_info[0])
                ):
                    continue
                # Loop through each bits in adjacent_info[r][c], if the bit is 0, check
                # if the [nr][nc] is in nearest_region.coords.
                # If so, set the bit to 1. At the same time, set the corresponding bit in the adjacent_info[nr][nc] to be 1.
                # For example, if adjacent_info[r][c] = 0000, and adjacent_info[nr][nc] = 0000, and nr, nc is in nearest_region.coords,
                # if nr = r - 1, nc = c, then set adjacent_info[r][c] = 0001, and adjacent_info[nr][nc] = 0100.
                # If nr = r, nc = c + 1, then set adjacent_info[r][c] = 0010, and adjacent_info[nr][nc] = 1000.
                if adjacent_info[r][c] & (1 << i) == 0:
                    hit = ((nr, nc) in nearest_region.coord_set) if FAST \
                        else ((nr, nc) in nearest_region.coords)
                    if hit:
                        adjacent_info[r][c] |= 1 << i
                        adjacent_info[nr][nc] |= 1 << ((i + 2) % 4)
        return adjacent_info

    # Must call get_adjacent_regions() before the check the length
    def update_nearest_region_id_and_dist(self, region_id, dist):
        self.nearest_region_id = region_id
        self.nearest_region_dist = dist

    def update_nearest_region_dist(self, dist):
        self.nearest_region_dist = dist

    def _initialize_bounding_box(self):
        min_x, min_y, max_x, max_y = (
            float("inf"),
            float("inf"),
            float("-inf"),
            float("-inf"),
        )
        for coord in self.coords:
            y, x = coord  # r, c
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            max_x = max(max_x, x)
            max_y = max(max_y, y)
        return (min_x, min_y, max_x, max_y)

    def find_nearest_region(self, region_map, adjacent_info, regions):
        nearest_region_dist = float("inf")
        nearest_region_id = 0
        nearest_regions = set()
        # Loop through r and c of the bounding box
        for r in range(self.bounding_box[1], self.bounding_box[3] + 1):
            for c in range(self.bounding_box[0], self.bounding_box[2] + 1):
                # If the pixel is not in the region, skip it
                if region_map[r][c] != self.id:
                    continue

                # If there is all 1's in the 4 bit of adjacent_info at (r, c), current pixel is an internal pixel, skip it
                if adjacent_info[r][c] == 0b1111:
                    continue

                # Loop through the 4 neighbors (N, E, S, W) of the pixel, check if the corresponding adjacency_info bit is 0,
                # If so, add the neighbor in that direction to the nearest_regions set.
                for i in range(4):
                    # Extract the i-th bit from the value (right shift by i, then bitwise AND with 1)
                    bit = (adjacent_info[r][c] >> i) & 1

                    # Check if the bit is equal to 0 and compute the neighbor's position
                    if bit == 0:
                        if i == 0:  # North
                            nr, nc = r - 1, c
                        elif i == 1:  # East
                            nr, nc = r, c + 1
                        elif i == 2:  # South
                            nr, nc = r + 1, c
                        elif i == 3:  # West
                            nr, nc = r, c - 1
                        if (
                            nr == len(region_map)
                            or nc == len(region_map[0])
                            or nc < 0
                            or nr < 0
                        ):
                            continue
                        nearest_regions.add(region_map[nr][nc])

        if len(nearest_regions) > 1000000:
            print(
                "Region {} has {} way too many nearest regions".format(
                    self.id, len(nearest_regions)
                )
            )
        # Loop through the nearest regions, find the region with the smallest centroid distance.
        # If there are more than 1 region with the same distance, randomly pick one.
        for region_id in _order(nearest_regions):
            if region_id not in regions:
                # We dont keep track of regions that are entirely masked out
                continue
            region = regions[region_id]
            # For some stage two layers, for example the species probability map, it is expected that some bands might have 0's,
            # which means that the probability of that species occuring in that pixel is 0.
            # We just need to ensure pixels with all zeros have been filtered out.
            assert sum(region.centroids) != 0
            dist = self.get_spectral_dist(region)
            STATS["cmp"] += 1
            if math.isclose(dist, nearest_region_dist, rel_tol=1e-6):
                STATS["ties"] += 1
                if _tie_takes():
                    nearest_region_id = region_id
            elif dist < nearest_region_dist:
                nearest_region_dist = dist
                nearest_region_id = region_id
        return nearest_region_id, nearest_region_dist

    # Check if coord1 and coord2 are 4-way adjacent
    # TODO: move this to a pixel class
    @staticmethod
    def is_adjacent(coord1, coord2):
        x1, y1 = coord1
        x2, y2 = coord2
        return abs(x1 - x2) + abs(y1 - y2) == 1

    def get_spectral_dist(self, region):
        c1 = self.centroids
        c2 = region.centroids
        dist = 0
        for i in range(len(c1)):
            dist += (float(c1[i]) - float(c2[i])) ** 2
        return dist
