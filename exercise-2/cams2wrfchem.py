"""Prepare WRF-chem emission input files using CAMS emission inventory.

Copyright (c) 2024 Regional-Modeling-LATMOS-IGE.

This software is released under the terms of the Expat (aka MIT) license:

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

-------------------------------------------------------------------------------

Notes:

 - This script is not operational yet -- this is work in progress!

 - Originally (things may have evolved since), this script was an adaptation of
   a Matlab(r) script written by Louis Marelle, timestamped 30th of May 2023
   and named "prep_anthro_emissions_mozartmosaic_cams53.m".

 - TODO (non-exhaustive):

   * Maybe libraries such as PySAL could be used to do the spatial attribution,
     instead of doing it by hand?

"""

import sys
if sys.version_info.major < 3:
    raise RuntimeError("Please use Python 3!")
from itertools import product
from collections import namedtuple
import numpy as np
try:
    from shapely import Polygon
except ImportError:
    from shapely.geometry import Polygon

## Data structures

MozartMozaicSpecies = namedtuple("MozartMozaicSpecies", "name_cams molmass")

CellWeight = namedtuple("CellWeight", "i j weight")

## Pure functions

def unique_index(vector, value):
    """Return index of value in vector, making sure it is unique."""
    if len(vector.shape) != 1:
        raise ValueError("Expecting a vector.")
    idx = np.where(vector == value)
    if len(idx[0]) != 1:
        raise ValueError("Could not determine unique index.")
    return int(idx[0][0])

def meshgrid(x, y):
    """Return meshgrid of x,y, using ij indexing convention."""
    return np.meshgrid(x, y, indexing="ij")

def check_shape_of_arrays(*arrays):
    """Return shape of arrays (raise ValueError if shapes differ)."""
    if len(arrays) == 0:
        return None
    shape = arrays[0].shape
    for array in arrays[1:]:
        if array.shape != shape:
            raise ValueError("Input arrays do not have the same shape.")
    return shape

def polygon_from_bounds(xmin, xmax, ymin, ymax):
    """Return rectangular polygon corresponding to given bounds."""
    return Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])

def polygon_from_center(xc, yc, dx, dy):
    """Return rectangular polygon corresponding to given center and size."""
    return polygon_from_bounds(xc-dx/2, xc+dx/2, yc-dy/2, yc+dy/2)

def surrounding_cells(i, j, nrows, ncols):
    """Return the list of indices of the cells surrounding [i,j]."""
    return [(x,y) for x, y in product(range(i-1,i+2), range(j-1,j+2))
            if x >= 0 and x < nrows and y >= 0 and y < ncols and
            (x != i or y != j)]

def calc_cellweights(x, y, x_bdy, y_bdy, polygon):
    """Return cell weights corresponding to polygon in given grid.

    For example, if this function returns [(1, 2, 0.25), (3, 4, 0.75)], it
    means that 25% of given polygon intersects with cell (1,2) of given grid
    and that 75% of given polygon intersects with cell (3,4) of given grid.

    """
    nrows, ncols = check_shape_of_arrays(x, y)
    check_shape_of_arrays(x_bdy, y_bdy)
    area_polygon = polygon.area
    xc, yc = polygon.centroid.coords.xy
    # The idea here is to find at least one grid cell that intersects with
    # given polygon. The implementation is straight-forward but not very
    # efficient. TODO improve this search
    dist2 = (x - xc)**2 + (y - yc)**2
    idx = np.where(dist2 == dist2.min())
    out, queue = [], [(idx[0][i], idx[1][i]) for i in range(len(idx[0]))]
    current = 0
    while current < len(queue):
        i, j = queue[current]
        xx, yy = x_bdy[i:i+2,j:j+2], y_bdy[i:i+2,j:j+2]
        cell = polygon_from_bounds(xx.min(), xx.max(), yy.min(), yy.max())
        area_intersection = cell.intersection(polygon).area
        if area_intersection > 0:
            out.append(CellWeight(i, j, area_intersection/area_polygon))
            surrounding = surrounding_cells(i, j, nrows, ncols)
            for cell in [cell for cell in surrounding if cell not in queue]:
                queue.append(cell)
        current += 1
    if abs(sum(cell.weight for cell in out)-1) > 1e-10:
        raise ValueError("Total weight is not equal to 1.")
    return out

def calc_mapping(x_from, y_from, x_bdy_from, y_bdy_from,
                 x_to, y_to, dx_to, dy_to):
    """Calculate the mapping between given grids.

    This function returns a mapping array (same size as x_to and y_to) where
    map[i][j] is a list of the cell weights corresponding to the projection of
    grid cell [i,j] of target grid onto the original grid.

    """
    check_shape_of_arrays(x_from, y_from)
    check_shape_of_arrays(x_bdy_from, y_bdy_from)
    nrows, ncols = check_shape_of_arrays(x_to, y_to)
    mapping = np.empty([nrows, ncols], dtype=list)
    for i, j in product(range(nrows), range(ncols)):
        polygon = polygon_from_center(x_to[i,j], y_to[i,j], dx_to, dy_to)
        mapping[i][j] = calc_cellweights(
            x_from, y_from, x_bdy_from, y_bdy_from, polygon)
    return mapping
