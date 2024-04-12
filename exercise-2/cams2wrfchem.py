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

 - References:

   * van der Gon, H. D., Hendriks, C., Kuenen, J., Segers, A., and Visschedijk,
     A.: Description of Current Temporal Emission Patterns and Sensitivity of
     Predicted AQ for Temporal Emission Patterns, EU FP7 MACC deliverable
     report D_D-EMIS_1.3, 2011.

 - TODO (non-exhaustive):

   * Maybe libraries such as PySAL could be used to do the spatial attribution,
     instead of doing it by hand?

"""

import sys
if sys.version_info.major < 3:
    raise RuntimeError("Please use Python 3!")
from argparse import ArgumentParser, RawDescriptionHelpFormatter as Formatter
from itertools import product
from collections import namedtuple
from datetime import datetime, timedelta
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
        cell = Polygon([(x_bdy[i,j], y_bdy[i,j]),
                        (x_bdy[i,j+1], y_bdy[i,j+1]),
                        (x_bdy[i+1,j+1], y_bdy[i+1,j+1]),
                        (x_bdy[i+1,j], y_bdy[i+1,j])])
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

## Impure functions

def time_as_datetime(nc):
    """Return time variable of given netCDF file as array of datetimes."""
    time = nc.variables["time"]
    if time.getncattr("calendar") != "gregorian":
        raise NotImplementedError("Only supported calendar is gregorian.")
    units = time.getncattr("units")
    dt = units.split()[0]
    start = datetime.strptime(units, dt + " since %Y-%m-%d %H:%M:%S")
    dt = timedelta(**{dt: 1})
    return np.full(time.shape, start) + np.full(time.shape, dt)*time

def ll2xy_wrf(nc):
    """Return the ll2xy function corresponding to open WRF file."""
    proj = nc.getncattr("MAP_PROJ_CHAR")
    if proj == "Polar Stereographic":
        crs = pyproj.CRS.from_dict(dict(
            proj = "stere",
            lat_0 = str(nc.getncattr("POLE_LAT")),
            lon_0 = str(nc.getncattr("CEN_LON")),
            lat_ts = str(nc.getncattr("TRUELAT1")),
        ))
    else:
        raise NotImplementedError('Projection type "%s" not supported.' % proj)
    return pyproj.Transformer.from_crs(crs.geodetic_crs, crs).transform

def create_wrfchemi_file(nc, filepath):
    """Create WRF-Chem emission file with same grid as in nc, return handle."""
    new = netCDF4.Dataset(filepath, mode="w")
    new.createDimension("Time", 1)
    new.createDimension("DateStrLen", 19)
    new.createDimension("emissions_zdim_stag", 1)
    for dim in ("west_east", "south_north", "bottom_top"):
        size = [1, nc.dimensions[dim].size][dim != "bottom_top"]
        new.createDimension(dim, size)
        setattr(new, dim.upper()+"_GRID_DIMENSION", size+(dim!="bottom_top"))
    for attr in ("DX", "DY", "CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2",
                 "MOAD_CEN_LAT", "MAP_PROJ", "MMINLU"):
        setattr(new, attr, getattr(nc, attr))
    new.TITLE = "Created by " + __file__ + " on " + str(datetime.today())
    var = new.createVariable("Times", "c", ("Time", "DateStrLen"))
    var.datetime_format = "%Y-%m-%d_%H:%M:%S"
    return new

if __name__ == "__main__":

    ## Command-line arguments

    # Read and pre-process arguments
    parser = ArgumentParser(description=__doc__, formatter_class=Formatter)
    parser.add_argument("--start", required=True)
    parser.add_argument("--end", required=True)
    parser.add_argument("--period", choices=("hour", "day"), default="day")
    parser.add_argument("--CAMS-version", default="5.3")
    parser.add_argument("--ndomains", default="1", type=int)
    parser.add_argument("--dir-wrf-in", default="./")
    args = parser.parse_args()
    start = datetime.strptime(args.start, "%Y-%m-%d")
    end = datetime.strptime(args.end, "%Y-%m-%d")
    period = timedelta(**{args.period+"s": 1})

    # Quality controls on arguments
    if start.year != 2019 or end.year != 2019:
        # TODO I don't think this limit applies here
        raise NotImplementedError("Only year 2019-emissions are supported.")
    if start >= end:
        raise ValueError("Start date must be before end date.")
    if args.CAMS_version != "5.3":
        raise NotImplementedError("Only v5.3 of CAMS emissions is supported.")
    if args.ndomains < 0:
        raise ValueError("The number of domains must be a positive integer.")

    ## Hard-coded parameters

    dir_data = "/data/marelle/EMISSIONS/CAMS/v%s" % args.CAMS_version

    species_info = dict(
        CO = MozartMozaicSpecies("carbon-monoxide", 28),
        # CO_A = MozartMozaicSpecies("carbon-monoxide", 28),
        # NH3 = MozartMozaicSpecies("ammonia", 17),
        # SO2 = MozartMozaicSpecies("sulphur-dioxide", 64),
        # NO = MozartMozaicSpecies("nitrogen-oxides", 30),
        # NO2 = MozartMozaicSpecies("nitrogen-oxides", 30),
        # HCL = MozartMozaicSpecies(None, 36.46),
        # BIGALK = MozartMozaicSpecies("non-methane-vocs", None),
        # BIGENE = MozartMozaicSpecies("non-methane-vocs", None),
        # C2H4 = MozartMozaicSpecies("non-methane-vocs", None),
        # C2H5OH = MozartMozaicSpecies("non-methane-vocs", None),
        # C2H6 = MozartMozaicSpecies("non-methane-vocs", None),
        # C3H6 = MozartMozaicSpecies("non-methane-vocs", None),
        # C3H8 = MozartMozaicSpecies("non-methane-vocs", None),
        # CH2O = MozartMozaicSpecies("non-methane-vocs", None),
        # CH3CHO = MozartMozaicSpecies("non-methane-vocs", None),
        # CH3COCH3 = MozartMozaicSpecies("non-methane-vocs", None),
        # CH3OH = MozartMozaicSpecies("non-methane-vocs", None),
        # C2H2 = MozartMozaicSpecies("non-methane-vocs", None),
        # MEK = MozartMozaicSpecies("non-methane-vocs", None),
        # TOLUENE = MozartMozaicSpecies("non-methane-vocs", None),
        # BENZENE = MozartMozaicSpecies("non-methane-vocs", None),
        # XYLENE = MozartMozaicSpecies("non-methane-vocs", None),
        # ORGJ = MozartMozaicSpecies("organic-carbon", None),
        # ECJ = MozartMozaicSpecies("black-carbon", None),
        # SO4J = MozartMozaicSpecies("sulphur-dioxide", 64),
    )
    species_mozmoz = sorted(species_info.keys())
    species_cams = sorted(set(sp.name_cams for _, sp in species_info.items()
                              if sp.name_cams is not None))
    sectors_cams = ("ene", "ind", "res", "tro", "tnr", "ags", "agl", "swd",
                    "slv", "fef", "shp")

    # Daily factors (Monday to Sunday) from TNO-MACC (van der Gon et al. 2011)
    daily_factors = dict(ene = (1.06, 1.06, 1.06, 1.06, 1.06, 0.85, 0.85),
                         ind = (1.08, 1.08, 1.08, 1.08, 1.08, 0.80, 0.80),
                         res = (1.08, 1.08, 1.08, 1.08, 1.08, 0.80, 0.80),
                         tro = (1.02, 1.06, 1.08, 1.10, 1.14, 0.81, 0.79),
                         ags = tuple([1]*7), swd = tuple([1]*7),
                         slv = tuple([1]*7), fef = tuple([1]*7),
                         shp = tuple([1]*7))
    daily_factors["tnr"] = daily_factors["tro"]
    daily_factors["agl"] = daily_factors["ags"]

    # Hourly factors (0h to 23h) from TNO-MACC (van der Gon et al. 2011)
    hourly_factors = dict(
        ene = [0.88, 0.79, 0.72, 0.72, 0.71, 0.74, 0.80, 0.92,
               1.08, 1.19, 1.22, 1.21, 1.21, 1.17, 1.15, 1.14,
               1.13, 1.10, 1.07, 1.04, 1.02, 1.02, 1.01, 0.96],
        ind = [0.75, 0.75, 0.75, 0.78, 0.82, 0.88, 0.95, 1.02,
               1.09, 1.16, 1.22, 1.28, 1.30, 1.22, 1.24, 1.25,
               1.16, 1.08, 1.01, 0.95, 0.90, 0.85, 0.81, 0.78],
        res = [0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.50, 1.20,
               1.50, 1.60, 1.60, 1.40, 1.20, 1.10, 1.10, 1.00,
               1.00, 1.00, 1.10, 1.40, 1.50, 1.40, 1.40, 1.00],
        tro = [0.44, 0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86,
               1.84, 1.86, 1.41, 1.24, 1.20, 1.32, 1.44, 1.45,
               1.59, 2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61],
        ags = [0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75,
               0.90, 1.10, 1.25, 1.45, 1.60, 1.80, 1.75, 1.70,
               1.55, 1.35, 1.10, 0.90, 0.75, 0.65, 0.60, 0.60],
        swd = tuple([1]*24), slv = tuple([1]*24),
        fef = tuple([1]*24), shp = tuple([1]*24)
    )
    hourly_factors["tnr"] = hourly_factors["tro"]
    hourly_factors["agl"] = hourly_factors["ags"]
