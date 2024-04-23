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

 - This script is not fully implemented nor tested yet -- work in progress!

 - Originally (things may have evolved since), this script was an adaptation of
   a Matlab(r) script written by Louis Marelle, timestamped 30th of May 2023
   and named "prep_anthro_emissions_mozartmosaic_cams53.m".

 - This script currently only supports the creation of emission files for the
   Mozart Mozaic chemical mechanism.

 - In this script, x and y are generally the arrays of coordinates for the
   centers of grid cells, while x_bdy and y_bdy are the arrays of coordinates
   for the boundaries of grid cells. If x and y have shape (n,m), then x_bdy
   and y_bdy have shape (n+1,m+1).

 - References:

   * Alexander, B.; Park, R. J.; Jacob, D. J.; Gong, S.: Transition
     metal-catalyzed oxidation of atmospheric sulfur: Global implications for
     the sulfur budget. Journal of Geophysical Research, 114, D02309, 2009.

   * Chin, M.; Rood, R. B.; Lin, S.-J.; Muller, J.-F., Thompson,
     A. M. Atmospheric sulfur cycle simulated in the global model GOCART: Model
     description and global properties. Journal of Geophysical Research, 105
     (D20), 24671-24687, 2009, 2000.

   * Finlayson-Pitts, B. J.; Pitts, J. N. Jr. Chemistry of the Upper and Lower
     Atmosphere: Theory, Experiments, and Applications. Elsevier Inc, 2000.

   * van der Gon, H. D.; Hendriks, C.; Kuenen, J.; Segers, A.; and Visschedijk,
     A.: Description of Current Temporal Emission Patterns and Sensitivity of
     Predicted AQ for Temporal Emission Patterns, EU FP7 MACC deliverable
     report D_D-EMIS_1.3, 2011.

 - TODO (non-exhaustive):

   * Maybe libraries such as PySAL could be used to do the spatial attribution,
     instead of doing it by hand?

   * The fraction of sulfate in anthropogenic sulfur emissions in Min et
     al. (2000) is actually "5% for Europe and 3% for elsewhere". This script
     uses 3% in all situations. Does this need to be adapted ?

   * Fix attribution of cell weights when current method gives a total weight
     that is < 1 (use brute-force in this case?)

   * For the NO/NOx ratio in anthropogenic emissions, Louis mentions Kaynak et
     al. 2009 as a reference, along with Finlayson-Pitts & Pitts (2000). Should
     we add this reference if the information is already avalaible in
     Finlayson-Pitts & Pitts (2000)?

   * In Louis' original script, there is a variable defining the ratio NO/NOx
     specifically for ship emissions, but this variable does not seem to be
     used anywhere, although CAMS emissions do have a shipping sector. Should
     we implement the use of this ratio?

   * In Louis' original script, there was a check that stopped the script if
     selected start and end dates were not in 2019. It is also implemented
     here, but is it necessary? I am guessing that the "2019" in CAMS'
     filenames probably refer to the year of the method used to
     # create the files. The files I looked at actually contained emissions for
     years 2000 through 2013.

   * It would probably be useful to write somewhere in this script what the
     CAMS emission sector abbreviations stand for.

   * For WRF species with no CAMS equivalent (eg. HCL), we should create the
     emissions and set them to zero. For now the emissions are simply not
     created.

"""

import sys
if sys.version_info.major < 3:
    raise RuntimeError("Please use Python 3!")
from os.path import join
from argparse import ArgumentParser, RawDescriptionHelpFormatter as Formatter
from itertools import product
from collections import namedtuple
import json
from datetime import datetime, timedelta
import numpy as np
import netCDF4
import pyproj
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

def shape_of_arrays(*arrays):
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
    means that 25% of given polygon intersects with cell [1,2] of given grid
    and that 75% of given polygon intersects with cell [3,4] of given grid.

    """
    nrows, ncols = shape_of_arrays(x, y)
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
            for cell in surrounding_cells(i, j, nrows, ncols):
                if cell not in queue:
                    queue.append(cell)
        current += 1
    weight = sum(cell.weight for cell in out)
    if abs(weight-1) > 1e-10:
        # TODO fix cases when total weight is < 1
        print("Total weight is not equal to 1 (%.2f)." % weight)
    return out

def calc_mapping(x_from, y_from, x_bdy_from, y_bdy_from,
                 x_to, y_to, dx_to, dy_to):
    """Calculate the mapping between given grids.

    This function returns a mapping array (same shape as x_to and y_to) where
    mapping[i][j] is a list of the cell weights corresponding to the projection
    of grid cell [i,j] of target grid onto the original grid.

    """
    nrows, ncols = shape_of_arrays(x_to, y_to)
    mapping = np.empty([nrows, ncols], dtype=list)
    for i in range(nrows):
        print("Spatial mapping progress = %.1f%%" % (100*i/nrows))
        for j in range(ncols):
            polygon = polygon_from_center(x_to[i,j], y_to[i,j], dx_to, dy_to)
            mapping[i][j] = calc_cellweights(
                x_from, y_from, x_bdy_from, y_bdy_from, polygon)
    return mapping

def remap(data, mapping):
    """Remap given data with given mapping."""
    nrows, ncols = mapping.shape
    out = np.zeros([nrows, ncols])
    for i, j in product(range(nrows), range(ncols)):
        for cw in mapping[i][j]:
            out[i,j] += data[cw.i,cw.j] * cw.weight
    return out

## Impure functions which do not use the parent environment

def jsonifile_mapping(mapping, filepath):
    """Write mapping to given file in JSON format."""
    cell2json = lambda c: ("%d" % c.i, "%d" % c.j, "%.10e" % c.weight)
    towrite = dict(nrows = mapping.shape[0], ncols = mapping.shape[1],
                   values = [[] if v is None else list(map(cell2json, v))
                             for v in mapping.flatten()])
    with open(filepath, mode="w") as f:
        json.dump(towrite, f)

def unjsonifile_mapping(filepath):
    """Return mapping read from given file (written in JSON format)."""
    with open(filepath, mode="r") as f:
        data = json.load(f)
    nrows, ncols = data["nrows"], data["ncols"]
    mapping = np.empty([nrows, ncols], dtype=list)
    counter = -1
    for i, j in product(range(nrows), range(ncols)):
        counter += 1
        mapping[i,j] = [CellWeight(int(d[0]), int(d[1]), float(d[2]))
                        for d in data["values"][counter]]
    return mapping

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
    """Return the ll2xy function corresponding to given WRF file."""
    proj = nc.getncattr("MAP_PROJ_CHAR")
    if proj == "Polar Stereographic":
        crs = pyproj.CRS.from_dict(dict(
            proj = "stere",
            lat_0 = str(nc.getncattr("POLE_LAT")),
            lon_0 = str(nc.getncattr("CEN_LON")),
            lat_ts = str(nc.getncattr("TRUELAT1")),
        ))
    elif proj == "Lambert Conformal":
        crs = pyproj.CRS.from_dict(dict(
            proj = "lcc",
            lat_0 = str(nc.getncattr("CEN_LAT")),
            lon_0 = str(nc.getncattr("CEN_LON")),
            lat_1 = str(nc.getncattr("TRUELAT1")),
            lat_2 = str(nc.getncattr("TRUELAT2")),
        ))
    else:
        raise NotImplementedError('Projection type "%s" not supported.' % proj)
    return pyproj.Transformer.from_crs(crs.geodetic_crs, crs).transform

def get_wrf_emissions_file(domain, time, ncs, nc_grid):
    """Return handle to opened WRF emissions file (create it if needed).

    Parameters "ncs" is a dictionary containing the handles of the WRF emission
    files that have already been created (the keys are the filepaths). If a new
    file is created, its handle is added to this dictionary.

    When creating a new file, the relevant grid information is read from
    nc_grid.

    """
    time_f = time.strftime("%Y-%m-%d_%H:%M:%S")
    filepath = "wrfchemi_d%s_%s" % (domain, time_f)
    try:
        nc = ncs[filepath]
    except KeyError:
        nc = netCDF4.Dataset(filepath, mode="w")
        nc.createDimension("Time", 1)
        nc.createDimension("DateStrLen", 19)
        nc.createDimension("emissions_zdim_stag", 1)
        for dim in ("west_east", "south_north", "bottom_top"):
            size = [nc_grid.dimensions[dim].size, 1][dim == "bottom_top"]
            nc.createDimension(dim, size)
            setattr(nc, dim.upper()+"_GRID_DIMENSION", size+(dim!="bottom_top"))
        for attr in ("DX", "DY", "CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2",
                     "MOAD_CEN_LAT", "MAP_PROJ", "MMINLU"):
            setattr(nc, attr, getattr(nc_grid, attr))
        nc.TITLE = "Created by " + __file__ + " on " + str(datetime.today())
        var = nc.createVariable("Times", "c", ("Time", "DateStrLen"))
        var.datetime_format = "%Y-%m-%d_%H:%M:%S"
        var[0,:] = time_f
        ncs[filepath] = nc
    return nc

def get_wrf_emissions_variable(nc, species, units):
    """Get emissions variable for species from file (create it if needed)."""
    varname = "E_" + species
    try:
        var = nc.variables[varname]
    except KeyError:
        dims = ("Time", "emissions_zdim_stag", "south_north", "west_east")
        var = nc.createVariable(varname, np.float32, dims)
        var.FieldType = 104
        var.MemoryOrder = "XYZ"
        var.description = "EMISSIONS"
        var.stagger = "Z"
        var.units = units
        var[:] = 0
    else:
        if var.units != units:
            raise ValueError("Given units do not match those on record.")
    return var

def close_ncfiles_in_dict(ncs):
    """Close of the files which handles are stored in given dictionary."""
    for nc in list(ncs.keys()):
        ncs.pop(nc).close()

## Here starts the "script part" of this script

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
    parser.add_argument("--dir-em-in",
                        default="/bettik/PROJECTS/pr-regionalchem/laperer/WRFCHEM_INPUT/CAMS_EMISSIONS/v5.3")
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
    if args.ndomains <= 0:
        raise ValueError("The number of domains must be a positive integer.")

    ## Hard-coded parameters

    species_info = dict(
        CO = MozartMozaicSpecies("carbon-monoxide", 28),
        # CO_A = MozartMozaicSpecies("carbon-monoxide", 28),
        # NH3 = MozartMozaicSpecies("ammonia", 17),
        SO2 = MozartMozaicSpecies("sulphur-dioxide", 64),
        NO = MozartMozaicSpecies("nitrogen-oxides", 30),
        NO2 = MozartMozaicSpecies("nitrogen-oxides", 30),
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
        ORGJ = MozartMozaicSpecies("organic-carbon", None),
        ECJ = MozartMozaicSpecies("black-carbon", None),
        SO4J = MozartMozaicSpecies("sulphur-dioxide", 64),
    )
    all_species_wrf = sorted(species_info.keys())
    all_species_cams = sorted(set(sp.name_cams for _, sp in species_info.items()
                                  if sp.name_cams is not None))
    all_sectors_cams = ("ene", "ind", "res", "tro", "tnr", "ags", "agl", "swd",
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

    # Fractions of anthropogenic NOx emissions released as NO and NO2
    # (Finlayson-Pitts & Pitts, 2000)
    frac_NO_anthropo = 0.9
    frac_NO2_anthropo = 1 - frac_NO_anthropo

    # Fractions of anthropogenic sulfur emissions released as SO2 and sulfate
    # (Alexander et al. 2009, citing Chin et al. 2000)
    frac_SO4_anthropo = 0.03
    frac_SO2_anthropo = 1 - frac_SO4_anthropo

    ## Prepare input files and grid information (open all the connections)

    ncs_cams = dict()
    for i, species in enumerate(all_species_cams):
        filename = "CAMS-GLOB-ANT_v%s_%s_2019.nc" % (args.CAMS_version, species)
        filepath = join(args.dir_em_in, filename)
        nc = ncs_cams[species] = netCDF4.Dataset(filepath)
        if nc.getncattr("projection") != "latlon":
            raise ValueError("Expecting lat-lon grid.")
        if i == 0:
            lat_cams = nc["lat"]
            lon_cams = nc["lon"]
            times_cams = time_as_datetime(nc)
            if any(t.strftime("%d%H%M%S") != "01000000" for t in times_cams):
                raise ValueError("Expecting monthly data in input files.")
        else:
            lat_ok = np.all(nc["lat"][:] == lat_cams[:])
            lon_ok = np.all(nc["lon"][:] == lon_cams[:])
            time_ok = np.all(time_as_datetime(nc) == times_cams)
            if not all([lat_ok, lon_ok, time_ok]):
                raise ValueError("Inconsistent grid data across input files.")
    times_cams = np.array([t.strftime("%Y%m") for t in times_cams])
    lat_cams_bdy, lon_cams_bdy = meshgrid(
        [-90] + list((lat_cams[1:] + lat_cams[:-1]) / 2) + [90],
        [-180] + list((lon_cams[1:] + lon_cams[:-1]) / 2) + [180],
    )
    lat_cams, lon_cams = meshgrid(lat_cams, lon_cams)

    # Process emissions
    for domain in range(args.ndomains):

        # Get basic information about WRF grid
        domain = ("%2d" % (domain+1)).replace(" ", "0")
        filename = "wrfinput_d%s" % domain
        nc_grid = netCDF4.Dataset(join(args.dir_wrf_in, filename))
        dx_wrf, dy_wrf = nc_grid.getncattr("DX"), nc_grid.getncattr("DY")
        ll2xy = ll2xy_wrf(nc_grid)
        x_wrf, y_wrf = ll2xy(nc_grid["XLONG"][0,:,:], nc_grid["XLAT"][0,:,:])
        x_cams, y_cams = ll2xy(lon_cams, lat_cams)
        x_cams_bdy, y_cams_bdy = ll2xy(lon_cams_bdy, lat_cams_bdy)

        # Calculate the spatial mapping (this is the time consuming operation)
        filepath = __file__[:-3] + "_mapping_d" + domain + ".json"
        try:
            mapping = unjsonifile_mapping(filepath)
        except FileNotFoundError:
            print("Calculating spatial mapping...")
            mapping = calc_mapping(x_cams, y_cams, x_cams_bdy, y_cams_bdy,
                                   x_wrf, y_wrf, dx_wrf, dy_wrf)
            jsonifile_mapping(mapping, filepath)
            print("Spatial mapping saved in %s." % filepath)
        else:
            print("Spatial mapping read from %s." % filepath)

        hourly_factors_cache = dict()
        def get_hourly_factors(hour, sector):
            """Return hourly factors for wrf grid for given hour and sector.

            We cache the results of this function for efficiency.

            """
            key = str(hour) + sector
            if key not in hourly_factors_cache:
                solar = np.round((hour + nc_grid["XLONG"][0,:,:]/180*12) % 24)
                solar = np.array(solar, dtype=int)
                solar[solar==24] = 0
                get_function = np.vectorize(hourly_factors[sector].__getitem__)
                hourly_factors_cache[key] = get_function(solar)
            return hourly_factors_cache[key]

        def get_factor(species, units, sector, time):
            """Return multiplicative factor for given emissions and time."""
            # Unit conversion
            if units == "kg m-2 s-1":
                factor = 1e9 * 3600 / species_info[species_wrf].molmass
                units = "mol km-2 h-1"
             if np.isin(species_info[species_wrf].name_cams,['organic-carbon','black-carbon']):
                factor = 1e9
                units = "ug m-2 s-1"
            else:
                raise ValueError("Unexpected units: %s." % units)
            # Speciation
            try:
                factor *= dict(NO=frac_NO_anthropo,
                               NO2=frac_NO2_anthropo,
                               SO4J=frac_SO4_anthropo,
                               SO2=frac_SO2_anthropo)[species]
            except KeyError:
                pass
            # Temporal projection (time of day and day of week)
            factor *= daily_factors[sector][time.weekday()]
            if args.period == "hour":
                factor = get_hourly_factors(time.hour, sector) * factor
            return factor, units

        def process_emissions(species_wrf, sector, ncs_wrf):
            """Process emissions for given species and sectors."""
            species_cams = species_info[species_wrf].name_cams
            nc_cams = ncs_cams[species_cams]
            if sector not in nc_cams.variables:
                return
            units_cams = nc_cams[sector].units
            # Quality check on the variable's dimensions in input file
            vdims = tuple(d.name for d in nc_cams[sector].get_dims())
            if vdims != ("time", "lat", "lon"):
                raise ValueError("Unexpected dimensions.")
            # Process each time step
            f_datestring = "%Y-%m-%d_%H:%M:%S"
            idx_old = None
            time = start
            while time <= end:
                time_f = time.strftime(f_datestring)
                filepath = "wrfchemi_d%s_%s" % (domain, time_f)
                factor, units_wrf = get_factor(species_wrf, units_cams,
                                               sector, time)
                nc_wrf = get_wrf_emissions_file(domain, time, ncs_wrf, nc_grid)
                var = get_wrf_emissions_variable(nc_wrf, species_wrf, units_wrf)
                idx = unique_index(times_cams, time.strftime("%Y%m"))
                if idx != idx_old:
                    em_in = nc_cams[sector][idx,:,:]
                    # Fix unrealistic NH3 high-lat agricultural emissions
                    if sector == "agl" and species_cams == "ammonia":
                        em_in[lat_in>72] = 0
                    em_out = remap(em_in, mapping)
                var[0,:,:] += em_out * factor
                idx_old = idx
                time += period

        # Here is the loop that actually processes emissions
        ncs_wrf = dict()
        for species_wrf, sector in product(all_species_wrf, all_sectors_cams):
            species_cams = species_info[species_wrf].name_cams
            if species_cams is None:
                # TODO There is no available data, so set emissions to zero
                pass
            else:
                print("Process species %s, sector %s" % (species_wrf, sector))
                process_emissions(species_wrf, sector, ncs_wrf)

        nc_grid.close()
        close_ncfiles_in_dict(ncs_wrf)

    close_ncfiles_in_dict(ncs_cams)
