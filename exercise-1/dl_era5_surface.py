import os
import cdsapi
from calendar import monthrange


#-------- Input --------
STARTYEAR=2019
ENDYEAR=2019
STARTMONTH=1
ENDMONTH=1
MAXLAT=75.0
MINLON=-10.0
MINLAT=10.0
MAXLON=50.0


#-------- Download ERA5 data from CDS --------
for iyear in range(STARTYEAR,ENDYEAR+1):
  era5year = str(iyear)
  if(not os.path.isdir('ERA5_grib1_'+era5year)):
    os.mkdir('ERA5_grib1_'+era5year)
  os.chdir('ERA5_grib1_'+era5year)
  for imonth in range(STARTMONTH,ENDMONTH+1):
    era5month = str(imonth).zfill(2)
    ndays_in_month = monthrange(int(era5year),imonth)[1]
    for iday in range(1,ndays_in_month+1):
      era5day = str(iday).zfill(2)
      era5filename = 'e5.sfc.'+era5year+era5month+era5day+'.grib'
      if (not os.path.isfile(era5filename)):
        c = cdsapi.Client()
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                    '2m_temperature', 'mean_sea_level_pressure', 'sea_ice_cover',
                    'sea_surface_temperature', 'skin_temperature', 'snow_depth',
                    'soil_temperature_level_1', 'soil_temperature_level_2', 'soil_temperature_level_3',
                    'soil_temperature_level_4', 'surface_pressure', 'volumetric_soil_water_layer_1',
                    'volumetric_soil_water_layer_2', 'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4',
                ],
                'year': era5year,
                'month': era5month,
                'day': era5day,
                'time': [
                    '00:00', '06:00', '12:00',
                    '18:00',
                ],
                'area': [
                    MAXLAT, MINLON, MINLAT,
                    MAXLON,
                ],
                'format': 'grib',
            },
            era5filename)
  os.chdir('..')
