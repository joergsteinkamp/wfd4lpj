#!/usr/bin/env python

import os
import sys
import numpy as np
import netCDF4 as nc
#import tables

wfd_base_dir  = "/data/external/global/Climate/WFD"

res = 0.5

# global longitudes/latitudes
out_lat = np.arange(90 - res/2, -90 - res/2, -res)
out_lon = np.arange(-180 + res/2, 180 + res/2, res)

# initialize the input. Dont't close the netCDF file,
# data is read, when it is used and not immediately!

# read the WFD land surface file
wfd_ncin    = nc.Dataset(os.path.join(wfd_base_dir, 'WFD-land-lat-long-z.nc'), 'r')
wfd_land_id = wfd_ncin.variables["land"]
wfd_lat     = wfd_ncin.variables["Latitude"]
wfd_lon     = wfd_ncin.variables["Longitude"]

lat_id = (max(out_lat[:]) - wfd_lat[:]) * 2
lon_id = (wfd_lon[:] - min(out_lon[:])) * 2
lat_id = np.array(lat_id.round(0), 'i')
lon_id = np.array(lon_id.round(0), 'i')

# read and write the climate data
for v in ['Rainf', 'Tair', 'SWdown']:
   # process the WFD data from 1901-1978
  for y in range(1979, 2002):
    for m in range(1, 13):

      if (v == "Rainf"):
        in_file = "%s_daily_WFD_GPCC_%04i%02i.nc" % (v, y, m)
      else:
        in_file = "%s_daily_WFD_%04i%02i.nc" % (v, y, m)
      in_file      = os.path.join(wfd_base_dir, v, in_file)
      #print(in_file)
      clim_ncin    = nc.Dataset(in_file, 'r')
      clim_land_id = clim_ncin.variables["land"]
      clim         = clim_ncin.variables[v]
      if (v=="Rainf"):
        clim = clim[:] * 86400.0

      if (sum(clim_land_id[:] == wfd_land_id[:]) != len(clim_land_id)):
        print("shape of ID not matching in '%s'" % in_file)
        sys.exit(99)
      out_clim = np.zeros((clim.shape[0], len(out_lat), len(out_lon)) , "f") - 999.9
      out_clim[:, lat_id, lon_id] = clim[:]

      # open netcdf output file
      if (v == "Rainf"):
        out_file = "%s_daily_WFD_GPCC_%04i%02i.nc" % (v, y, m)
      else:
        out_file = "%s_daily_WFD_%04i%02i.nc" % (v, y, m)
      out_file = os.path.join(wfd_base_dir, 'decomp', out_file)
      ncout = nc.Dataset(out_file, 'w')

      ncout_dim_lat  = ncout.createDimension('lat', len(out_lat))
      ncout_dim_lon  = ncout.createDimension('lon', len(out_lon))
      ncout_dim_time = ncout.createDimension('time', 0)

      ncout_lat            = ncout.createVariable('lat',  'f4', ('lat'))
      ncout_lat.long_name  = "latitude"
      ncout_lat.units      = "degrees_north"
      ncout_lon            = ncout.createVariable('lon', 'f4', ('lon'))
      ncout_lon.long_name  = "longitude"
      ncout_lon.units      = "degrees_east"
      ncout_time           = ncout.createVariable('time', 'i', ('time'))
      ncout_time.long_name = "time"
      ncout_time.units     = "days since %04i-%02i-01 00:00:00" % (y, m)
      ncout_time.calendar  = "standard"

      ncout_clim                      = ncout.createVariable(v, 'f4', ('time', 'lat', 'lon'), fill_value=-999.9)
      if (v == "Rainf"):
        ncout_clim.long_name          = "precipitation_amount"
        ncout_clim.units                = "kg m-2"
      elif (v == "Tair"):
        ncout_clim.long_name          = "air_temperature"
        ncout_clim.units                = "K"
      elif (v == "SWdown"):
        ncout_clim.long_name          = "surface_downwelling_shortwave_flux"
        ncout_clim.units                = "W m-2"
      ncout_clim.missing_value        = -999.9

      ncout_lat[:]         = out_lat
      ncout_lon[:]         = out_lon
      ncout_time[:]        = range(1, clim.shape[0] + 1)
      if (v=="Rainf"):
        ncout_clim[:]        = out_clim
      else:
        ncout_clim[:]        = out_clim
      ncout.sync()
      ncout.close()
