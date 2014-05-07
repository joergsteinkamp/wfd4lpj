#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import netCDF4 as nc
#import tables

# adjust the paths here, where your files are located
wfd_base_dir  = "/data/external/global/Climate/WFD"
lpj_soil_data = "/data/LPJ/input/soils_lpj.dat"

# spatial resolution in degree
res = 0.5

# global longitudes/latitudes; 
# lon from west to east
# lat from north to south
out_lat = np.arange(90 - res/2, -90 - res/2, -res)
out_lon = np.arange(-180 + res/2, 180 + res/2, res)

# initialize the input. Dont't close the netCDF file,
# data is read, when it is used and not immediately!

# read the WFD land surface file
wfd_ncin    = nc.Dataset(os.path.join(wfd_base_dir, 'WFD-land-lat-long-z.nc'), 'r')
wfd_land_id = wfd_ncin.variables["land"]
wfd_lat     = wfd_ncin.variables["Latitude"]
wfd_lon     = wfd_ncin.variables["Longitude"]
wfd_z       = wfd_ncin.variables["Z"]

wfd_lat_id = (max(out_lat[:]) - wfd_lat[:]) * 2
wfd_lon_id = (wfd_lon[:] - min(out_lon[:])) * 2
wfd_lat_id = np.array(wfd_lat_id.round(0), 'i')
wfd_lon_id = np.array(wfd_lon_id.round(0), 'i')

# convert to 2D
wfd_z_2d = np.zeros((len(out_lat), len(out_lon)) , "f") - 999.9
wfd_z_2d[wfd_lat_id, wfd_lon_id] = wfd_z[:]

# read the LPJ soil data
lpj_soil = pd.io.parsers.read_csv(lpj_soil_data,
                                      skipinitialspace=True, sep=" ", 
                                      header=None, names=["lon", "lat", "type"],
                                      dtype={'lon':np.float32,'lat':np.float32,'type':np.int32})
# lpj coordinates are at the south-western edge of a gridcell
lpj_soil.lon += res/2.0
lpj_soil.lat += res/2.0

lpj_lat_id = (max(out_lat[:]) - lpj_soil.lat[:]) * 2
lpj_lon_id = (lpj_soil.lon[:] - min(out_lon[:])) * 2
lpj_lat_id = np.array(lpj_lat_id.round(0), 'i')
lpj_lon_id = np.array(lpj_lon_id.round(0), 'i')

# convert to 2D
lpj_soil_2d = np.zeros((len(out_lat), len(out_lon)) , "i")
lpj_soil_2d[lpj_lat_id, lpj_lon_id] = lpj_soil.type[:]

# WFD and WFDEI do not fully agree (figured out manually)
# set these 10 locations to 0
lon_exclude = np.array([-106.25, 49.75, 49.25, 46.25, 46.25, 45.75, 48.25, 44.75, -67.25, -67.25])
lat_exclude = np.array([ 77.25,  49.25, 48.75, 47.25, 46.75, 45.25, 40.25,  1.25, -49.25, -49.75])
for i in range(len(lon_exclude)):
  lpj_soil_2d[(np.max(out_lat) - lat_exclude[i]) * 2, (lon_exclude[i] - np.min(out_lon)) * 2] = 0

# take the smallest overlap of WFD and LPJ soil data
lpj_soil_2d[wfd_z_2d == -999.9] = 0
wfd_z_2d[lpj_soil_2d == 0] = -999.9

# create a new index (count along longitude and from north to south)
index_2d = np.arange(1, len(out_lat) * len(out_lon) + 1).reshape(len(out_lat), len(out_lon))
index_2d = np.where(lpj_soil_2d == 0, 0, index_2d)

# create 2D longitude and latitude field
lon_2d = [out_lon] * len(out_lat)
lon_2d = np.array(lon_2d)
lat_2d = [out_lat] * len(out_lon)
lat_2d = np.array(lat_2d)
lat_2d = lat_2d[:].transpose()

#########################
### write the 2D data
#########################
ncout = nc.Dataset(os.path.join(wfd_base_dir, 'LPJ', 'spatial_2D.nc'), 'w')

ncout_dim_lat = ncout.createDimension('lat', len(out_lat))
ncout_dim_lon = ncout.createDimension('lon', len(out_lon))

ncout_lat                    = ncout.createVariable('lat',  'f4', ('lat'))
ncout_lat.standard_name      = "latitude" 
ncout_lat.long_name          = "latitude"
ncout_lat.units              = "degrees_north"
ncout_lon                    = ncout.createVariable('lon', 'f4', ('lon'))
ncout_lon.standard_name      = "longitude" 
ncout_lon.long_name          = "longitude"
ncout_lon.units              = "degrees_east"
ncout_z                      = ncout.createVariable('altitude', 'f4', ('lat', 'lon'), fill_value=-999.9)
ncout_z.long_name            = "altitude"
ncout_z.units                = "m"
ncout_z.missing_value        = -999.9
ncout_soil                   = ncout.createVariable('soiltype', 'i', ('lat', 'lon'), fill_value=0)
ncout_soil.long_name         = "LPJ soilcode"
ncout_soil.units             = "-"
ncout_soil.missing_value     = 0
ncout_index_2d               = ncout.createVariable('index', 'i', ('lat', 'lon'), fill_value=0)
ncout_index_2d.long_name     = "index of gridcell"
ncout_index_2d.units         = "-"
ncout_index_2d.missing_value = 0
ncout_lat[:]                 = out_lat
ncout_lon[:]                 = out_lon
ncout_z[:]                   = wfd_z_2d[:]
ncout_soil[:]                = lpj_soil_2d[:]
ncout_index_2d[:]            = index_2d[:]

ncout_lon_2d = ncout.createVariable('lon_2d', 'f4', ('lat', 'lon'), fill_value=-999.9)
ncout_lat_2d = ncout.createVariable('lat_2d', 'f4', ('lat', 'lon'), fill_value=-999.9)
ncout_lon_2d[:] = lon_2d
ncout_lat_2d[:] = lat_2d
ncout.close()

out_lon_1d = lon_2d[index_2d != 0]
out_lat_1d = lat_2d[index_2d != 0]
out_soil   = lpj_soil_2d[index_2d != 0]
out_index  = index_2d[index_2d != 0]

# write the gridlist file for LPJ input
fout = open(os.path.join(wfd_base_dir, 'LPJ', "gridlist_cf.txt"), 'w')
for i in range(len(out_index)):
  fout.write("%7.2f %7.2f\n" % (out_lon_1d[i], out_lat_1d[i]))
fout.close()

fout = open(os.path.join(wfd_base_dir, 'LPJ', "landid_cf.txt"), 'w')
for i in range(len(out_index)):
  fout.write("%i\n" % (out_index[i]))
fout.close()

#######################################################
### create the NetCDF files for each input valiable
#######################################################
for v in ['Rainf', 'Tair', 'SWdown']:
  if (v == "Rainf"):
    ncout = nc.Dataset(os.path.join(wfd_base_dir, 'LPJ', 'prec.nc'), 'w')
  elif (v == "Tair"):
    ncout = nc.Dataset(os.path.join(wfd_base_dir, 'LPJ', 'temp.nc'), 'w')
  elif (v == "SWdown"):
    ncout = nc.Dataset(os.path.join(wfd_base_dir, 'LPJ', 'insol.nc'), 'w')

  ncout.Conventions = "CF-1.4"

  ncout_dim_id   = ncout.createDimension('land_id', len(out_index))
  ncout_dim_time = ncout.createDimension('time', None)

  ncout_id    = ncout.createVariable('land_id',  'i', ('land_id',))
  ncout_id[:] = out_index

  ncout_time               = ncout.createVariable('time', 'i', ('time',))
  ncout_time.standard_name = "time"
  ncout_time.long_name     = "time"
  ncout_time.units         = "days since 1901-01-01 00:00:00"
  ncout_time.calendar      = "365_day"

  ncout_lon                = ncout.createVariable('lon', 'f', ('land_id',))
  ncout_lon.standard_name  = "longitude" 
  ncout_lon.long_name      = "longitude"
  ncout_lon.units          = "degrees_east"
  ncout_lat                = ncout.createVariable('lat', 'f', ('land_id',))
  ncout_lat.standard_name  = "latitude" 
  ncout_lat.long_name      = "latitude"
  ncout_lat.units          = "degrees_north"
  ncout_soil               = ncout.createVariable('soil', 'f', ('land_id',))
  ncout_soil.long_name     = "LPJ soilcode"
  ncout_soil.units         = "-"

  ncout_lon[:] = out_lon_1d
  ncout_lat[:] = out_lat_1d
  ncout_soil[:] = out_soil

##############################################
### read and write the compressed WFD data
##############################################
  if (v == "Rainf"):
    ncout_clim               = ncout.createVariable('prec', 'f', ('time', 'land_id',))
    ncout_clim.standard_name = "precipitation_amount"
    ncout_clim.long_name     = "Daily precipitation amount"
    ncout_clim.units         = "kg m-2"
  elif (v == "Tair"):
    ncout_clim               = ncout.createVariable('temp', 'f', ('time', 'land_id',))
    ncout_clim.standard_name = "air_temperature"
    ncout_clim.long_name     = "Near surface air temperature at 2m"
    ncout_clim.units         = "K"
  elif (v == "SWdown"):
    ncout_clim               = ncout.createVariable('insol', 'f', ('time', 'land_id',))
    ncout_clim.standard_name = "surface_downwelling_shortwave_flux"
    ncout_clim.long_name     = "Mean daily surface incident shortwave radiation"
    ncout_clim.units         = "W m-2"
  ncout_clim.coordinates = "lat lon"

  days = 0
  for y in range(1901, 1979):
    print(v, y)
    for m in range(1, 13):
      if (v == "Rainf"):
        in_file = "%s_daily_WFD_GPCC_%04i%02i.nc" % (v, y, m)
      else:
        in_file = "%s_daily_WFD_%04i%02i.nc" % (v, y, m)
      in_file = os.path.join(wfd_base_dir, v, in_file)
      clim_ncin    = nc.Dataset(in_file, 'r')
      clim_land_id = clim_ncin.variables["land"]
      clim         = clim_ncin.variables[v]
      if (v=="Rainf"):
        clim         = clim[:] * 86400.0
      else:
        clim         = clim[:]

      if (clim.shape[0] == 29):
        clim = clim[0:28,:]

      if (sum(clim_land_id[:] == wfd_land_id[:]) != len(clim_land_id)):
        print("shape of ID not matching in '%s'" % in_file)
        sys.exit(99)

      clim_2d = np.zeros((len(out_lat), len(out_lon)) , "f") - 999.9
      
      for d in np.arange(clim.shape[0]):
        clim_2d[wfd_lat_id, wfd_lon_id] = clim[d,:]
        ncout_clim[days,:] = np.array(clim_2d[index_2d != 0], 'f')
        days += 1
    ncout.sync()

##############################################
### read and write the WFDEI data
##############################################

  for y in range(1979, 2001):
    print(v, y)
    for m in range(1, 13):
      if (v == "Rainf"):
        in_file = "%s_daily_WFDEI_GPCC_%04i%02i.nc" % (v, y, m)
      else:
        in_file = "%s_daily_WFDEI_%04i%02i.nc" % (v, y, m)
      in_file = os.path.join(wfd_base_dir, "WFDEI", v, in_file)
      clim_ncin    = nc.Dataset(in_file, 'r')
      clim         = clim_ncin.variables[v]
      if (v=="Rainf"):
        clim         = clim[:] * 86400.0
      else:
        clim         = clim[:]

      if (clim.shape[0] == 29):
        clim = clim[0:28,:,:]

      for d in np.arange(clim.shape[0]):
        clim_2d = clim[d,:,:]
        clim_2d[clim_2d>1.e19] = -999.9
        clim_2d = np.flipud(clim_2d)
        clim_out = np.array(clim_2d[index_2d != 0], 'f')
        ncout_clim[days,:] = clim_out
        days += 1
    ncout.sync()

##############################################
### write the time var and close
##############################################

  ncout_time[:] = np.arange(days) + 1
  ncout.close()
