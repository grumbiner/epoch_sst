from math import *
import os
import datetime
import copy

import numpy as np
import numpy.ma as ma
import netCDF4

#------------------------------------------------
def applymask(mask, grid, indices):
  for k in range(0, len(indices[0])):
    i = indices[1][k]
    j = indices[0][k]
    grid[j,i] = 0.   
#-------------------------------------------------
import ncoutput

def writeout(sumx1, sumx2, sumx3, sumx4, tmax, tmin, base, tag, n = 30):
  #RG: write out mean, max, min to save file
  mask =  ma.masked_array(sumx1 < -900.*n)
  indices = mask.nonzero()
  
  applymask(mask, sumx1, indices)
  applymask(mask, sumx2, indices)
  applymask(mask, sumx3, indices)
  applymask(mask, sumx4, indices)
  
  print("sumx1", sumx1.max(), sumx1.min() )
  print("sumx2", sumx2.max(), sumx2.min() )
  print("sumx3", sumx3.max(), sumx3.min() )
  print("sumx4", sumx4.max(), sumx4.min() )
  print("tmax", tmax.max() , tmax.min() )
  print("tmin", tmin.max() , tmin.min() )

  name = base+"traditional_"+tag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('sumx1', dtype = sumx1.dtype)
  foroutput.addvar('mean', dtype = sumx1.dtype)
  foroutput.addvar('sumx2', dtype = sumx2.dtype)
  foroutput.addvar('sumx3', dtype = sumx3.dtype)
  foroutput.addvar('sumx4', dtype = sumx4.dtype)
  foroutput.addvar('tmax', dtype = tmax.dtype)
  foroutput.addvar('tmin', dtype = tmin.dtype)
  
  mean = sumx1/n
  applymask(mask, mean, indices)
  
  foroutput.encodevar(sumx1, 'sumx1')
  foroutput.encodevar(mean,  'mean')
  foroutput.encodevar(sumx2, 'sumx2')
  foroutput.encodevar(sumx3, 'sumx3')
  foroutput.encodevar(sumx4, 'sumx4')
  foroutput.encodevar(tmin,   'tmin')
  foroutput.encodevar(tmax,   'tmax')
  
  tmask = np.zeros((ny,nx))
  for k in range(0, len(indices[0]) ):
      i = indices[1][k]
      j = indices[0][k]
      tmask[j,i] = 1.0
  foroutput.addvar('mask', dtype = tmask.dtype)
  foroutput.encodevar(tmask, 'mask')
  
  foroutput.close()
#------------------ End writeout ---------- --------------------------


#-------------------------------------------------
#  Compute a traditional style climatology, day by day for 30 years
#-------------------------------------------------
# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"

# Defining the quarter degree grid
nx = 1440
ny = 720

# Start-finish, but will be iterating through next 30 years
start = datetime.datetime(1981,9,1)
end = datetime.datetime(1982,8,31)

dt = datetime.timedelta(1)

#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = start
count = 0

while (tag <= end ):
  print("tag =",tag, flush=True)
  # Initialize files for accumulations
  sst = np.zeros((ny,nx)) # temporary file for reading in data
  tmp = np.zeros((ny,nx))
  
  # for accumulating moments:
  sumx1 = np.zeros((ny,nx))
  sumx2 = np.zeros((ny,nx))
  sumx3 = np.zeros((ny,nx))
  sumx4 = np.zeros((ny,nx))

  # extrema
  tmax = np.zeros((ny,nx))
  tmin = np.zeros((ny,nx))
  tmax.fill(-3.0)
  tmin.fill(45.0)

  # Iterate over the 30 years for this day
  for yy in range(0, 30):
    tagyy = datetime.datetime(tag.year+yy, tag.month, tag.day)

# Get the day's data:
    fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
    tmpnc = netCDF4.Dataset(fbase + fname)
    sst = tmpnc.variables['sst'][0,0,:,:]
    if ( count ==  0 ):
        lons = tmpnc.variables['lon'][:]
        lats = tmpnc.variables['lat'][:]
    tmpnc.close()

# Accumulate moments:
    tmp = copy.deepcopy(sst)
    sumx1 += tmp
    tmp *= sst
    sumx2 += tmp
    tmp *= sst
    sumx3 += tmp
    tmp *= sst
    sumx4 += tmp

# Find extrema:
    tmax = np.fmax(tmax, sst)
    tmin = np.fmin(tmin, sst)
    
  writeout(sumx1, sumx2, sumx3, sumx4, tmax, tmin, fbase, tag)

  count += 1   # number of days' data
  tag   += dt
  
#-------------------------------------------------
