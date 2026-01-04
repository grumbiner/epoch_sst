'''
First pass of traditional climatology
'''

#from math import *
import datetime
import copy

import numpy as np
from numpy import ma
import netCDF4

from functions import *
import ncoutput

#------------------------------------------------
def writeout(fsumx1, flats, flons, f2base, ftag, fn = 30):
  '''
  #RG: write out mean to save file
  '''
  mask =  ma.masked_array(fsumx1 < -900.*fn)
  indices = mask.nonzero()

  name = f2base+"traditional_"+ftag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, flats, flons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('mean', dtype = fsumx1.dtype)

  mean = fsumx1/fn
  applymask(mean, indices)

  foroutput.encodevar(mean,  'mean')

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
fbase = "/Volumes/Data2/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"

# Defining the quarter degree grid
nx = 1440
ny = 720

dt = datetime.timedelta(1)
# Start-finish, but will be iterating through next 30 years
epoch = datetime.datetime(1981,9,1)
end = epoch + 364*dt

#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = epoch
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

  writeout(sumx1, lats, lons, fbase, tag, fn = 30)

  count += 1   # number of days' data
  tag   += dt

#-------------------------------------------------
