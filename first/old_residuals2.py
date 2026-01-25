'''
Compute the residuals for a decade following the 30 years used for the traditional climatology
'''

#from math import *
import datetime
import copy

import numpy as np
from numpy import ma
import netCDF4

from functions import applymask
import ncoutput

#------------------------------------------------
def writeout(flons, flats, fanomaly, base, ftag, n = 28):
  '''
  #RG: write out traditional mean, max, min to save file
  '''
  mask =  ma.masked_array(fanomaly < -900.*n)
  indices = mask.nonzero()

  applymask(fanomaly, indices)

  print("fanomaly", fanomaly.max(), fanomaly.min() )

  name = base+"res_traditional_"+ftag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, flats, flons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('fanomaly', dtype = fanomaly.dtype)

  foroutput.encodevar(fanomaly, 'fanomaly')

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
# location of data files
fbase = "/Volumes/Data2/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"

# Defining the quarter degree grid
nx = 1440
ny = 720

dt = datetime.timedelta(1)

# Start-finish for traditional climatology, write out anomalies for 30 years
#   and writing out residuals
epoch = datetime.datetime(1981,9,1)
end = epoch + 364*dt

#---------------------------------------------

tag = epoch
count = 0
while (tag <= end ):
  print("tag =",tag, flush=True)

# Get the day's climatological data:
  fname = "traditional_" + tag.strftime("%Y%m%d") + ".nc"
  tmpnc = netCDF4.Dataset(fbase + fname)
  mean = tmpnc.variables['mean'][:,:]
  tmpnc.close()

  for nn in range(0,30):
    tagyy = datetime.datetime(tag.year+nn, tag.month, tag.day)
    fname = "oisst-avhrr-v02r01." + tagyy.strftime("%Y%m%d")+".nc"
    tmpnc = netCDF4.Dataset(fbase+fname)
    sst   = tmpnc.variables['sst'][0,0,:,:]
    if ( count ==  0 ):
        lons = tmpnc.variables['lon'][:]
        lats = tmpnc.variables['lat'][:]
    tmpnc.close()
    count += 1

    sst -= mean
    if (sst.max() == sst.min() == 0.0):
      print("day is exactly equal to climatology, exiting",tagyy)

    writeout(lons, lats, sst, fbase, tagyy, n = 10)

  tag   += dt

#-------------------------------------------------
