'''
Compute the residuals from the traditional climatology
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
def writeout(flons, flats, fsumx1, fsumx2, fsumx3, fsumx4, base, ftag, n = 28):
  '''
  #RG: write out traditional mean, max, min to save file
  '''
  mask =  ma.masked_array(fsumx1 < -900.*n)
  indices = mask.nonzero()

  applymask(fsumx1, indices)
  applymask(fsumx2, indices)
  applymask(fsumx3, indices)
  applymask(fsumx4, indices)

  print("sumx1", fsumx1.max(), fsumx1.min() )
  print("sumx2", fsumx2.max(), fsumx2.min() )
  print("sumx3", fsumx3.max(), fsumx3.min() )
  print("sumx4", fsumx4.max(), fsumx4.min() )

  name = base+"res_traditional_"+ftag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, flats, flons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('sumx1', dtype = fsumx1.dtype)
  foroutput.addvar('mean', dtype = fsumx1.dtype)
  foroutput.addvar('sumx2', dtype = fsumx2.dtype)
  foroutput.addvar('sumx3', dtype = fsumx3.dtype)
  foroutput.addvar('sumx4', dtype = fsumx4.dtype)

  fmean = fsumx1/n
  applymask(fmean, indices)

  foroutput.encodevar(fsumx1, 'sumx1')
  foroutput.encodevar(fmean,  'mean')
  foroutput.encodevar(fsumx2, 'sumx2')
  foroutput.encodevar(fsumx3, 'sumx3')
  foroutput.encodevar(fsumx4, 'sumx4')

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

# Start-finish for traditional climatology, but will be iterating through following decade
start = datetime.datetime(1981,9,1)
end = datetime.datetime(1982,8,31)

dt = datetime.timedelta(1)

#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = start
count = 0

# for accumulating moments through whole year:
yrsumx1 = np.zeros((ny,nx))
yrsumx2 = np.zeros((ny,nx))
yrsumx3 = np.zeros((ny,nx))
yrsumx4 = np.zeros((ny,nx))

while (tag <= end ):
  print("tag =",tag, flush=True)

# Get the day's climatological data:
  fname = "traditional/traditional_" + tag.strftime("%Y%m%d") + ".nc"
  tmpnc = netCDF4.Dataset(fbase + fname)
  mean = tmpnc.variables['mean'][:,:]
  tmpnc.close()

  # for accumulating moments:
  sumx1 = np.zeros((ny,nx))
  sumx2 = np.zeros((ny,nx))
  sumx3 = np.zeros((ny,nx))
  sumx4 = np.zeros((ny,nx))

  # Iterate over the 10 years for this day
  for yy in range(0, 10):

    tagyy = datetime.datetime(tag.year+30+yy, tag.month, tag.day)
    fname = "oisst-avhrr-v02r01." + tagyy.strftime("%Y%m%d")+".nc"
    tmpnc = netCDF4.Dataset(fbase+fname)
    sst   = tmpnc.variables['sst'][0,0,:,:]
    if ( count ==  0 ):
        lons = tmpnc.variables['lon'][:]
        lats = tmpnc.variables['lat'][:]
    tmpnc.close()

    sst -= mean

# Accumulate moments:
    tmp = copy.deepcopy(sst)
    sumx1 += tmp
    tmp *= sst
    sumx2 += tmp
    tmp *= sst
    sumx3 += tmp
    tmp *= sst
    sumx4 += tmp


  writeout(lons, lats, sumx1, sumx2, sumx3, sumx4, fbase, tagyy, n = 10)
  yrsumx1 += sumx1
  yrsumx2 += sumx2
  yrsumx3 += sumx3
  yrsumx4 += sumx4

  count += 1
  tag   += dt

tagyy = datetime.datetime(2025,1,1)
writeout(lons, lats, yrsumx1, yrsumx2, yrsumx3, yrsumx4, fbase, tagyy, n = 10*365)
#-------------------------------------------------
