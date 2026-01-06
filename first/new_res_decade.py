'''
new
'''

import datetime
import copy
from math import pi

import numpy as np
from numpy import ma
import netCDF4

from functions import *
import ncoutput

#------------------------------------------------
def writeout(flons, flats, fsumx1, fsumx2, fsumx3, fsumx4, f2base, ftag, n = 28):
  '''
  #RG: write out mean, max, min to save file
  '''
  fmask =  ma.masked_array(fsumx1 < -900.*n)
  indices = fmask.nonzero()

  applymask(fsumx1, indices)
  applymask(fsumx2, indices)
  applymask(fsumx3, indices)
  applymask(fsumx4, indices)

  print("sumx1", fsumx1.max(), fsumx1.min() )
  print("sumx2", fsumx2.max(), fsumx2.min() )
  print("sumx3", fsumx3.max(), fsumx3.min() )
  print("sumx4", fsumx4.max(), fsumx4.min() )

  name = f2base+"res_new_decade_"+ftag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, flats, flons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('sumx1', dtype = fsumx1.dtype)
  foroutput.addvar('mean',  dtype = fsumx1.dtype)
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
#  Work from new climatology, first pass
#-------------------------------------------------
# location of data files
fbase = "/Volumes/Data2/qdoi/v2.1.nc/"

# Defining the quarter degree grid
nx = 1440
ny = 720
loy       = 365.2422 #days, tropical year
freq_base = 2.*pi/loy

dt = datetime.timedelta(1)

# Read in new style climatology -------------------------
epoch = datetime.datetime(1981,9,1)

dset = netCDF4.Dataset(f"epoch{epoch.year:4d}.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]

mask      = dset.variables['mask'][:,:]
mean      = dset.variables['mean'][:,:]
slope     = dset.variables['slope'][:,:]
intercept = dset.variables['intercept'][:,:]

ampl = np.zeros((3,ny,nx))
phas = np.zeros((3, ny, nx))
freq = np.zeros((3))

ampl[0] = dset.variables['cpy1_amp'][:,:]
ampl[1] = dset.variables['cpy2_amp'][:,:]
ampl[2] = dset.variables['cpy3_amp'][:,:]
phas[0] = dset.variables['cpy1_pha'][:,:]
phas[1] = dset.variables['cpy2_pha'][:,:]
phas[2] = dset.variables['cpy3_pha'][:,:]
freq[0] = freq_base
freq[1] = freq_base*2
freq[2] = freq_base*3
phas *= pi/180.


# Ensuing decade
resstart = datetime.datetime(2011,1,1)
ressend  = datetime.datetime(2011,12,31)

#---------------------------------------------
# Now run through the data files and accumulate terms:

# for accumulating moments through whole year:
sumx1 = np.zeros((ny,nx))
sumx2 = np.zeros((ny,nx))
sumx3 = np.zeros((ny,nx))
sumx4 = np.zeros((ny,nx))

tag = resstart
count = 0

while (tag <= ressend ):
  #debug: print("tag =",tag, flush=True)
  # Initialize files for accumulations
  sst = np.zeros((ny,nx)) # temporary file for reading in data
  mean = np.zeros((ny,nx))
  tmp = np.zeros((ny,nx))

# Get the day's climatological data:
  tclim = climo(intercept, slope, ampl, phas, freq, epoch, tag)

  fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d")+".nc"
  tmpnc = netCDF4.Dataset(fbase+fname)
  sst   = tmpnc.variables['sst'][0,0,:,:]
  tmpnc.close()

  print(tag.strftime("%Y%m%d"),sst.max(), sst.min(), tclim.max(), tclim.min(), ' ', end="")
  sst -= tclim
  print(sst.max(), sst.min() )

# Accumulate moments:
  tmp = copy.deepcopy(sst)
  sumx1 += tmp
  tmp *= sst
  sumx2 += tmp
  tmp *= sst
  sumx3 += tmp
  tmp *= sst
  sumx4 += tmp

  count += 1
  tag   += dt

tagyy = datetime.datetime(2025,1,1)
writeout(lons, lats, sumx1, sumx2, sumx3, sumx4, fbase, tagyy, n = 10*365)
#-------------------------------------------------
