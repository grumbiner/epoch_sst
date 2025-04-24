import copy
import datetime
from math import *

import numpy as np
import numpy.ma as ma

import netCDF4 as nc

"""
#Offline:
#  Compute + map:
#    %variance explained by mean, trend, harmonics, Nino3.4
#    Magnitude residual variance

"""
 
#----------------------------------------------------------------------

from functions import *
import ncoutput

def writeout(tsst, mask, nx, ny, lats, lons, tag):
  #indices = mask.nonzero()
  #applymask(mask, tsst, indices)
  print("tsst ",tag, tsst.max(), tsst.min(), tsst.mean() )

  name = "v2.1.nc/ninores1_"+tag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('ninores1', dtype = tsst.dtype)
  foroutput.encodevar(tsst, 'ninores1')

  foroutput.close()

#----------------------------------------------------------------------
nx = 1440
ny =  720
loy = 365.2422 # tropical year
freq_base = 2.*pi/loy

dset = nc.Dataset("first_pass.nc", "r")
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


dset = nc.Dataset("newres1_30.nc","r")
sumx1 = dset.variables['sumx1'][:,:]
sumx2 = dset.variables['sumx2'][:,:]
sumx3 = dset.variables['sumx3'][:,:]
sumx4 = dset.variables['sumx4'][:,:]
sumxn = dset.variables['sumxn'][:,:]
sumn  = dset.variables['sumn'][:,:]
sumn2 = dset.variables['sumn2'][:,:]

days = 10957

tmpn = days*sumn2 - sumn*sumn
print(tmpn.min(), tmpn.max(), tmpn.mean() )

nino_slope     = (days*sumxn - sumx1*sumn) / tmpn
print('nino_slope',nino_slope.max(), nino_slope.min(), nino_slope.mean(), flush=True )

epoch = datetime.datetime(1981,9,1)

#-------------------------------------------------
fbase = "/Volumes/Data/qdoi/v2.1.nc/"

# Initialize files for accumulations
sst = np.zeros((ny,nx)) # temporary file for reading in data
sumx1 = np.zeros((ny,nx))

# Original span:
#start = datetime.datetime(1981,9,1)
#debug: end   = datetime.datetime(1981,9,18)
#end   = datetime.datetime(2011,8,31)
# Next Decade
start = datetime.datetime(2011,9,1)
end   = datetime.datetime(2021,8,31)

dt = datetime.timedelta(1)
tag = start
count = 0
while (tag <= end):
  if (count % 90 == 0):
    print(tag, flush=True)

# Get the day's data:
  fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
  tmpnc = nc.Dataset(fbase + fname)
  sst = tmpnc.variables['sst'][0,0,:,:]
  if ( count ==  0 ):
      lons = tmpnc.variables['lon'][:]
      lats = tmpnc.variables['lat'][:]
      zlon = np.logical_and(lons > 190, lons < 240)
      zlat = np.logical_and(lats > -5., lats < 5.)
      nino34 = np.zeros((ny,nx),dtype=bool)
      nino34[zlat,:] = True
      nino34[:,zlon] &= True
      del zlon, zlat
  tmpnc.close()

# Accumulate moments:
  tsst = copy.deepcopy(sst)
  tclim = climo(intercept, slope, ampl, phas, freq, epoch, tag)
  tsst -= tclim

# Nino3.4 orthogonalization info
  tnino34 = tsst[nino34].mean() 
  tnino34 *= nino_slope
  tsst -= tnino34

  sumx1 += tsst

  count += 1   # number of days' data
  tag   += dt 
#------------------------------------------------
days = count

name = "ninores_30.nc"
sumx1 /= days
foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('mean', dtype = sumx1.dtype)
foroutput.encodevar(sumx1,  'mean')
foroutput.encodescalar(count, 'days')

foroutput.close()


# ---- .nc encoding --------------------------------------------------
import ncoutput



#------------------ End of second pass --------------------------
