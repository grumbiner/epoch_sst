import copy
import datetime
from math import *

import numpy as np
import numpy.ma as ma

import netCDF4 as nc

"""
Read in first pass -- 
  intercept, trend, harmonics 1-3 and their phase
fn to compute Tclim(tau) given the above (tau = days since 1 Sep 1981)

Second pass -- 
  Read in daily analyses
    subtract climatology
    accumulate stats on residuals
    accumulate terms for orthogonalizing w.r.t. Nino3.4
    write out residual field for the day
  Write out statistics on residuals
  Write out statistics on Nino3.4 orthogonalizing
  Maps of deviation norms
  Maps of correlations to Nino3.4

Third pass --
  Read in daily analyses
    subtract climatology
    accumulate stats on residuals
    write out residual field for the day
  Write out stats on residuals
  Maps of deviation norms

Offline:
  Compute + map:
    %variance explained by mean, trend, harmonics, Nino3.4
    Magnitude residual variance

"""
 
#----------------------------------------------------------------------

from functions import *

#----------------------------------------------------------------------
nx = 1440
ny =  720
loy = 365.2422 # tropical year
freq_base = 2.*pi/loy

dset = nc.Dataset("first_pass.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]

fmask     = dset.variables['mask'][:,:]
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

print(ampl[0].max(), phas[0].max() )

epoch = datetime.datetime(1981,9,1)
#tag   = datetime.datetime(2021,9,1)
tag   = datetime.datetime(1981,9,1)

#debug: sst = climo(intercept, slope, ampl, phas, freq, epoch, tag)
#debug: print(sst.max(), sst.min(), sst.mean() )
#debug: print(sst[sst < -1.8])

'''

Second pass --
  Read in first pass and prep for climatology computation (above)
  Read in daily analyses
    subtract climatology
    accumulate stats on residuals
    accumulate terms for orthogonalizing w.r.t. Nino3.4
    write out residual field for the day
  Write out statistics on residuals
  Write out statistics for Nino3.4 orthogonalizing
  Maps of deviation norms
  Maps of correlations to Nino3.4

'''

#-------------------------------------------------
fbase = "/Volumes/Data/qdoi/v2.1.nc/"

# Initialize files for accumulations
sst = np.zeros((ny,nx)) # temporary file for reading in data
tmp = np.zeros((ny,nx))

# for accumulating moments:
sumx1 = np.zeros((ny,nx))
sumx2 = np.zeros((ny,nx))
sumx3 = np.zeros((ny,nx))
sumx4 = np.zeros((ny,nx))
#debug: print('dtype for sumx1 ',sumx1.dtype, flush=True)

# for nino3.4 orthogonalization
sumxn = np.zeros((ny,nx))
sumn  = 0.0
sumn2 = 0.0

start = datetime.datetime(2011,9,1)
#end   = datetime.datetime(2024,12,25)
#end   = datetime.datetime(2011,12,25)
end   = datetime.datetime(2021,8,31)

dt = datetime.timedelta(1)
tag = start
count = 0
while (tag <= end):
  if (count % 30 == 0):
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
  
  sumx1 += tsst
  sumx2 += (tsst*tsst)
  sumx3 += (tsst*tsst*tsst)
  sumx4 += (tsst*tsst)*(tsst*tsst)

# Accumulate Nino3.4 orthogonalization info
  tnino34 = tsst[nino34].mean() 
  sumxn += tnino34*tsst
  sumn2 += tnino34*tnino34
  sumn  += tnino34

  del tclim
  count += 1   # number of days' data
  tag   += dt 
#------------------------------------------------
indices = fmask.nonzero()
applymask(fmask, sumx1, indices)
applymask(fmask, sumx2, indices)
applymask(fmask, sumx3, indices)
applymask(fmask, sumx4, indices)
# orthog1
# orthog2

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )
print("sumxn", sumxn.max(), sumxn.min() )
mean = sumx1 / count
print("mean", mean.max(), mean.min() )


# ---- .nc encoding --------------------------------------------------
import ncoutput

name = "second_pass.nc"

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('mean', dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)
foroutput.addvar('sumxn', dtype = sumxn.dtype)

foroutput.addvar('mask', dtype = fmask.dtype)
foroutput.encodevar(fmask, 'mask')

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(mean,  'mean')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')
foroutput.encodevar(sumxn, 'sumxn')

print("number of days = ",count)
foroutput.encodescalar(count, 'days')
foroutput.encodescalar(sumn2, 'sumn2')
foroutput.encodescalar(sumn, 'sumn')

foroutput.close()
#------------------ End of second pass --------------------------
