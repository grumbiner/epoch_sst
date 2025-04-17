import copy
import datetime
from math import *

import numpy as np
import numpy.ma as ma

import netCDF4 as nc

"""
#Read in first pass -- 
#  intercept, trend, harmonics 1-3 and their phase

#fn to compute Tclim(tau) given the above (tau = days since 1 Sep 1981)

!Second pass -- 
!  Read in daily analyses
!    subtract climatology
!    write out residual field for the day
!    accumulate stats on residuals
!    accumulate terms for orthogonalizing w.r.t. Nino3.4
!  Write out statistics on residuals
!  Write out statistics on Nino3.4 orthogonalizing
!  Maps of deviation norms
!  Maps of correlations to Nino3.4

#Third pass --
#  Read in daily analyses
#    subtract climatology
#    accumulate stats on residuals
#    write out residual field for the day
#  Write out stats on residuals
#  Maps of deviation norms

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

  name = "v2.1.nc/newres1_"+tag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
  foroutput.ncoutput(name)
  foroutput.addvar('newres1', dtype = tsst.dtype)
  foroutput.encodevar(tsst, 'newres1')

  foroutput.close()

#----------------------------------------------------------------------
nx = 1440
ny =  720
loy = 365.2422 # tropical year
freq_base = 2.*pi/loy

dset = nc.Dataset("first_pass.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]

mask     = dset.variables['mask'][:,:]
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
sumn  = np.zeros((ny,nx))
sumn2 = np.zeros((ny,nx))

# Original span:
start = datetime.datetime(1981,9,1)
#debug: end   = datetime.datetime(1981,9,18)
end   = datetime.datetime(2011,8,31)
# Next Decade
#start = datetime.datetime(2011,9,1)
#end   = datetime.datetime(2021,8,31)

dt = datetime.timedelta(1)
tag = start
count = 0
while (tag <= end):
  if (count % 30 == 0):
    print(tag, flush=True)

# Get the day's data:
  fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
  #fname = "new_residual1/newres1_" + tag.strftime("%Y%m%d") + ".nc"
  tmpnc = nc.Dataset(fbase + fname)
  sst = tmpnc.variables['sst'][0,0,:,:]
  #sst = tmpnc.variables['newres1'][:,:]
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
# write out residual, tsst
  writeout(tsst, mask, nx, ny, lats, lons, tag)
  
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
days = count

indices = mask.nonzero()
applymask(mask, sumx1, indices)
applymask(mask, sumx2, indices)
applymask(mask, sumx3, indices)
applymask(mask, sumx4, indices)
applymask(mask, sumxn, indices)
applymask(mask, sumx2, indices)
applymask(mask, sumn,  indices)
# orthog1
# orthog2

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )
print("sumxn", sumxn.max(), sumxn.min() )
print("sumn2", sumn2.max(), sumn2.min() )
print("sumn", sumn.max(), sumn.min() )
mean = sumx1 / count
print("mean", mean.max(), mean.min() )


# ---- .nc encoding --------------------------------------------------
import ncoutput

name = "newres1_30.nc"

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('mean', dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)
foroutput.addvar('sumxn', dtype = sumxn.dtype)
foroutput.addvar('sumn2', dtype = sumxn.dtype)
foroutput.addvar('sumn', dtype = sumxn.dtype)

foroutput.addvar('mask', dtype = mask.dtype)
foroutput.encodevar(mask, 'mask')

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(mean,  'mean')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')
foroutput.encodevar(sumxn, 'sumxn')
foroutput.encodevar(sumn2, 'sumn2')
foroutput.encodevar(sumn, 'sumn')

tmpn = days*sumn2 - sumn*sumn
tmpx = days*sumx2 - sumx1*sumx1
tmpx[tmpx == 0 ] = 1
print(tmpn.min(), tmpn.max(), tmpn.mean() )
print(tmpx.min(), tmpx.max(), tmpx.mean() )

slope     = (days*sumxn - sumx1*sumn) / tmpn
print('slope',slope.max(), slope.min(), slope.mean(), flush=True )

intercept = (sumx1/days - slope*sumn/days)
print('intercept',intercept.max(), intercept.min(), intercept.mean(), flush=True )

correl    = (days*sumxn - sumx1*sumn) / np.sqrt(tmpn) / np.sqrt(tmpx)
print('correl',correl.max(), correl.min(), correl.mean(), flush=True )



print("number of days = ",count)
foroutput.encodescalar(count, 'days')

foroutput.close()
#------------------ End of second pass --------------------------
