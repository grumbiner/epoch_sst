"""
Doesn't really fit its name. Computes a second pass estimating moments 
  (sumx2 etc.) using deviations from climatology, hopefully more accurately 
  than first_pass working with full temperatures.
"""

import copy
import datetime

import numpy as np
import netCDF4 as nc

#----------------------------------------------------------------------
from functions import *
import ncoutput

#----------------------------------------------------------------------

def writeout(ftsst, fnx, fny, flats, flons, ftag):
  ''' writing out the sst -- writeout(sst, nx, ny, lats, lons, tag) '''
  print("tsst ",tag.strftime("%Y%m%d"), ftsst.max(), ftsst.min(), ftsst.mean() )

  f2name = "v2.1.nc/oldres1_"+ftag.strftime("%Y%m%d")+".nc"

  foroutput = ncoutput.ncoutput(fnx, fny, flats, flons, f2name)
  foroutput.ncoutput(f2name)
  foroutput.addvar('oldres1', dtype = ftsst.dtype)
  foroutput.encodevar(ftsst, 'oldres1')

  foroutput.close()


#----------------------------------------------------------------------
nx = 1440
ny =  720
dt = datetime.timedelta(1)

epoch = datetime.datetime(1981,9,1)
#end   = datetime.datetime(2020,12,31)
end   = datetime.datetime(2011,8,31)

dset = nc.Dataset(f"epoch{epoch.year:4d}.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
fmask  = dset.variables['mask'][:,:]
dset.close()

#-------------------------------------------------
fbase = "/Volumes/Data2/qdoi/v2.1.nc/"

# Initialize files for accumulations
sst = np.zeros((ny,nx)) # temporary file for reading in data
tmp = np.zeros((ny,nx))

# for accumulating moments:
sumx1 = np.zeros((ny,nx))
sumx2 = np.zeros((ny,nx))
sumx3 = np.zeros((ny,nx))
sumx4 = np.zeros((ny,nx))

tag = epoch
count = 0
while (tag <= end):
  if (count % 30 == 0):
    print(tag, flush=True)

  tclim = old_climo(epoch, tag)

# Get the day's data:
  fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
  tmpnc = nc.Dataset(fbase + fname)
  sst = tmpnc.variables['sst'][0,0,:,:]
  tmpnc.close()

# Accumulate moments:
  tsst = copy.deepcopy(sst)
  tsst -= tclim

  sumx1 += tsst
  sumx2 += (tsst*tsst)
  sumx3 += (tsst*tsst*tsst)
  sumx4 += (tsst*tsst)*(tsst*tsst)
  writeout(tsst, nx, ny, lats, lons, tag)

  del tclim, tsst
  count += 1   # number of days' data
  tag   += dt
#------------------------------------------------
indices = fmask.nonzero()
applymask(sumx1, indices)
applymask(sumx2, indices)
applymask(sumx3, indices)
applymask(sumx4, indices)
# orthog1
# orthog2

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )
mean = sumx1 / count
print("mean", mean.max(), mean.min() )


# ---- .nc encoding --------------------------------------------------

name = "second_pass.nc"

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('mean', dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)

foroutput.addvar('mask', dtype = fmask.dtype)
foroutput.encodevar(fmask, 'mask')

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(mean,  'mean')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')

print("number of days = ",count)
foroutput.encodescalar(count, 'days')

foroutput.close()
#------------------ End of second pass --------------------------
