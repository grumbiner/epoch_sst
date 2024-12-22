from math import *
import os
import datetime
import copy

import numpy as np
import numpy.ma as ma
import netCDF4

#-------------------------------------------------
# Defining the quarter degree grid
nx = 1440
ny = 720

# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"
start = datetime.datetime(1981,9,1)
#debug: end = datetime.datetime(1981,9,30)
#debug: end = datetime.datetime(1982,8,31)
#ops: 
end = datetime.datetime(2010,8,31)

dt = datetime.timedelta(1)
tag = start

# quick check that all data files exist
errcount = 0
while (tag <= end and errcount < 90 ):
    fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
    if (not os.path.exists(fbase+fname)):
        print("no file for ",fbase+fname)
        errcount += 1
    tag += dt

#debug: print("error count in running over target period:",errcount)
if (errcount != 0):
    print("find the missing data files!")
    exit(1)

#-------------------------------------------------
# Initialize files for accumulations
sst = np.zeros((ny,nx)) # temporary file for reading in data
tmp = np.zeros((ny,nx))

# for accumulating moments:
sumx1 = np.zeros((ny,nx))
sumx2 = np.zeros((ny,nx))
sumx3 = np.zeros((ny,nx))
sumx4 = np.zeros((ny,nx))
#debug: print('dtype for sumx1 ',sumx1.dtype, flush=True)

sumt  = np.zeros((ny,nx))
sumxt = np.zeros((ny,nx))
sumt2 = np.zeros((ny,nx))

# extrema
tmax = np.zeros((ny,nx))
tmin = np.zeros((ny,nx))
tmax.fill(-3.0)
tmin.fill(45.0)

#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = start
days = 0
while (tag <= end ):
    if (days % 30 == 0):
      print(tag, flush=True)

# Get the day's data:
    fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
    tmpnc = netCDF4.Dataset(fbase + fname)
    sst = tmpnc.variables['sst'][0,0,:,:]
    if ( days == 0 ):
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

# Accumulate trend:
    tmp    = copy.deepcopy(sst)
    sumt  += float(days)
    sumt2 += float(days*days)
    sumxt += tmp*float(days)

# Find extrema:
    tmax = np.fmax(tmax, sst)
    tmin = np.fmin(tmin, sst)
    
    days += 1
    tag += dt

#------------------------------------------------
def applymask(mask, grid, indices):
  for k in range(0, len(indices[0])):
    i = indices[1][k]
    j = indices[0][k]
    grid[j,i] = 0.   

#------------------------------------------------
#RG: write out mean, max, min to save file
mask =  ma.masked_array(sumx1 < -900.*days)
indices = mask.nonzero()

applymask(mask, sumx1, indices)
applymask(mask, sumx2, indices)
applymask(mask, sumx3, indices)
applymask(mask, sumx4, indices)
applymask(mask, sumt , indices)
applymask(mask, sumxt, indices)
applymask(mask, sumt2, indices)

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )
print("sumt ", sumt.max(), sumt.min() )
print("sumxt", sumxt.max(), sumxt.min() )
print("sumt2", sumt2.max(), sumt2.min() )
print("tmax", tmax.max() , tmax.min() )
print("tmin", tmin.max() , tmin.min() )

#-------------------------------------------------
import ncoutput

name = "first_pass.nc"

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)
foroutput.addvar('sumt', dtype = sumt.dtype)
foroutput.addvar('sumxt', dtype = sumxt.dtype)
foroutput.addvar('sumt2', dtype = sumt2.dtype)
foroutput.addvar('tmax', dtype = tmax.dtype)
foroutput.addvar('tmin', dtype = tmin.dtype)

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')
foroutput.encodevar(sumt, 'sumt')
foroutput.encodevar(sumxt, 'sumxt')
foroutput.encodevar(sumt2, 'sumt2')
foroutput.encodevar(tmin, 'tmin')
foroutput.encodevar(tmax, 'tmax')

tmask = np.zeros((ny,nx))
for k in range(0, len(indices[0]) ):
    i = indices[1][k]
    j = indices[0][k]
    tmask[j,i] = 1.0
foroutput.addvar('mask', dtype = tmask.dtype)
foroutput.encodevar(tmask, 'mask')

print("days = ",days)
foroutput.encodescalar(days, 'days')

foroutput.close()
#------------------ End of first pass --------------------------
