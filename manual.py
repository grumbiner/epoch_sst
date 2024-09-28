from math import *
import os
import datetime
import copy

import numpy as np
import numpy.ma as ma
import netCDF4

import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt

# Defining the quarter degree grid
nx = 1440
ny = 720

# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"
start = datetime.datetime(1981,9,1)
#ops: end = datetime.datetime(2010,8,31)
#debug: 
end = datetime.datetime(1981,9,30)
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

# Initialize files for accumulations
sst = np.zeros((ny,nx)) # temporary file for reading in data
tmp = np.zeros((ny,nx))
t0 = np.zeros((ny,nx)) # first sst field, which will be subtracted

# for accumulating moments:
sumx1 = np.zeros((ny,nx))
sumx2 = np.zeros((ny,nx))
sumx3 = np.zeros((ny,nx))
sumx4 = np.zeros((ny,nx))
#debug: print('dtype for sumx1 ',sumx1.dtype, flush=True)

# for trend computation
sumt = 0
sumxt = np.zeros((ny,nx))
sumt2 = np.zeros((ny,nx))
days = 0 # starting time = 0 at first date of the record


# extrema
tmax = np.zeros((ny,nx))
tmin = np.zeros((ny,nx))
tmax.fill(-3.0)
tmin.fill(45.0)
# numpy.fmax(f1,f2)

# for accumulating harmonics

# For masking:

#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = start
#debug: end = datetime.datetime(1981,9,30)
days = 0
while (tag <= end ):
    if (days % 30 == 0):
      print(tag)

# Get the day's data:
    fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
    tmpnc = netCDF4.Dataset(fbase + fname)
    sst = tmpnc.variables['sst'][0,0,:,:]
    if ( days == 0 ):
        t0 = copy.deepcopy(sst)
        lons = tmpnc.variables['lon'][:]
        lats = tmpnc.variables['lat'][:]
    tmpnc.close()

    # this messes up because of masking: sst -= t0

# Accumulate moments:
    tmp = copy.deepcopy(sst)
    #debug: print("moment1", tag, tmp.max(), tmp.min() )
    sumx1 += tmp
    tmp *= sst
    #debug: print("moment2",tag, tmp.max(), tmp.min() )
    sumx2 += tmp
    tmp *= sst
    sumx3 += tmp
    tmp *= sst
    sumx4 += tmp
    #debug: print("moment4",sumx4.max(), sumx4.min() )

# Accumulate for trends:
    tmp = copy.deepcopy(sst)
    #debug: print("trends ",tmp.max(), tmp.min() )
    tmp *= days
    sumxt += tmp
    tmp.fill(days*days)
    #debug: print(tmp.max(), tmp.min() )
    sumt2 += tmp
    sumt += days

# Find extrema:
    tmax = np.fmax(tmax, sst)
    tmin = np.fmin(tmin, sst)
    #debug: print("tmax ",tmax.max(), tmax.min() )
    #debug: print("tmin ",tmin.max(), tmin.min() )
    
    days += 1
    tag += dt


#------------------------------------------------
def findmask(mask, sst):
    #debug: print("mask ",mask.max(), mask.min(), sst.max(), sst.min(), flush = True)
    tmask =  ma.masked_array(sst < -900.*30)
    mask = ma.mask_or(mask, tmask)

    indices = tmask.nonzero()
    #debug: print("later-t ",  len(indices[0]) )
    del tmask

    indices = mask.nonzero()
    #debug: print("later ",  len(indices[0]) )
    return indices

def applymask(mask, grid, indices):
  #debug: print("applymask", indices[0])
  for k in range(0, len(indices[0])):
    i = indices[1][k]
    j = indices[0][k]
    grid[j,i] = 0.   

print("tmax", tmax.max() , tmax.min() )
print("tmin", tmin.max() , tmin.min() )

#RG: write out mean, max, min to save file
mask =  ma.masked_array(sumx1 < -900.*days)
indices = findmask(mask, sumx1)
indices = findmask(mask, sumx2)
indices = findmask(mask, sumx3)
indices = findmask(mask, sumx4)

applymask(mask, sumx1, indices)
applymask(mask, sumx2, indices)
applymask(mask, sumx3, indices)
applymask(mask, sumx4, indices)
applymask(mask, sumxt, indices)
applymask(mask, sumt2, indices)

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )

#
print("sumxt", sumxt.max(), sumxt.min() )
print("sumt2", sumt2.max(), sumt2.min() )
print("sumt",sumt, sumt/days)

#-------------------------------------------------
import ncoutput

#-------------------------------------------------
name = "first_pass.nc"
sumx1 /= float(days)
sumx2 /= float(days)
sumx3 /= float(days)
sumx4 /= float(days)
sumxt /= float(days)
sumt2 /= float(days)

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)
foroutput.addvar('sumxt', dtype = sumxt.dtype)
foroutput.addvar('sumt2', dtype = sumt2.dtype)

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')
foroutput.encodevar(sumxt, 'sumxt')
foroutput.encodevar(sumt2, 'sumt2')

foroutput.close()

print("days = ",days, "sumt = ",sumt, sumt/days)

