from math import *
import numpy as np
import numpy.ma as ma

import netCDF4 as nc

# given the mean values of sumx1, sumx2, sumx3, sumx4, sumxt, sumt2, sumt, compute
#    mean, sd, skew, kurtosis, and trend

nx = 1440
ny =  720

dset = nc.Dataset("first_pass.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
sumx1 = dset.variables['sumx1'][:,:]
sumx2 = dset.variables['sumx2'][:,:]
sumx3 = dset.variables['sumx3'][:,:]
sumx4 = dset.variables['sumx4'][:,:]
sumxt = dset.variables['sumxt'][:,:]
fmask = dset.variables['mask'][:,:]

# scalars saved in global attributes
days = getattr(dset, 'days')
sumt = getattr(dset, 'sumt')
sumt2 = getattr(dset, 'sumt2')

# Proceed to statistical summary
mean1 = sumx1 / days
print("mean1 ",mean1.max(), mean1.min() )
print("sumx2 ",sumx2.max(), sumx2.min() )
print("sumx3 ",sumx3.max(), sumx3.min() )
print("sumx4 ",sumx4.max(), sumx4.min() )
print("sumxt ",sumxt.max(), sumxt.min() )
print("days ", days)
meant = sumt / days
meant2 = sumt2 / days
print("sumt ", sumt)
print("sumt2 ", sumt2)

mask = ma.masked_array(fmask > 0)
indices = mask.nonzero()
print("fmask > 0 ",len(indices[0]) )

sumx2 = ma.masked_array(sumx2, mask)
sumxt = ma.masked_array(sumxt, mask)

#----------------------------------------------------------------------
#avg = (masked array)(mean, mask) 
#var = 
#skew = 
#kurtosis = 

zeros     = np.zeros((ny,nx))
slope     = np.zeros((ny,nx))
intercept = np.zeros((ny,nx))
r         = np.zeros((ny,nx))

n = days

var = ma.masked_array(zeros, mask)
var = (sumx2 - sumx1*sumx1/n)/n
tmask = ma.masked_array(var <= 0)
indices = tmask.nonzero()
print("var nonpositive ",len(indices[0]) )
del tmask
hist,binedges = np.histogram(var, bins = [0, 0.25, 0.5, 1, 2, 3, 4, 6, 9, 16, 25, 100] )
print("variance\n",hist)
print(binedges)
print(hist.sum() )

#rg: slope = ma.masked_array( (sumxt - sumx1*sumt/n) / var, mask)
slope = (n*sumxt - sumx1*sumt) / (n*sumx2-sumx1*sumx1)
print("slope ",slope.max(), slope.min() )

hist,binedges = np.histogram(slope, bins = [-1500, -10, -1, 0, 1, 10, 700])
print(hist, hist.sum() )
print(binedges)


smask = ma.masked_array(np.abs(slope) > 100) 
indices = smask.nonzero()
print(len(indices[0]))
for k in range(0, len(indices[0])):
  i = indices[1][k]
  j = indices[0][k]
  print(i,j, slope[j,i], var[j,i])

exit(0)

sdyy = sqrt(sumyy - days*ym*ym)

for k in range(0, len(indices[0]) ):
    if (k%300 == 0):
      print("k = ",k)
    i = indices[1][k]
    j = indices[0][k]
    if (sigmax2[j,i] != 0):
      slope = (sumxt[j,i] - n*xm[j,i]*ym) / sigmax2[j,i]
      r = (sumxt[j,i] - n*xm[j,i]*ym)/ np.sqrt(sigmax2[j,i])/ sdyy
    intercept = (ym - slope*xm[j,i])
 
print("slope    ", slope.max(), slope.min() )
print("intercept", intercept.max(), intercept.min() )
print("r        ", r.max(), r.min() )


