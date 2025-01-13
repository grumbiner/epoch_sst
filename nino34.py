from math import *
import os
import datetime
import copy

import numpy as np
import numpy.ma as ma
import netCDF4

#-------------------------------------------------
from harmonic_grid import * 

# Define harmonic frequencies
loy = 365.2422 #days, tropical year
lom = 29.53    #days, lunar synodic month
#lom = 27.322    = 1/(1/lom + 1/loy)

nfreq = 6
omega = np.zeros((nfreq))
omega[0] = 2.*pi/loy
omega[1] = 2.*pi/loy*2
omega[2] = 2.*pi/loy*3
omega[3] = 2.*pi/loy*4
omega[4] = 2.*pi/loy*5
omega[5] = 2.*pi/loy*6

#-------------------------------------------------
# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"

# Defining the quarter degree grid
nx = 1440
ny = 720

start = datetime.datetime(1981,9,1)
end = datetime.datetime(2024,12,8)

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

# for trend
sumt  = np.zeros((ny,nx))
sumxt = np.zeros((ny,nx))
sumt2 = np.zeros((ny,nx))

# for harmonic summing
hsum1    = np.zeros((ny, nx, nfreq))
hsum2    = np.zeros((ny, nx, nfreq))

# extrema
tmax = np.zeros((ny,nx))
tmin = np.zeros((ny,nx))
tmax.fill(-3.0)
tmin.fill(45.0)

#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = start
days = (tag - datetime.datetime(1981,9,1) ).days
days  = 0
count = 0
n0   = days
while (tag <= end ):
    if (count % 30 == 0):
      print(tag, flush=True)

# Get the day's data:
    fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
    tmpnc = netCDF4.Dataset(fbase + fname)
    sst = tmpnc.variables['sst'][0,0,:,:]
    if ( count ==  0 ):
        lons = tmpnc.variables['lon'][:]
        #debug: print(lons, "\nlen(lons)",len(lons),  flush=True)
        zlon = np.logical_and(lons > 190, lons < 240)
        #debug: print(zlon[720:960], flush=True)
        lats = tmpnc.variables['lat'][:]
        zlat = np.logical_and(lats > -5., lats < 5.)
        #debug: print(zlat[300:420], flush=True)
        nino34 = np.zeros((ny, nx),dtype=bool)
        nino34[zlat,:] = True
        nino34[:,zlon] &= True
        print(nino34)
        del zlon, zlat
    tmpnc.close()
    print(count, sst[nino34].mean() )
    #debug: del sst
    #debug: days  += 1   # days since epoch
    #debug: count += 1   # number of days' data
    #debug: tag   += dt
    #debug: continue

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

# Accumulate harmonics:
    tmp = copy.deepcopy(sst)
    for j in range(0, nfreq):
      hsum1[:,:,j] += tmp * cos(omega[j]*days) 
      hsum2[:,:,j] += tmp * sin(omega[j]*days) 

# Find extrema:
    tmax = np.fmax(tmax, sst)
    tmin = np.fmin(tmin, sst)
    
    days  += 1   # days since epoch
    count += 1   # number of days' data
    tag   += dt

#------------------------------------------------
def applymask(mask, grid, indices):
  for k in range(0, len(indices[0])):
    i = indices[1][k]
    j = indices[0][k]
    grid[j,i] = 0.   

#------------------------------------------------
lda = 2*nfreq
coeff = np.zeros((lda, lda))
harmsums = np.zeros((ny, nx, nfreq*2))
alpha    = np.zeros((ny, nx, nfreq))
beta     = np.zeros((ny, nx, nfreq))

#harmonic_coeffs(coeff, omega, days, nfreq)
harmonic_coeffs(coeff, omega, count, nfreq, n0 = n0) # rg: probably need n0 here, too.

mean = sumx1/count
for j in range(0, nfreq):
  harmsums[:,:,2*j  ]  = hsum1[:,:,j ] 
  harmsums[:,:,2*j+1]  = hsum2[:,:,j ] 

# solve for harmonics-only
harmonic_solve(coeff, harmsums, alpha, beta, nfreq)

#trend_harmonic_solve(coeff2, sumx1, sumxt, sumycos, sumysin,
#                      coeff, harmsums, alpha, beta, nfreq)

#------------------------------------------------
#RG: write out mean, max, min to save file
mask =  ma.masked_array(sumx1 < -900.*days)
indices = mask.nonzero()

applymask(mask, sumx1, indices)
applymask(mask, sumx2, indices)
applymask(mask, sumx3, indices)
applymask(mask, sumx4, indices)
applymask(mask, sumxt, indices)
applymask(mask, sumt, indices)
applymask(mask, sumt2, indices)
applymask(mask, alpha, indices)
applymask(mask, beta , indices)

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )
print("tmax", tmax.max() , tmax.min() )
print("tmin", tmin.max() , tmin.min() )
print("alpha", alpha.max(), alpha.min(), alpha.mean() )
print("beta ", beta.max(), beta.min(), beta.mean() )

ampls = np.zeros((ny, nx, nfreq))
phase = np.zeros((ny, nx, nfreq))

for j in range(0, nfreq):
  #debug: print(j, "alpha", alpha[:,:,j].max(), alpha[:,:,j].min(), alpha[:,:,j].mean() )
  #debug: print(j, "beta ", beta[:,:,j].max(), beta[:,:,j].min(), beta[:,:,j].mean() )
  ampls[:,:,j] = np.sqrt(alpha[:,:,j]**2 + beta[:,:,j]**2)
  phase[:,:,j] = np.arctan2(beta[:,:,j], alpha[:,:,j])*180./pi
  print(j, "ampls", ampls[:,:,j].max(), ampls[:,:,j].min(), ampls[:,:,j].mean() )
  #debug: print(j, "phase", phase[:,:,j].max() )
  for k in range(0,ny):
    for l in range(0,nx):
      if (phase[k,l,j] < -180.): phase[k,l,j] += 360.
      if (phase[k,l,j] <  180.): phase[k,l,j] -= 360.
  
#-------------------------------------------------
import ncoutput

name = "first_pass.nc"

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('mean',  dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)
foroutput.addvar('sumt',  dtype = sumt.dtype)
foroutput.addvar('sumxt', dtype = sumxt.dtype)
foroutput.addvar('sumt2', dtype = sumt2.dtype)
foroutput.addvar('intercept', dtype = sumx1.dtype)
foroutput.addvar('slope',    dtype = sumx1.dtype)
foroutput.addvar('correl'    dtype = sumx1.dtype)
foroutput.addvar('tstat',    dtype = sumx1.dtype)
foroutput.addvar('tmax',     dtype = tmax.dtype)
foroutput.addvar('tmin',     dtype = tmin.dtype)
foroutput.addvar('cpy1_amp', dtype = ampls.dtype)
foroutput.addvar('cpy1_pha', dtype = phase.dtype)
foroutput.addvar('cpy2_amp', dtype = ampls.dtype)
foroutput.addvar('cpy2_pha', dtype = phase.dtype)
foroutput.addvar('cpy3_amp', dtype = ampls.dtype)
foroutput.addvar('cpy3_pha', dtype = phase.dtype)
foroutput.addvar('cpy4_amp', dtype = ampls.dtype)
foroutput.addvar('cpy4_pha', dtype = phase.dtype)
foroutput.addvar('cpy5_amp', dtype = ampls.dtype)
foroutput.addvar('cpy5_pha', dtype = phase.dtype)
foroutput.addvar('cpy6_amp', dtype = ampls.dtype)
foroutput.addvar('cpy6_pha', dtype = phase.dtype)

mean = sumx1/count
applymask(mask, mean, indices)
applymask(mask, sumx1, indices)
applymask(mask, sumxt, indices)
applymask(mask, sumt, indices)
applymask(mask, sumt2, indices)

tmpt = days*sumt2 - sumt*sumt
tmpx = days*sumx2 - sumx1*sumx1
print(tmpt.min(), tmpt.max(), tmpt.mean() )
tcount = 0
xcount = 0
tlim = tmpt.max()
xlim = tmpx.max()
for j in range (0, ny):
  for i in range (0, nx):
    if (tmpt[j,i] == 0): 
      tmpt[j,i] = tlim
      tcount += 1
    if (tmpx[j,i] == 0): 
      tmpx[j,i] = xlim
      xcount += 1
print("count of zero tmpt",tcount, "tmpx",xcount) 
print(tmpt.min(), tmpt.max(), tmpt.mean() )
slope = (days*sumxt - sumx1*sumt) / tmpt #RG: should be from trend solver
applymask(mask, slope, indices)

intercept = (sumx1/days - slope*sumt/days) #RG: ditto
applymask(mask, slope, indices)

correl = (days*sumxt - sumx1*sumt ) / (np.sqrt(tmpt) * np.sqrt(tmpx))
applymask(mask, correl, indices)

tstat = correl*sqrt(days) / (1. - correl*correl)
applymask(mask, tstat, indices)

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(mean, 'mean')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')
foroutput.encodevar(sumt, 'sumt')
foroutput.encodevar(sumxt, 'sumxt')
foroutput.encodevar(sumt2, 'sumt2')
foroutput.encodevar(slope, 'slope')
foroutput.encodevar(intercept, 'intercept')
foroutput.encodevar(correl, 'correl')
foroutput.encodevar(tstat, 'tstat')
foroutput.encodevar(tmin, 'tmin')
foroutput.encodevar(tmax, 'tmax')
foroutput.encodevar(ampls[:,:,0], 'cpy1_amp')
foroutput.encodevar(phase[:,:,0], 'cpy1_pha')
foroutput.encodevar(ampls[:,:,1], 'cpy2_amp')
foroutput.encodevar(phase[:,:,1], 'cpy2_pha')
foroutput.encodevar(ampls[:,:,2], 'cpy3_amp')
foroutput.encodevar(phase[:,:,2], 'cpy3_pha')
foroutput.encodevar(ampls[:,:,3], 'cpy4_amp')
foroutput.encodevar(phase[:,:,3], 'cpy4_pha')
foroutput.encodevar(ampls[:,:,4], 'cpy5_amp')
foroutput.encodevar(phase[:,:,4], 'cpy5_pha')
foroutput.encodevar(ampls[:,:,5], 'cpy6_amp')
foroutput.encodevar(phase[:,:,5], 'cpy6_pha')

tmask = np.zeros((ny,nx))
for k in range(0, len(indices[0]) ):
    i = indices[1][k]
    j = indices[0][k]
    tmask[j,i] = 1.0
foroutput.addvar('mask', dtype = tmask.dtype)
foroutput.encodevar(tmask, 'mask')

print("number of days = ",count)
foroutput.encodescalar(days, 'days')

foroutput.close()
#------------------ End of first pass --------------------------
