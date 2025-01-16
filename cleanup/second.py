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
#loy = 27.322   #days, lunar sidereal month
#loy = 29.53    #days, lunar synodic month

nfreq = 6
omega = np.zeros((nfreq))
omega[0] = 2.*pi/loy
omega[1] = 2.*pi/loy*2
omega[2] = 2.*pi/loy*3
omega[3] = 2.*pi/loy*4
omega[4] = 2.*pi/loy*5
omega[5] = 2.*pi/loy*6

#-------------------------------------------------
# Defining the quarter degree grid
nx = 1440
ny = 720

# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"
start = datetime.datetime(1981,9,1)
#ops: 
end = datetime.datetime(2010,8,31)

#debug: end = datetime.datetime(1981,9,30)
#debug: end = datetime.datetime(1982,8,31)
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

# Accumulate harmonics:
    tmp = copy.deepcopy(sst)
    for j in range(0, nfreq):
      hsum1[:,:,j] += tmp * cos(omega[j]*days) 
      hsum2[:,:,j] += tmp * sin(omega[j]*days) 

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
lda = 2*nfreq
coeff = np.zeros((lda, lda))
harmsums = np.zeros((ny, nx, nfreq*2))
alpha    = np.zeros((ny, nx, nfreq))
beta     = np.zeros((ny, nx, nfreq))

harmonic_coeffs(coeff, omega, days, nfreq)

for j in range(0, nfreq):
  harmsums[:,:,2*j  ] = hsum1[:,:,j ]
  harmsums[:,:,2*j+1] = hsum2[:,:,j ]

harmonic_solve(coeff, harmsums, alpha, beta, nfreq)

#------------------------------------------------
#RG: write out mean, max, min to save file
mask =  ma.masked_array(sumx1 < -900.*days)
indices = mask.nonzero()

applymask(mask, sumx1, indices)
applymask(mask, sumx2, indices)
applymask(mask, sumx3, indices)
applymask(mask, sumx4, indices)
applymask(mask, alpha, indices)
applymask(mask, beta , indices)

var = (sumx2 - sumx1*sumx1/days)/days
applymask(mask, var , indices)

print("sumx1", sumx1.max(), sumx1.min() )
print("sumx2", sumx2.max(), sumx2.min() )
print("sumx3", sumx3.max(), sumx3.min() )
print("sumx4", sumx4.max(), sumx4.min() )
print("tmax", tmax.max() , tmax.min() )
print("tmin", tmin.max() , tmin.min() )
print("alpha", alpha.max(), alpha.min(), alpha.mean() )
print("beta ", beta.max(), beta.min(), beta.mean() )
print("var  ", var.max(), var.min(), var.mean() )

ampls = np.zeros((ny, nx, nfreq))
phase = np.zeros((ny, nx, nfreq))

for j in range(0, 6):
  #debug: print(j, "alpha", alpha[:,:,j].max(), alpha[:,:,j].min(), alpha[:,:,j].mean() )
  #debug: print(j, "beta ", beta[:,:,j].max(), beta[:,:,j].min(), beta[:,:,j].mean() )
  ampls[:,:,j] = np.sqrt(alpha[:,:,j]**2 + beta[:,:,j]**2)
  phase[:,:,j] = np.arctan2(beta[:,:,j], alpha[:,:,j])
  print(j, "ampls", ampls[:,:,j].max(), ampls[:,:,j].min(), ampls[:,:,j].mean() )
  print(j, "phase", phase[:,:,j].max() )

  
harm_var = 0.5*ampls*ampls
tot_var  = np.zeros((ny,nx))
for i in range(0,nfreq):
  tot_var += harm_var[:,:,i]

frac_var   = np.zeros((ny,nx))
frac_var1  = np.zeros((ny,nx))
for j in range(0, ny):
  for i in range(0, nx):
    if (var[j,i] > 0):
      frac_var[j,i] = tot_var[j,i] / var[j,i]
      frac_var1[j,i] = harm_var[j,i,0] / var[j,i]

print("total harmonic fraction of variance ",frac_var.max(), frac_var.min(), frac_var.mean() )
print("first harmonic fraction of variance ",frac_var1.max(), frac_var1.min(), frac_var1.mean() )

#-------------------------------------------------
import ncoutput

name = "first_pass.nc"

foroutput = ncoutput.ncoutput(nx, ny, lats, lons, name)
foroutput.ncoutput(name)
foroutput.addvar('sumx1', dtype = sumx1.dtype)
foroutput.addvar('sumx2', dtype = sumx2.dtype)
foroutput.addvar('sumx3', dtype = sumx3.dtype)
foroutput.addvar('sumx4', dtype = sumx4.dtype)
foroutput.addvar('tmax', dtype = tmax.dtype)
foroutput.addvar('tmin', dtype = tmin.dtype)
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

foroutput.encodevar(sumx1, 'sumx1')
foroutput.encodevar(sumx2, 'sumx2')
foroutput.encodevar(sumx3, 'sumx3')
foroutput.encodevar(sumx4, 'sumx4')
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

print("days = ",days)
foroutput.encodescalar(days, 'days')

foroutput.close()
#------------------ End of first pass --------------------------

"""
Second pass:
consider harmonics and how much variance they explain (doable in first)
for a given harmonic, % = (1/2*ampl**2) / sigma**2
May sum them
"""
