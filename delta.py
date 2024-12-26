import sys
from math import *
import numpy as np
import netCDF4 as nc

d1 = sys.argv[1]
d2 = sys.argv[2]
first  = nc.Dataset(d1+"/first_pass.nc","r")
second = nc.Dataset(d2+"/first_pass.nc", "r")
lons = first.variables['lon'][:]
lats = first.variables['lat'][:]
nx = len(lons)
ny = len(lats)

field = sys.argv[3]

x = first.variables[field][:,:]
y = second.variables[field][:,:]

delta = (y-x)
print(delta.max(), delta.min(), delta.mean() )
# For phase differences
print("nx, ny ",nx, ny)
delta *= 180./pi
for k in range(0,ny):
  for l in range(0,nx):
    if (delta[k,l] < -180.): delta[k,l] += 360.
    if (delta[k,l] >  180.): delta[k,l] -= 360.

dmin = -180.
dmax =  180.
space = 15.
#dmin = -0.75
#dmax =  2.00
#space = 0.25
#dmin = -8.
#dmax =  6.
#space = 1.
#dmin  = -0.5*10592
#dmax  =  2.0*10592
#space =  0.1*10592
bins = np.linspace(dmin,dmax,int((dmax-dmin)/space+1) )
hist, binedges =  np.histogram(delta, bins = bins)
#hist, binedges =  np.histogram(delta)
print(hist, hist.sum() )
print(binedges,"\n\n")

import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')
#=================================================================
def show(bins, x, title, fbase):
  bounds = np.array(bins)
  norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

  fig = plt.figure(figsize=(12, 9))
  ax  = fig.add_subplot(1, 1, 1, projection = proj)
  ax.coastlines(resolution='10m')
  ax.gridlines()

  cs = ax.pcolormesh(lons, lats, x, norm = norm, transform=ccrs.PlateCarree() )
  cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
  cbarlabel = '%s' % title
  cb.set_label(cbarlabel, fontsize=12)

  plt.savefig(fbase+".png")
  plt.close()

  hist, binedges =  np.histogram(x, bins = bins)
  print(hist, hist.sum() )
  #debug: print(binedges)


#------------------------------
proj  = ccrs.PlateCarree()
colors = matplotlib.colormaps.get_cmap('terrain')

#tmax
#bins = [ -6.0, -2.5, -1.0, -0.5, -0.25, 0, 0.25, 0.5, 1.0, 2.5, 8.0 ]
#tmin
#bins = [ -7.5, -2.5, -1.0, -0.5, -0.25, 0, 0.25, 0.5, 1.0, 2.5, 7.5 ]
# 1,2 cpy:
#bins = [ -0.75, -0.5, -0.3, -0.2,  -0.1, -0.05, 0, 0.05, 0.1, 0.20, 0.3, 0.5, 1, 2 ]
#bins = [ -0.75, -0.5, -0.3, -0.2,  -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.20, 0.3, 0.5, 1, 2 ]
# 1cpy phase
bins = [ -180., -90., -30., -15., -10., -5., 0, 5, 10, 15, 30, 90, 180 ]

bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, delta, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % (field)
cb.set_label(cbarlabel, fontsize=12)

#plt.show()
plt.savefig(field+".png")
plt.close()

hist, binedges =  np.histogram(delta, bins = bins)
print(hist, hist.sum() )
print(binedges)
