''' legacy -- plot mean from a.nc and print its histogram '''

#from math import *
import numpy as np

import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt

from functions import *

matplotlib.use('Agg')

#=================================================================

nx = 1440
ny =  720

dset = nc.Dataset("a.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
mean = dset.variables['mean'][:,:]
print("mean ",mean.max(), mean.min(), mean.mean() )

# scalars saved in global attributes
days = getattr(dset, 'days')

#----------------------------------------------------------------------
zeros = np.zeros((ny,nx))
proj  = ccrs.PlateCarree()

#----------------------------------------------------------------------
#RG: binning of residuals
bins = np.linspace(-6.,6,25)

#RG: plot geographic map in ~deciles ----------------

#colors = matplotlib.colormaps.get_cmap('terrain')
#colors = matplotlib.colormaps.get_cmap('seismic')
colors = matplotlib.colormaps.get_cmap('bwr')

# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, mean, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("Residual after Niño 3.4 extraction")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("finalres.png")
plt.close()

hist, binedges =  np.histogram(mean, bins = bins)
print(hist, hist.sum() )
print(binedges)

#----------------------------------------------------------------------
