import sys
from math import *
import numpy as np
import numpy.ma as ma

import netCDF4 as nc

from functions import *

#=================================================================
# given the mean values of sumx1, sumx2, sumx3, sumx4, compute
#    mean, sd, skew, kurtosis

nx = 1440
ny =  720

dset = nc.Dataset(sys.argv[1], "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
fmask = dset.variables['mask'][:,:]
sumx1 = dset.variables['sumx1'][:,:]
sumx2 = dset.variables['sumx2'][:,:]
sumx3 = dset.variables['sumx3'][:,:]
sumx4 = dset.variables['sumx4'][:,:]
mean = dset.variables['mean'][:,:]

# scalars saved in global attributes
#days = getattr(dset, 'days')
days=3650
#days = 10

mask = ma.masked_array(fmask > 0)
indices = mask.nonzero()
print("fmask > 0 ",len(indices[0]) )

#----------------------------------------------------------------------
zeros = np.zeros((ny,nx))
proj  = ccrs.PlateCarree()

#----------------------------------------------------------------------

#------------- ---------------------------------------------
colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')
#colors = matplotlib.colormaps.get_cmap('viridis')

bins = find_bins(mean, 36)
bins = np.linspace(-6., 6., 25)
#RG: plot histogram too
show(bins, lons, lats, mean, "decadal mean residual", "mean", cmap = colors, proj = proj)

sumx1 /= days
bins = find_bins(sumx1, 36)
#bins = np.linspace(-5,5,21)
show(bins, lons, lats, sumx1, "s1", "s1", cmap = colors)

sumx2 /= days
bins = find_bins(sumx2, 32)
#bins = [0., 0.25, 1.0, 2.25, 4, 9, 16, 25, 36, 49. ]
show(bins, lons, lats, sumx2, "s2", "s2", cmap = colors)

sumx3 /= days
bins = find_bins(sumx3, 32)
#bins = [-512, -343, -216, -125, -64, -27, -8, -1, 0, 1, 8, 27, 64, 125, 216, 343, 512 ]
show(bins, lons, lats, sumx3, "s3", "s3", cmap = colors)

sumx4 /= days
bins = find_bins(sumx4, 32)
#bins = [0, 1, 16, 81, 625, 1296, 2401, 4096, 10000 ]
show(bins, lons, lats, sumx4, "s4", "s4", cmap = colors)

var = (sumx2*days - mean*mean)
var = np.maximum(var, 0)
print('var ',var.min(), var.max() )
bins = find_bins(var, 32)
#bins = [0, 0.25, 1, 2.25, 4, 6.25, 9, 12.25, 16, 20.25, 25, 30.25, 36 ]
show(bins, lons, lats, var, "var", "var", cmap = colors )

sdev = np.sqrt(var)
bins = find_bins(sdev, 32)
#bins = [0., 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5 ]
show(bins, lons, lats, sdev, "sdev", "sdev", cmap = colors)

