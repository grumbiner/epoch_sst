from math import *
import numpy as np
import numpy.ma as ma

import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt

from functions import *

matplotlib.use('Agg')

#=================================================================
# given the mean values of sumx1, sumx2, sumx3, sumx4, compute
#    mean, sd, skew, kurtosis

nx = 1440
ny =  720

dset = nc.Dataset("second_pass.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
fmask = dset.variables['mask'][:,:]
sumx1 = dset.variables['sumx1'][:,:]
sumx2 = dset.variables['sumx2'][:,:]
sumx3 = dset.variables['sumx3'][:,:]
sumx4 = dset.variables['sumx4'][:,:]
mean = dset.variables['mean'][:,:]
sumxn = dset.variables['sumxn'][:,:]

# scalars saved in global attributes
days = getattr(dset, 'days')
sumn = getattr(dset, 'sumn')
sumn2 = getattr(dset, 'sumn2')

print("mean nino3.4 ",sumn/days)

mask = ma.masked_array(fmask > 0)
indices = mask.nonzero()
print("fmask > 0 ",len(indices[0]) )

#----------------------------------------------------------------------
zeros = np.zeros((ny,nx))
proj  = ccrs.PlateCarree()

#----------------------------------------------------------------------

#------------- ---------------------------------------------
colors = matplotlib.colormaps.get_cmap('bwr')

print("\n\nmean")
bins = find_bins(mean, 34)
bins = np.linspace(-8,8,65)
#RG: show histogram too
show(bins, mean, "mean", "mean")

colors = matplotlib.colormaps.get_cmap('bwr')
bins = find_bins(sumx1, 32)
show(bins, sumx1, "s1", "s1", cmap = colors)
bins = find_bins(sumx2, 32)
show(bins, sumx2, "s2", "s2", cmap = colors)
bins = find_bins(sumx3, 32)
show(bins, sumx3, "s3", "s3", cmap = colors)
bins = find_bins(sumx4, 32)
show(bins, sumx4, "s4", "s4", cmap = colors)

var = (sumx2 - sumx1*sumx1/days)/days
bins = find_bins(var, 32)
show(bins, var, "var", "var", cmap = colors)

sdev = np.sqrt(var)
#bins = find_bins(sdev, 32)
bins = [0., 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5 ]
show(bins, sdev, "sdev", "sdev", cmap = colors)

ninocorr = sumxn / sqrt(sumn2) / np.sqrt(sumx2)
print("ninocorr ",ninocorr.max(), ninocorr.min(), ninocorr.mean() )
bins = np.linspace(-1,1,17)
show(bins, ninocorr, "nino", "nino", cmap = colors)
 
gram = sumxn / sumn2
print("gram ", gram.max(), gram.min(), gram.mean() )
gram = np.maximum(gram, -6.0)
#bins = find_bins(gram, 33)
bins = np.linspace(-6., 10., 33)
show(bins, gram, "gram", "gram")

cmap = matplotlib.colormaps.get_cmap('Grays')
bins = np.linspace(-0., 10., 23)
show(bins, gram, "gram2", "gram2", cmap = cmap)

cmap = matplotlib.colormaps.get_cmap('gist_gray')
bins = np.linspace(-6., 0., 13)
show(bins, gram, "gram3", "gram3", cmap = cmap)
