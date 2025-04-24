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

dset = nc.Dataset("first_pass.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
fmask = dset.variables['mask'][:,:]
sumx1 = dset.variables['sumx1'][:,:]
sumx2 = dset.variables['sumx2'][:,:]
mean      = dset.variables['mean'][:,:]
slope     = dset.variables['slope'][:,:]
intercept = dset.variables['intercept'][:,:]
correl    = dset.variables['correl'][:,:]
tstat     = dset.variables['tstat'][:,:]

harm1 = dset.variables['cpy1_amp'][:,:]
harm2 = dset.variables['cpy2_amp'][:,:]
harm3 = dset.variables['cpy3_amp'][:,:]
harm4 = dset.variables['cpy4_amp'][:,:]
harm5 = dset.variables['cpy5_amp'][:,:]
harm6 = dset.variables['cpy6_amp'][:,:]
#harm7 = dset.variables['cpy7_amp'][:,:]
#harm8 = dset.variables['cpy8_amp'][:,:]
#harm9 = dset.variables['cpy9_amp'][:,:]
#harm10 = dset.variables['cpy10_amp'][:,:]
#harm11 = dset.variables['cpy11_amp'][:,:]
#harm12 = dset.variables['cpy12_amp'][:,:]

# scalars saved in global attributes
days = getattr(dset, 'days')

mask = ma.masked_array(fmask > 0)
indices = mask.nonzero()
print("fmask > 0 ",len(indices[0]) )

#----------------------------------------------------------------------
zeros = np.zeros((ny,nx))
proj  = ccrs.PlateCarree()


#----------------------------------------------------------------------
#1 CPY:

#RG: binning of annual cycle amplitudes
bins = [0, 0.01, 0.1, 0.25, 0.5, 0.75, 1., 1.5, 2, 3, 5, 7.5, 10, 15 ]

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

cs = ax.pcolormesh(lons, lats, harm1, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("first harmonic temperature")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("cpy1_geo.png")
plt.close()

hist, binedges =  np.histogram(harm1, bins = bins)
print(hist, hist.sum() )

#----------------------------------------------------------------------
#2 CPY:

#RG: binning of annual cycle amplitudes
bins = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1., 1.5, 2, 3, 3.5, 5 ]

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

cs = ax.pcolormesh(lons, lats, harm2, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("second harmonic temperature")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("cpy2_geo.png")
plt.close()

hist, binedges =  np.histogram(harm2, bins = bins)
print(hist, hist.sum() )
#----------------------------------------------------------------------
#3 CPY:

#RG: binning of harmonic cycle amplitudes
# Use same bins + colors as 2 cpy for the rest
#bins = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1., 2, 3  ]
#colors = matplotlib.colormaps.get_cmap('terrain')
#colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')

# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, harm3, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("third harmonic temperature")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("cpy3_geo.png")
plt.close()

hist, binedges =  np.histogram(harm3, bins = bins)
print(hist, hist.sum() )
#----------------------------------------------------------------------
#4 CPY:

#RG: binning of annual cycle amplitudes
#bins = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.1, 0.25, 0.5  ]
bins = np.linspace(0,3.5,14)
#colors = matplotlib.colormaps.get_cmap('terrain')
#colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')

# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, harm4, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("fourth harmonic temperature")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("cpy4_geo.png")
plt.close()

hist, binedges =  np.histogram(harm4, bins = bins)
print(hist, hist.sum() )
#----------------------------------------------------------------------
#5 CPY:

#RG: binning of harmonic amplitudes
#bins = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1., 2, 2.5  ]
#colors = matplotlib.colormaps.get_cmap('terrain')
#colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')

# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, harm5, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("fifth harmonic temperature")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("cpy5_geo.png")
plt.close()

hist, binedges =  np.histogram(harm5, bins = bins)
print(hist, hist.sum() )
#----------------------------------------------------------------------
#6 CPY:

#RG: binning of annual cycle amplitudes
#bins = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1., 2, 2.5  ]
#colors = matplotlib.colormaps.get_cmap('terrain')
#colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')

# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, harm6, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("sixth harmonic temperature")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("cpy6_geo.png")
plt.close()

hist, binedges =  np.histogram(harm6, bins = bins)
print(hist, hist.sum() )
#----------------------------------------------------------------------

#Fraction of variance explained
colors = matplotlib.colormaps.get_cmap('bwr')

var = (sumx2 - sumx1*sumx1/days)/days
tot_var  = np.zeros((ny,nx))
tot_var += harm1*harm1/2.
tot_var += harm2*harm2/2.
tot_var += harm3*harm3/2.
tot_var += harm4*harm4/2.
tot_var += harm5*harm5/2.
tot_var += harm6*harm6/2.

frac_var   = np.zeros((ny,nx))
frac_var1  = np.zeros((ny,nx))
res_var    = np.zeros((ny,nx))
for j in range(0, ny):
  for i in range(0, nx):
    if (var[j,i] > 0):
      frac_var[j,i]  = tot_var[j,i] / var[j,i]
      frac_var1[j,i] = 0.5*harm1[j,i]*harm1[j,i] / var[j,i]
hist,binedges = np.histogram(frac_var)
print("total harmonic fraction ")
print(hist, hist.sum() )

print("1st harmonic fraction ")
hist,binedges = np.histogram(frac_var1)
print(hist, hist.sum() )

# -- plot tot
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, frac_var, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("harmonic fraction of variance")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("harm_var.png")
plt.close()

# -- plot 1st harmonic 
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, frac_var1, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("first harmonic fraction of variance")
cb.set_label(cbarlabel, fontsize=12)

plt.savefig("harm1_var.png")
plt.close()

#------------- ---------------------------------------------
colors = matplotlib.colormaps.get_cmap('bwr')

print("\n\nmean")
bins = find_bins(mean, 34)
show(bins, lons, lats, mean, "mean", "mean")

print("\n\nintercept")
bins = find_bins(intercept, 36)
show(bins, lons, lats, intercept, "intercept", "intercept")

print("\n\nslope")
slope *= 365.2422*10.
print("slope stats",slope.max(), slope.min(), slope.mean() )
bins = [ -2, -1, -.5, -.2, -.1, 0, .1, .2, .5, 1.0, 2.00 ]
#bins = find_bins(slope, 64)
show(bins, lons, lats, slope, "trend (K/decade)", "slope")

print("\n\ncorrel")
#bins = find_bins(correl, 32)
bins = np.linspace(-1., 1., 33 )
show(bins, lons, lats, correl, "correl", "correl")

#bins = find_bins(tstat, 32)
print("tstat stats",tstat.max(), tstat.min(), tstat.mean() )
bins = [-400, -63.7, -6.3, -3.078, 3.078, 6.3, 63.7, 400 ]
show(bins, lons, lats, tstat, "tstat", "tstat")



res_var = (1.-frac_var1)*var
bins = [0, 0.01, 0.0625, 0.25, 0.5, 1, 4, 9, 25 ]
show(bins, lons, lats, res_var, "residual variance", "res_var")

colors = matplotlib.colormaps.get_cmap('bwr')
bins = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1., 1.5, 2, 3, 4, 5, 7.5 ]
show(bins, lons, lats, harm1, "h1", "h1", cmap = colors)
show(bins, lons, lats, harm2, "h2", "h2", cmap = colors)
show(bins, lons, lats, harm3, "h3", "h3", cmap = colors)

for i in range (0, len(bins)):
  bins[i] /= 10.
show(bins, lons, lats, harm4, "h4", "h4", cmap = colors)
show(bins, lons, lats, harm5, "h5", "h5", cmap = colors)
show(bins, lons, lats, harm6, "h6", "h6", cmap = colors)
