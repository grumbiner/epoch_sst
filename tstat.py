from math import *
import numpy as np
import numpy.ma as ma

import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt



# given the mean values of sumx1, sumx2, sumx3, sumx4, compute
#    mean, sd, skew, kurtosis

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
sumt  = dset.variables['sumt'][:,:]
sumt2 = dset.variables['sumt2'][:,:]
fmask = dset.variables['mask'][:,:]

# scalars saved in global attributes
days = getattr(dset, 'days')

# Proceed to statistical summary
mean1 = sumx1 / days
print("mean1 ",mean1.max(), mean1.min() )
print("sumx2 ",sumx2.max(), sumx2.min() )
print("sumx3 ",sumx3.max(), sumx3.min() )
print("sumx4 ",sumx4.max(), sumx4.min() )
print("days ", days)

mask = ma.masked_array(fmask > 0)
indices = mask.nonzero()
print("fmask > 0 ",len(indices[0]) )

sumx2 = ma.masked_array(sumx2, mask)

#----------------------------------------------------------------------
zeros = np.zeros((ny,nx))
proj  = ccrs.PlateCarree()


#----------------------------------------------------------------------
slope = (days*sumxt - sumx1*sumt) / (days*sumt2 - sumt*sumt)
print(slope.max(), slope.min(), slope.mean() )
intercept = (sumx1/days - slope*sumt/days)
print(intercept.max(), intercept.min(), intercept.mean() )

correl = (days*sumxt - sumx1*sumt ) / (np.sqrt(days*sumt2 - sumt*sumt) * np.sqrt(days*sumx2 - sumx1*sumx1))
print(correl.max(), correl.min(), correl.mean() )
tstat = correl*sqrt(days) / (1. - correl*correl)
print(tstat.max(), tstat.min(), tstat.mean() )
#----------------------------------------------------------------------
#Plot the trend stats

#colors = matplotlib.colormaps.get_cmap('terrain')
colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')

#bins = np.linspace(-1,2,30)
bins = [ -0.75, -0.5, -0.25, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.5, 0.75, 1.0, 1.5, 1.75 ]

slope *= (365.2422*10)
hist,binedges = np.histogram(slope, bins = bins)
print(hist, hist.sum() )
print(binedges)


# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

#cs = ax.pcolormesh(lons, lats, slope, vmin=vmin,vmax=vmax,cmap=colors, transform=ccrs.PlateCarree() )
cs = ax.pcolormesh(lons, lats, slope, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("decadal slope")
cb.set_label(cbarlabel, fontsize=12)

#plt.show()
plt.savefig("slope.png")
plt.close()

#----------------------------------------------------------------------
# Plot correlations
colors = matplotlib.colormaps.get_cmap('seismic')
bins = [ -0.75, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 0.75 ]
hist,binedges = np.histogram(correl, bins = bins)
print(hist, hist.sum() )
print(binedges)

bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, correl, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("correlation")
cb.set_label(cbarlabel, fontsize=12)

#plt.show()
plt.savefig("correl.png")
plt.close()

#----------------------------------------------------------------------
# Plot tstat
colors = matplotlib.colormaps.get_cmap('seismic')
bins = [ -100., -25., -5., -2., 0, 2, 5, 25, 100 ]
hist,binedges = np.histogram(tstat, bins = bins)
print(hist, hist.sum() )
print(binedges)

bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')
ax.gridlines()

cs = ax.pcolormesh(lons, lats, tstat, norm = norm, cmap = colors, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("tstatistic")
cb.set_label(cbarlabel, fontsize=12)

#plt.show()
plt.savefig("tstatistic.png")
plt.close()

 
exit(0)
#----------------------------------------------------------------------




#----------------------------------------------------------------------
#avg: 
avg = ma.masked_array(sumx1/days, mask) 

#RG: find and plot histogram -- area-weighted
#hist,binedges = np.histogram(avg, bins = 80, range = [-2,38] )
bins = [-2., -1.5, -1., -0.5, 0., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.,
        13., 14., 15., 16., 17., 18., 19., 20., 
        20.5, 21., 21.5, 22., 22.5, 23., 23.5, 24., 24.5, 25., 25.5, 
        26., 26.5, 27., 27.5, 28., 28.5, 29., 29.5, 30., 30.5 ]
hist,binedges = np.histogram(avg, bins = bins)

#debug: 
print(hist, hist.sum() )
#debug: 
print(binedges, flush=True)

print("plotting average histogram", flush=True)

fig,ax = plt.subplots()
plt.grid(visible=True)
ax.plot(binedges[0:-1], hist)
#plt.show()
plt.savefig("avg_hist.png")
plt.close()

#RG: plot geographic map in ~deciles ----------------


#----------------------------------------------------------------------
#var = 
var = ma.masked_array(zeros, mask)
var = (sumx2 - sumx1*sumx1/days)/days
tmask = ma.masked_array(var <= 0)
indices = tmask.nonzero()
print("var nonpositive ",len(indices[0]) )
del tmask

bins = [0, 0.25, 0.5, 1, 2, 3, 4, 6, 9, 16, 25, 100]
hist,binedges = np.histogram(var, bins = bins)
print("variance\n",hist)
print(binedges)
print(hist.sum() )
print(len(binedges), len(hist)) 

print("plotting histogram", flush=True)
#plot histogram
fig,ax = plt.subplots()
plt.grid(visible=True)
ax.plot(binedges[0:-1], hist)
#plt.show()
plt.savefig("var.png")
plt.close()

#RG:plot geographic map with these bounds
#colors = matplotlib.colormaps.get_cmap('terrain')
#colors = matplotlib.colormaps.get_cmap('seismic')
#colors = matplotlib.colormaps.get_cmap('bwr')

# For nonuniform binning:
bounds = np.array(bins)
norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.coastlines(resolution='10m')

cs = ax.pcolormesh(lons, lats, var, norm = norm, transform=ccrs.PlateCarree() )
cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
cbarlabel = '%s' % ("var temperature")
cb.set_label(cbarlabel, fontsize=12)

#plt.show()
plt.savefig("var_geo.png")
plt.close()

#----------------------------------------------------------------------
#skew = 

#----------------------------------------------------------------------
#kurtosis = 


