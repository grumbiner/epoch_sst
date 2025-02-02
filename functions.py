import copy
from math import *
import numpy as np
import numpy.ma as ma

import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

# collection bin for miscellaneous functions

#=================================================================
def show(bins, lons, lats, x, title, fbase, cmap = matplotlib.colormaps.get_cmap('bwr'), proj = ccrs.PlateCarree() ):
  bounds = np.array(bins)
  norm = matplotlib.colors.BoundaryNorm(boundaries = bounds, ncolors = 256)

  fig = plt.figure(figsize=(12, 9))
  ax  = fig.add_subplot(1, 1, 1, projection = proj)

  #WNA: ax.set_extent((-95, -15, 0, 75),crs=proj)
  #Nino 3.4: ax.set_extent((-170, -120, -5, 5), crs=proj)

  ax.coastlines(resolution='10m')
  ax.gridlines()

  cs = ax.pcolormesh(lons, lats, x, norm = norm, cmap = cmap, transform=ccrs.PlateCarree() )
  cb = plt.colorbar(cs, extend='both', orientation='horizontal', shrink=0.5, pad=.04)
  cbarlabel = '%s' % title
  cb.set_label(cbarlabel, fontsize=12)

  plt.savefig(fbase+".png")
  plt.close()

  hist, binedges =  np.histogram(x, bins = bins)
  print(title)
  print(hist, hist.sum() )
  print(binedges,"\n\n")
  #debug: print(binedges)


#------------------------------------------------
#=================================================================
def applymask(mask, grid, indices):
  for k in range(0, len(indices[0])):
    i = indices[1][k]
    j = indices[0][k]
    grid[j,i] = 0.   

#=================================================================
def find_bins(x, nbin):
  vmin = floor(x.min())
  vmax = ceil(x.max())
  return np.linspace(vmin, vmax, nbin)

#=================================================================
def climo(intercept, slope, ampl, phase, freq, epoch, tag):
  delta = (tag - epoch).days
  sst = copy.deepcopy(intercept)
  sst += slope*delta
  for j in range(0,3):
    sst += ampl[j]*np.cos(phase[j] + freq[j]*delta)

  return sst

#----------------------------------------------------------------------

def old_climo(epoch, tag):
  fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#RG: This is hard wiring somewhat to epoch 1 sep 1981
  if (tag.month == 2 and tag.day == 29):
    ref = datetime.datetime(epoch.year, tag.month, 28)
  else:
    ref = datetime.datetime(epoch.year, tag.month, tag.day)

  if (ref < epoch):
    ref = datetime.datetime(epoch.year+1, ref.month, ref.day)
  fname = "traditional_" + ref.strftime("%Y%m%d") + ".nc"

  tmpnc = nc.Dataset(fbase + fname)
  sst = tmpnc.variables['mean'][:,:]
  tmpnc.close()

  return sst



