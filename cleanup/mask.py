import sys
from math import *
import datetime
import copy

import numpy as np
import numpy.ma as ma
import netCDF4

#-------------------------------------------------
# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"
#file name format: "oisst-avhrr-v02r01.YYYYMMDD.nc"

dt = datetime.timedelta(1)

# Defining the quarter degree grid
nx = 1440
ny = 720

start = datetime.datetime(1981,9,1)

end  = datetime.datetime(2024,12,8)
#end  = datetime.datetime(1981,9,30)
#end  = datetime.datetime(1986,8,31)

allcount = np.zeros((ny,nx))

#-------------------------------------------------
tag = start
count = 0
while (tag <= end ):
    if (count % 30 == 0):
      print(tag, flush=True, file = sys.stderr )

# Get the day's data:
    fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
    tmpnc = netCDF4.Dataset(fbase + fname)
    sst   = tmpnc.variables['sst'][0,0,:,:]
    tmpnc.close()
    tmp   = copy.deepcopy(sst)

    allcount[tmp.mask] += 1
    # allcount[~tmp.mask] += 1 -- for the unmasked points

    #indices = tmp.nonmask()
    #for k in range(0,len(indices[0])):
    #  j = indices[0][k]
    #  i = indices[1][k]
    #  allcount[j,i] += 1

# Advance
    count += 1   # number of days' data
    tag   += dt

#------------------------------------------------
print(allcount.max(), allcount.min(), allcount.mean() )
print(np.count_nonzero(allcount))
#------------------ End of first pass --------------------------
cmax = allcount.max()
for j in range(0, ny):
  for i in range(0, nx):
    if (allcount[j,i] < cmax and allcount[j,i] > 0):
      print(i,j,allcount[j,i], file = sys.stdout )

