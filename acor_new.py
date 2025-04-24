from math import *
import os
import datetime
import copy

import numpy as np
import numpy.ma as ma
import netCDF4

from functions import *
import ncoutput

#------------------------------------------------

#-------------------------------------------------
#  Compute a traditional style climatology, day by day for 30 years
#-------------------------------------------------
# location of data files
fbase = "/Volumes/Data/qdoi/v2.1.nc/"

# Defining the quarter degree grid
nx = 1440
ny = 720

# Start-finish, but will be iterating through next 30 years
start = datetime.datetime(1981,9,1)
#end = datetime.datetime(1982,8,31)
end = datetime.datetime(2011,8,31)
nt = (end - start).days + 1
print("days ",nt)
stride = 4

dt = datetime.timedelta(1)
series = np.zeros(nt)
#---------------------------------------------
# Now run through the data files and accumulate terms:

tag = start

count = 0

sst = np.zeros((ny,nx)) # temporary file for reading in data

while (tag <= end ):
  if (count % 90 == 0):
    print("tag =",tag, flush=True)

# Get the day's data:
  fname = "newres1_" + tag.strftime("%Y%m%d") + ".nc"
#  fname = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"
#  fname = "res_traditional_" + tag.strftime("%Y%m%d") + ".nc"

  if not (tag.month == 2 and tag.day == 29):
    tmpnc = netCDF4.Dataset(fbase + fname)
    sst = tmpnc.variables['newres1'][:,:]
    #sst = tmpnc.variables['sst'][0,0,:,:]
    #sst = tmpnc.variables['tradres'][:,:]
    tmpnc.close()
  # sst persists for leap days

  # 45.125 N, 315.125E (44.875 W)
  series[count] = sst[int(ny*3/4), 1260 ]

  count += 1   # number of days' data
  tag   += dt
  
#-------------------------------------------------
print('done')

series -= np.mean(series)
for i in range(0, len(series) ):
  print(i, series[i])

import scipy
auto = scipy.signal.correlate(series, series, mode='full')
auto /= np.max(auto)
for i in range(0,len(auto)):
  print(i-10956, auto[i])

y = scipy.fft.fft(series)
yf = scipy.fft.fftfreq(len(series), 1)
y = np.abs(y)*2/len(series)
for i in range(0, len(y)):
  print(i, y[i], yf[i])

