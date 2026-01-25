'''
#  Compute autocorrelations of residuals from newstyle climatology
  # 45.125 N, 315.125E (44.875 W)
  series[count] = sst[int(ny*3/4), 1260 ]
'''
import datetime

import numpy as np
import netCDF4
import scipy

from functions import *

#-------------------------------------------------
# location of data files
fbase = "/Volumes/Data2/qdoi/v2.1.nc/"

# Defining the quarter degree grid
nx = 1440
ny = 720
  # 45.125 N, 315.125E (44.875 W)
nj = int(ny*3/4)
ni = int(1260)
  # equator, 315.125E
nj = int(ny/2)
print("nj ni: ",nj,ni)

# Start-finish, but will be iterating through next 30 years
epoch = datetime.datetime(1981,9,1)
#end = datetime.datetime(1982,8,31)
end = datetime.datetime(2011,8,31)
nt = (end - epoch).days + 1
#debug: print("days ",nt, flush=True)

dt = datetime.timedelta(1)
#---------------------------------------------
# Now run through the data files and accumulate terms:

series1 = np.zeros(nt)
series2 = np.zeros(nt)
series3 = np.zeros(nt)
series4 = np.zeros(nt)
tag = epoch
count = 0

while (tag <= end ):
  if (count % 90 == 0):
    print("tag =",tag, flush=True)

# Get the day's data:
  fname1 = "ninores1_" + tag.strftime("%Y%m%d") + ".nc"
  fname2 = "newres1_" + tag.strftime("%Y%m%d") + ".nc"
  fname3 = "res_traditional_" + tag.strftime("%Y%m%d") + ".nc"
  fname4 = "oisst-avhrr-v02r01." + tag.strftime("%Y%m%d") + ".nc"

  # sst persists for leap days
  if not (tag.month == 2 and tag.day == 29):
    tmpnc1 = netCDF4.Dataset(fbase + fname1)
    tmpnc2 = netCDF4.Dataset(fbase + fname2)
    tmpnc3 = netCDF4.Dataset(fbase + fname3)
    tmpnc4 = netCDF4.Dataset(fbase + fname4)

    sst1 = tmpnc1.variables['ninores1'][:,:]
    sst2 = tmpnc2.variables['newres1'][:,:]
    sst3 = tmpnc3.variables['fanomaly'][:,:]
    sst4 = tmpnc4.variables['sst'][0,0,:,:]
    tmpnc1.close()
    tmpnc2.close()
    tmpnc3.close()
    tmpnc4.close()

  # 45.125 N, 315.125E (44.875 W)
  series1[count] = sst1[nj, ni]
  series2[count] = sst2[nj, ni]
  series3[count] = sst3[nj, ni]
  series4[count] = sst4[nj, ni]

  count += 1   # number of days' data
  tag   += dt

#-------------------------------------------------
def show_series(fseries, ftag):

  # Print out the anomaly series
  #series -= np.mean(series)
  #for i in range(0, len(series) ):
  #  print(i, series[i])

  fname = ftag+".acor"
  fout = open(fname,"w",encoding = "utf-8")
  auto = scipy.signal.correlate(fseries, fseries, mode='full')
  print("max autocovariance ",np.max(auto) )
  auto /= np.max(auto)
  for i in range(0,len(auto)):
    print(i-10956, auto[i], file = fout)
  fout.close()

  fname = ftag+".fft"
  fout = open(fname,"w",encoding = "utf-8")
  y = scipy.fft.fft(fseries)
  yf = scipy.fft.fftfreq(len(fseries), 1)
  y = np.abs(y)*2/len(fseries)
  for i in range(0, len(y)):
    print(i, y[i], yf[i], file = fout)
  fout.close()
#-------------------------------------------------
show_series(series1, "nino34")
show_series(series2, "new")
show_series(series3, "old")
show_series(series4, "orig")
