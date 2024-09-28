

import netCDF4 as nc

# given the mean values of sumx1, sumx2, sumx3, sumx4, sumxt, sumt2, sumt, compute
#    mean, sd, skew, kurtosis, and trend

dset = nc.Dataset("first_pass.nc", "r")
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
mean = dset.variables['sumx1'][:,:]
mean2 = dset.variables['sumx2'][:,:]
mean3 = dset.variables['sumx3'][:,:]
mean4 = dset.variables['sumx4'][:,:]
meanxt = dset.variables['sumxt'][:,:]
meant2 = dset.variables['sumt2'][:,:]
sumt = 182.0 #from 1 yr
days = 365

print("mean ",mean.max(), mean.min() )
print("mean2 ",mean2.max(), mean2.min() )
print("mean3 ",mean3.max(), mean3.min() )
print("mean4 ",mean4.max(), mean4.min() )
print("meanxt ",meanxt.max(), meanxt.min() )
print("meant2 ",meant2.max(), meant2.min() )

#--------------------------------
#var = 
#skew = 
#kurtosis = 


