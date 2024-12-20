"""
@author: Jasmine
"""

import numpy as np
import math
import sys
import csv
from scipy.stats.stats import pearsonr
from scipy import interpolate
from sklearn.metrics import pairwise_distances_argmin
import os
import xarray as xr
import pandas as pd
import tarfile
import netCDF4 as nc
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
from datetime import datetime as dt
from datetime import timedelta as td
import matplotlib.dates as mdates
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, ScalarFormatter,LogFormatter)

# data = Dataset('ows_papa.nc')

# PAPAtemp = netCDF4.Dataset('../ows_papa/t50n145w_hr.cdf', 'r')
# TempPAPAtemp = xr.open_dataset('../ows_papa/t50n145w_hr.cdf', decode_times = False)
TempPAPAtemp = xr.open_dataset('../ows_papa/t50n145w_10m.cdf', decode_times = False)
# TempPAPAtemp = nc.Dataset('../ows_papa/t50n145w_10m.cdf', decode_times = False)

# owsPAPA = Dataset('ows_papa.nc')
PAPAtemp = TempPAPAtemp.T_20; PAPAdepth = np.negative(TempPAPAtemp.depth);  # access a variable in the file 
# SPAPAtemp = STempPAPAtemp.T_20; SPAPAdepth = STempPAPAtemp.depth;  # access a variable in the file 
# print('depth', PAPAdepth)
# PAPAtemp = TempPAPAtemp.variables['T_20'];
# PAPAtemp = PAPAtemp.fillna(PAPAtemp.mean())

# # convert 3D to 2D
PAPAtemp = PAPAtemp.isel(lat=0,lon=0)
## # convert to pandas dataframe
#PAPAtemp1 = PAPAtemp.to_dataframe()
## # change the missing values to the mean value
#PAPAtemp2 = PAPAtemp1.fillna(PAPAtemp1.mean()).reset_index()
## # convert back to xarray
#PAPAtemp3 = PAPAtemp2.to_xarray()
#print(PAPAtemp3)
# print('temperature', PAPAtemp)

# mod_PAPAtemp = []

# for tme in PAPAtemp.values:
# 	for dpth in range(len(tme)):
# 		print(tme[dpth])
		# if tme[dpth] > 100:
		# 	previous = PAPAtemp[dpth-1]
		# 	mod_PAPAtemp.append(previous)
		# else:
		# 	mod_PAPAtemp.append(dpth)
# mod_PAPAtemp = np.asarray(mod_PAPAtemp)
# time array
PAPAtime = TempPAPAtemp.time; #2012-09-04 00:00:00 8856
PAPAtime2 = 248 + (PAPAtime[:]/1440); # 2012-09-04 00:00:00 -- 2013-09-07 00:00:00


# PAPAtemp2 = []
# for row in PAPAtemp.values:
# 	#d = row.length
# 	_row = row

# 	for i in range(len(row)):
# 		if row[i] > 100:
# 			if i > 0:
# 				_row[i]=row[i-1]
# 			else:
# 				_row[i]=row[i]
# 	# print(_row)
# 	# print(row)
# 	# break;
# 	PAPAtemp2.append(_row)
#         #plt.plot(row[0])

# PAPAtemp2 = xr.DataArray(PAPAtemp2)
# PAPAtemp2 = xr.DataArray(PAPAtemp2, dims=("time","depth","lat","lon"), coords={"time": np.linspace(0.0,531350.0,num=10), "depth": [1.0,5.0,10.0,14.0,20.0,25.0,30.0,37.0,45.0,60.0,80.0,100.0,120.0,150.0,175.0,200.0,300.0], "lat":[50.1],"lon":[215.1]})
# print(PAPAtemp2)
# print(PAPAtemp)
# print(PAPAtime)
# print(PAPAdepth)

gNaN = PAPAtemp.where(PAPAtemp.time > 100)

print(len(gNaN))
exit()
print('where >100', np.where(PAPAtemp.values == 1e+35))

print('shape of temp before: ',PAPAtemp.shape)
PAPAtemp = PAPAtemp.where(PAPAtemp.time==1e+35, drop=True)
print('shape of temp after: ',PAPAtemp.shape)

#print('where >100', np.where(PAPAtemp.values>100))

exit()

#  mixed layer depth
mld_tempPAPA = PAPAtemp.sel(depth=10.0, method='nearest') - 0.2; 
z_indexes = abs(PAPAtemp-mld_tempPAPA).argmin('depth'); # print(z_indexes)
mldPAPA = PAPAtemp.depth[z_indexes]; 
# mldPAPA = mldPAPA.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# print('index for MLD :', MLDindx, MLDindx.shape)
# print('temp:', PAPAtemp[9640,:].values)#, mldPAPA.shape)
# print('PAPA time:', PAPAtime2[53000:].values, PAPAtime.shape)#, PAPAtime.shape)
# exit()
# print('depth:', PAPAdepth[53000:].values)#, mldPAPA.shape)
# print('mldKPP:', mldPAPA[0:10].values, mldPAPA.shape)#, mldPAPA.shape)

# print('temp:', PAPAtemp[9640,0:10].values)#, mldPAPA.shape)
# print(np.where(PAPAtemp[9640,0:10]  == 1.00000004e+35))

# print(np.where(PAPAtemp >1000))


plt.plot(PAPAtemp.time, mldPAPA)
plt.ylabel('Mixed layer depth -m')
plt.xlabel('Time from 04-Sept-2012 until 30-Sep-2013')
plt.title('OWS Papa Mixed layer depth evolution')
plt.show()

# print('modified temperature', mod_PAPAtemp)

# MLD_PAPAindx = []; 

# mld_tempPAPA = PAPAtemp.sel(depth=10.0, method='nearest') - 0.2;
# mld_tempPAPA = (PAPAtemp-mld_tempPAPA).isel(lat=0,lon=0).values
# # print('the difference',(mld_tempKPP)[0:50,469:499])
# # print('the difference',(mld_tempKPP)[0:50,469:499].shape)

# for time in mld_tempPAPA:
#     for dpth in range(len(time)):
#         if time[dpth] > 0:
#             MLD_PAPAindx.append(dpth)
#             break
# MLD_PAPAindx = np.asarray(MLD_PAPAindx); mldPAPA = TempPAPAtemp.depth[MLD_PAPAindx];
 ## reduce eerrors


# for x in enumerate(PAPAtemp.values):
#     print(PAPAtemp[x])
# print(PAPAtemp.shape)

# PAPAtemp_imp = []

# for i, ele in enumerate(PAPAtemp):
     
#     # replace if greater than K
#     if ele > 100 :
#         previous = PAPAtemp[i-1]
#         PAPAtemp_imp.append(previous)
#     else :
#         PAPAtemp_imp.append(ele)

# print('improved dataset of OWS Papa temp :', PAPAtemp[9640,0:10].values)
