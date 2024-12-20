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
import netCDF4
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
import matplotlib.dates as mdates
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, ScalarFormatter,LogFormatter)

# plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

### ----------------------------- GLIDER ----------------------------- ##
# file1 = 'C:/home/users/mc837749/Documents/gotm-4.0.0/simulations/annualOSMOSISforcing_2502/glider_timeseries.nc'
# gliderdata= Dataset('../glider_timeseries.nc'); 
gliderdata= xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc'); 
glider = xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; GldrDens = glider.pot_den; 
Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen
# print('original Gldr temp', GldrTemp.shape)
# GldrTempNaN = GldrTemp.isnull(); print('number of missing values in GldrTemp', GldrTempNaN)

# GldrTemp = GldrTemp.to_dataframe().dropna(how='all'); print('Gldr temp with dropped NaN values?', GldrTemp, GldrTemp.shape)
# GldrTemp = GldrTemp.to_dataframe();

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days
glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})

# # # check variables in file
# print(gliderdata.variables.items())
# print(GldrDens.values)

### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
gotmkeps = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_keps.nc", decode_times=False)

# gotmkppnum2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_num15.nc",decode_times=False)
# gotmkppnuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nuh15.nc",decode_times=False)
## LaT addition ##
gotmkppLaRib = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat.nc",decode_times=False)
## h/L correction ##
# gotmkppLaRibhLPveP100 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP100.nc",decode_times=False)
# gotmkppLaRibhLPveP1500 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP1500.nc",decode_times=False)
# gotmkppLaRibhLPveP10000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP10000.nc",decode_times=False)

# # check variables in files
# print(gotmkpp.variables.items())

#######################
## collect variables ##
#######################

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMsal = gotmkpp.salt; GTMo2 = gotmkpp.o2_obs; GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMxflx = gotmkpp.gamu; GTMyflx = gotmkpp.gamv; time = gotmkpp.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time
GTMtempLaRib = gotmkppLaRib.temp; GTMdepthLaRib = gotmkppLaRib.z; GTMnumLaRib = gotmkppLaRib.num; GTMnuhLaRib = gotmkppLaRib.nuh; GTMuLaRib = gotmkppLaRib.u; GTMvLaRib = gotmkppLaRib.v; timeLaRib = gotmkppLaRib.time
# GTMtempnum2 = gotmkppnum2.temp; GTMdepthnum2 = gotmkppnum2.z; GTMnumnum2 = gotmkppnum2.num; GTMnuhnum2 = gotmkppnum2.nuh; GTMunum2 = gotmkppnum2.u; GTMvnum2 = gotmkppnum2.v; timenum2 = gotmkppnum2.time
# GTMtempnuh2 = gotmkppnuh2.temp; GTMdepthnuh2 = gotmkppnuh2.z; GTMnumnuh2 = gotmkppnuh2.num; GTMnuhnuh2 = gotmkppnuh2.nuh; GTMunuh2 = gotmkppnuh2.u; GTMvnuh2 = gotmkppnuh2.v; timenuh2 = gotmkppnuh2.time

# GTMtempLaRibhLPveP100 = gotmkppLaRibhLPveP100.temp; GTMdepthLaRibhLPveP100 = gotmkppLaRibhLPveP100.z; GTMnumLaRibhLPveP100 = gotmkppLaRibhLPveP100.num; GTMnuhLaRibhLPveP100 = gotmkppLaRibhLPveP100.nuh; GTMuLaRibhLPveP100 = gotmkppLaRibhLPveP100.u; GTMvLaRibhLPveP100 = gotmkppLaRibhLPveP100.v; timeLaRibhLPveP100 = gotmkppLaRibhLPveP100.time
# GTMtempLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.temp; GTMdepthLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.z; GTMnumLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.num; GTMnuhLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.nuh; GTMuLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u; GTMvLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.v; timeLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.time
# GTMtempLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.temp; GTMdepthLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.z; GTMnumLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.num; GTMnuhLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.nuh; GTMuLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u; GTMvLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.v; timeLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.time


##################
## convert time ##
##################

time2 = 267.75 + (time[:]/86400); 
ketime2 = 267.75 + (ketime[:]/86400);
time2LaRib = 267.75 + (timeLaRib[:]/86400); 
# time2num2 = 267.75 + (timenum2[:]/86400); 
# time2nuh2 = 267.75 + (timenuh2[:]/86400); 

# time2LaRibhLPveP100 = 267.75 + (timeLaRibhLPveP100[:]/86400); 
# time2LaRibhLPveP1500 = 267.75 + (timeLaRibhLPveP1500[:]/86400); 
# time2LaRibhLPveP10000 = 267.75 + (timeLaRibhLPveP10000[:]/86400); 

###################
# time allocation #
###################

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})
# gotmkppL = gotmkppL.assign(time2L=("time", time2L)); gotmkppL = gotmkppL.swap_dims({"time" : "time2L"})
# gotmkppLa = gotmkppLa.assign(time2La=("time", time2La)); gotmkppLa = gotmkppLa.swap_dims({"time" : "time2La"})
gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})
# gotmkppLaRibL = gotmkppLaRibL.assign(time2LaRibL=("time", time2LaRibL)); gotmkppLaRibL = gotmkppLaRibL.swap_dims({"time" : "time2LaRibL"})
# gotmkppnum2 = gotmkppnum2.assign(time2num2=("time", time2num2)); gotmkppnum2 = gotmkppnum2.swap_dims({"time" : "time2num2"})
# gotmkppnuh2 = gotmkppnuh2.assign(time2nuh2=("time", time2nuh2)); gotmkppnuh2 = gotmkppnuh2.swap_dims({"time" : "time2nuh2"})
# # gotmkppRibnum2 = gotmkppRibnum2.assign(time2Ribnum2=("time", time2Ribnum2)); gotmkppRibnum2 = gotmkppRibnum2.swap_dims({"time" : "time2Ribnum2"})
# gotmkppRib05num2 = gotmkppRib05num2.assign(time2Rib05num2=("time", time2Rib05num2)); gotmkppRib05num2 = gotmkppRib05num2.swap_dims({"time" : "time2Rib05num2"})
# # gotmkppRibnuh2 = gotmkppRibnuh2.assign(time2Ribnuh2=("time", time2Ribnuh2)); gotmkppRibnuh2 = gotmkppRibnuh2.swap_dims({"time" : "time2Ribnuh2"})
# gotmkppRib2nuh2 = gotmkppRib2nuh2.assign(time2Rib2nuh2=("time", time2Rib2nuh2)); gotmkppRib2nuh2 = gotmkppRib2nuh2.swap_dims({"time" : "time2Rib2nuh2"})

# gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.assign(time2LaRibhLPveP100=("time", time2LaRibhLPveP100)); gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.swap_dims({"time" : "time2LaRibhLPveP100"})
# # gotmkppLaRibhLPveP200 = gotmkppLaRibhLPveP200.assign(time2LaRibhLPveP200=("time", time2LaRibhLPveP200)); gotmkppLaRibhLPveP200 = gotmkppLaRibhLPveP200.swap_dims({"time" : "time2LaRibhLPveP200"})
# # gotmkppLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.assign(time2LaRibhLPveP1000=("time", time2LaRibhLPveP1000)); gotmkppLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.swap_dims({"time" : "time2LaRibhLPveP1000"})
# # gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.assign(time2LaRibhLPveP2000=("time", time2LaRibhLPveP2000)); gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.swap_dims({"time" : "time2LaRibhLPveP2000"})
# # gotmkppLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.assign(time2LaRibhLPveP1400=("time", time2LaRibhLPveP1400)); gotmkppLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.swap_dims({"time" : "time2LaRibhLPveP1400"})
# gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.assign(time2LaRibhLPveP1500=("time", time2LaRibhLPveP1500)); gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.swap_dims({"time" : "time2LaRibhLPveP1500"})
# # gotmkppLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.assign(time2LaRibhLPveP1600=("time", time2LaRibhLPveP1600)); gotmkppLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.swap_dims({"time" : "time2LaRibhLPveP1600"})
# # gotmkppLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.assign(time2LaRibhLPveP1300=("time", time2LaRibhLPveP1300)); gotmkppLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.swap_dims({"time" : "time2LaRibhLPveP1300"})
# # gotmkppLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.assign(time2LaRibhLPveP1750=("time", time2LaRibhLPveP1750)); gotmkppLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.swap_dims({"time" : "time2LaRibhLPveP1750"})
# # gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.assign(time2LaRibhLPveP2000=("time", time2LaRibhLPveP2000)); gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.swap_dims({"time" : "time2LaRibhLPveP2000"})
# # gotmkppLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.assign(time2LaRibhLPveP4000=("time", time2LaRibhLPveP4000)); gotmkppLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.swap_dims({"time" : "time2LaRibhLPveP4000"})
# # gotmkppLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.assign(time2LaRibhLPveP1350=("time", time2LaRibhLPveP1350)); gotmkppLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.swap_dims({"time" : "time2LaRibhLPveP1350"})
# # gotmkppLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.assign(time2LaRibhLPveP5000=("time", time2LaRibhLPveP5000)); gotmkppLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.swap_dims({"time" : "time2LaRibhLPveP5000"})
# gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.assign(time2LaRibhLPveP10000=("time", time2LaRibhLPveP10000)); gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.swap_dims({"time" : "time2LaRibhLPveP10000"})
# # gotmkppLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.assign(time2LaRibhLPveP10000shear=("time", time2LaRibhLPveP10000shear)); gotmkppLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.swap_dims({"time" : "time2LaRibhLPveP10000shear"})

# print(gotmkppLaRibhLPveP2000.variables.items()

## September 2012
kppSep12 = np.where((time2 > 244) & (time2 <= 274)); kppSep12_time = time2[kppSep12]; 
kppSep12_u = GTMu.isel(lat=0,lon=0)[0:899,:]; kppSep12_u = kppSep12_u.mean(axis=0)
kppSep12_v = GTMv.isel(lat=0,lon=0)[0:899,:]; kppSep12_v = kppSep12_v.mean(axis=0)

kepsSep12 = np.where((ketime2 > 244) & (ketime2 <= 274)); kepsSep12_time = ketime2[kepsSep12]; 
kepsSep12_u = GTMkeu.isel(lat=0,lon=0)[0:899,:]; kepsSep12_u = kepsSep12_u.mean(axis=0)
kepsSep12_v = GTMkev.isel(lat=0,lon=0)[0:899,:]; kepsSep12_v = kepsSep12_v.mean(axis=0)

## October 2012
kppOct12 = np.where((time2 > 274) & (time2 <= 305)); kppOct12_time = time2[kppOct12]; 
kppOct12_u = GTMu.isel(lat=0,lon=0)[899:5363,:]; kppOct12_u = kppOct12_u.mean(axis=0)
kppOct12_v = GTMv.isel(lat=0,lon=0)[899:5363,:]; kppOct12_v = kppOct12_v.mean(axis=0)

kepsOct12 = np.where((ketime2 > 274) & (ketime2 <= 305)); kepsOct12_time = ketime2[kepsOct12]; 
kepsOct12_u = GTMkeu.isel(lat=0,lon=0)[899:5363,:]; kepsOct12_u = kepsOct12_u.mean(axis=0)
kepsOct12_v = GTMkev.isel(lat=0,lon=0)[899:5363,:]; kepsOct12_v = kepsOct12_v.mean(axis=0)

## November 2012
kppNov12 = np.where((time2 > 305) & (time2 <= 335)); kppNov12_time = time2[kppNov12]; 
kppNov12_u = GTMu.isel(lat=0,lon=0)[5363:9683,:]; kppNov12_u = kppNov12_u.mean(axis=0)
kppNov12_v = GTMv.isel(lat=0,lon=0)[5363:9683,:]; kppNov12_v = kppNov12_v.mean(axis=0)

kepsNov12 = np.where((ketime2 > 305) & (ketime2 <= 335)); kepsNov12_time = ketime2[kepsNov12]; 
kepsNov12_u = GTMkeu.isel(lat=0,lon=0)[5363:9683,:]; kepsNov12_u = kepsNov12_u.mean(axis=0)
kepsNov12_v = GTMkev.isel(lat=0,lon=0)[5363:9683,:]; kepsNov12_v = kepsNov12_v.mean(axis=0)

## December 2012
kppDec12 = np.where((time2 > 335) & (time2 <= 366)); kppDec12_time = time2[kppDec12]; 
kppDec12_u = GTMu.isel(lat=0,lon=0)[9683:14147,:]; kppDec12_u = kppDec12_u.mean(axis=0)
kppDec12_v = GTMv.isel(lat=0,lon=0)[9683:14147,:]; kppDec12_v = kppDec12_v.mean(axis=0)

kepsDec12 = np.where((ketime2 > 335) & (ketime2 <= 366)); kepsDec12_time = ketime2[kepsDec12]; 
kepsDec12_u = GTMkeu.isel(lat=0,lon=0)[9683:14147,:]; kepsDec12_u = kepsDec12_u.mean(axis=0)
kepsDec12_v = GTMkev.isel(lat=0,lon=0)[9683:14147,:]; kepsDec12_v = kepsDec12_v.mean(axis=0)

## January 2013
kppJan13 = np.where((time2 > 366) & (time2 <= 397)); kppJan13_time = time2[kppJan13]; 
kppJan13_u = GTMu.isel(lat=0,lon=0)[14147:18611,:]; kppJan13_u = kppJan13_u.mean(axis=0)
kppJan13_v = GTMv.isel(lat=0,lon=0)[14147:18611,:]; kppJan13_v = kppJan13_v.mean(axis=0)

kepsJan13 = np.where((ketime2 > 366) & (ketime2 <= 397)); kepsJan13_time = ketime2[kepsJan13]; 
kepsJan13_u = GTMkeu.isel(lat=0,lon=0)[14147:18611,:]; kepsJan13_u = kepsJan13_u.mean(axis=0)
kepsJan13_v = GTMkev.isel(lat=0,lon=0)[14147:18611,:]; kepsJan13_v = kepsJan13_v.mean(axis=0)

## Februay 2013
kppFeb13 = np.where((time2 > 397) & (time2 <= 456)); kppFeb13_time = time2[kppFeb13]; 
kppFeb13_u = GTMu.isel(lat=0,lon=0)[18611:27107,:]; kppFeb13_u = kppFeb13_u.mean(axis=0)
kppFeb13_v = GTMv.isel(lat=0,lon=0)[18611:27107,:]; kppFeb13_v = kppFeb13_v.mean(axis=0)

kepsFeb13 = np.where((ketime2 > 397) & (ketime2 <= 456)); kepsFeb13_time = ketime2[kepsFeb13]; 
kepsFeb13_u = GTMkeu.isel(lat=0,lon=0)[18611:27107,:]; kepsFeb13_u = kepsFeb13_u.mean(axis=0)
kepsFeb13_v = GTMkev.isel(lat=0,lon=0)[18611:27107,:]; kepsFeb13_v = kepsFeb13_v.mean(axis=0)

## March 2013
kppMar13 = np.where((time2 > 456) & (time2 <= 487)); kppMar13_time = time2[kppMar13]; 
kppMar13_u = GTMu.isel(lat=0,lon=0)[27107:31571,:]; kppMar13_u = kppMar13_u.mean(axis=0)
kppMar13_v = GTMv.isel(lat=0,lon=0)[27107:31571,:]; kppMar13_v = kppMar13_v.mean(axis=0)

kepsMar13 = np.where((ketime2 > 456) & (ketime2 <= 487)); kepsMar13_time = ketime2[kepsMar13]; 
kepsMar13_u = GTMkeu.isel(lat=0,lon=0)[27107:31571,:]; kepsMar13_u = kepsMar13_u.mean(axis=0)
kepsMar13_v = GTMkev.isel(lat=0,lon=0)[27107:31571,:]; kepsMar13_v = kepsMar13_v.mean(axis=0)

## April 2013
kppApr13 = np.where((time2 > 487) & (time2 <= 517)); kppApr13_time = time2[kppApr13]; 
kppApr13_u = GTMu.isel(lat=0,lon=0)[31571:35891,:]; kppApr13_u = kppApr13_u.mean(axis=0)
kppApr13_v = GTMv.isel(lat=0,lon=0)[31571:35891,:]; kppApr13_v = kppApr13_v.mean(axis=0)

kepsApr13 = np.where((ketime2 > 487) & (ketime2 <= 517)); kepsApr13_time = ketime2[kepsApr13]; 
kepsApr13_u = GTMkeu.isel(lat=0,lon=0)[31571:35891,:]; kepsApr13_u = kepsApr13_u.mean(axis=0)
kepsApr13_v = GTMkev.isel(lat=0,lon=0)[31571:35891,:]; kepsApr13_v = kepsApr13_v.mean(axis=0)

## May 2013
kppMay13 = np.where((time2 > 517) & (time2 <= 548)); kppMay13_time = time2[kppMay13]; 
kppMay13_u = GTMu.isel(lat=0,lon=0)[35891:40355,:]; kppMay13_u = kppMay13_u.mean(axis=0)
kppMay13_v = GTMv.isel(lat=0,lon=0)[35891:40355,:]; kppMay13_v = kppMay13_v.mean(axis=0)

kepsMay13 = np.where((ketime2 > 517) & (ketime2 <= 548)); kepsMay13_time = ketime2[kepsMay13]; 
kepsMay13_u = GTMkeu.isel(lat=0,lon=0)[35891:40355,:]; kepsMay13_u = kepsMay13_u.mean(axis=0)
kepsMay13_v = GTMkev.isel(lat=0,lon=0)[35891:40355,:]; kepsMay13_v = kepsMay13_v.mean(axis=0)

## June 2013
kppJun13 = np.where((time2 > 548) & (time2 <= 578)); kppJun13_time = time2[kppJun13]; 
kppJun13_u = GTMu.isel(lat=0,lon=0)[40355:44675,:]; kppJun13_u = kppJun13_u.mean(axis=0)
kppJun13_v = GTMv.isel(lat=0,lon=0)[40355:44675,:]; kppJun13_v = kppJun13_v.mean(axis=0)

kepsJun13 = np.where((ketime2 > 548) & (ketime2 <= 578)); kepsJun13_time = ketime2[kepsJun13]; 
kepsJun13_u = GTMkeu.isel(lat=0,lon=0)[40355:44675,:]; kepsJun13_u = kepsJun13_u.mean(axis=0)
kepsJun13_v = GTMkev.isel(lat=0,lon=0)[40355:44675,:]; kepsJun13_v = kepsJun13_v.mean(axis=0)

## July 2013
kppJul13 = np.where((time2 > 578) & (time2 <= 609)); kppJul13_time = time2[kppJul13]; 
kppJul13_u = GTMu.isel(lat=0,lon=0)[44675:49139,:]; kppJul13_u = kppJul13_u.mean(axis=0)
kppJul13_v = GTMv.isel(lat=0,lon=0)[44675:49139,:]; kppJul13_v = kppJul13_v.mean(axis=0)

kepsJul13 = np.where((ketime2 > 578) & (ketime2 <= 609)); kepsJul13_time = ketime2[kepsJul13]; 
kepsJul13_u = GTMkeu.isel(lat=0,lon=0)[44675:49139,:]; kepsJul13_u = kepsJul13_u.mean(axis=0)
kepsJul13_v = GTMkev.isel(lat=0,lon=0)[44675:49139,:]; kepsJul13_v = kepsJul13_v.mean(axis=0)

## August 2013
kppAug13 = np.where((time2 > 609) & (time2 <= 640)); kppAug13_time = time2[kppAug13]; 
kppAug13_u = GTMu.isel(lat=0,lon=0)[49139:52595,:]; kppAug13_u = kppAug13_u.mean(axis=0)
kppAug13_v = GTMv.isel(lat=0,lon=0)[49139:52595,:]; kppAug13_v = kppAug13_v.mean(axis=0)

kepsAug13 = np.where((ketime2 > 609) & (ketime2 <= 640)); kepsAug13_time = ketime2[kepsAug13]; 
kepsAug13_u = GTMkeu.isel(lat=0,lon=0)[49139:52595,:]; kepsAug13_u = kepsAug13_u.mean(axis=0)
kepsAug13_v = GTMkev.isel(lat=0,lon=0)[49139:52595,:]; kepsAug13_v = kepsAug13_v.mean(axis=0)

# concatenate monthly data
GTMkpp_u = np.c_[kppSep12_u,kppOct12_u,kppNov12_u,kppDec12_u,kppJan13_u,kppFeb13_u,kppMar13_u,kppApr13_u,kppMay13_u,
    kppJun13_u,kppJul13_u,kppAug13_u]
GTMkeps_u = np.c_[kepsSep12_u,kepsOct12_u,kepsNov12_u,kepsDec12_u,kepsJan13_u,kepsFeb13_u,kepsMar13_u,kepsApr13_u,
    kepsMay13_u, kepsJun13_u,kepsJul13_u,kepsAug13_u]

GTMkpp_v = np.c_[kppSep12_v,kppOct12_v,kppNov12_v,kppDec12_v,kppJan13_v,kppFeb13_v,kppMar13_v,kppApr13_v,kppMay13_v,
    kppJun13_v,kppJul13_v,kppAug13_v]
GTMkeps_v = np.c_[kepsSep12_v,kepsOct12_v,kepsNov12_v,kepsDec12_v,kepsJan13_v,kepsFeb13_v,kepsMar13_v,kepsApr13_v,
    kepsMay13_v,kepsJun13_v,kepsJul13_v,kepsAug13_v]
# gliderOS_time = np.c_[Sep12_time,Oct12_time,Nov12_time,Dec12_time,
    # Jan13_time,Feb13_time,Mar13_time,Apr13_time,May13_time,Jun13_time,Jul13_time,Aug13_time]

GTMkpp_u = np.array(GTMkpp_u[0])
GTMkeps_u = np.array(GTMkeps_u[0])

GTMkpp_v = np.array(GTMkpp_v[0])
GTMkeps_v = np.array(GTMkeps_v[0])
# print('glider time', gliderOS_time)
# print('glider SST', GTM_u)
# exit()



#############################################################################################
### ----------------------------- Ocean Reanalysis System 5 ----------------------------- ###
#############################################################################################

# ORA5_temp = xr.open_dataset('votemper_combined.nc', decode_times=False);
# ORA5_salt = xr.open_dataset('vosaline_combined.nc', decode_times=False);
ORA5_wtrflx = xr.open_dataset('sowaflup_combined.nc', decode_times=False)
ORA5_zonvel = xr.open_dataset('vozocrtx_combined.nc', decode_times=False)
ORA5_merivel = xr.open_dataset('vomecrty_combined.nc', decode_times=False)

# # combine netcdf4 files with same predixes
# ORA5_zonvel = xr.open_mfdataset('vozocrtx_control_monthly_highres_3D_*.nc')
# ORA5_zonvel = ORA5_zonvel.to_netcdf('vozocrtx_combined.nc')

# ORA5_merivel = xr.open_mfdataset('vomecrty_control_monthly_highres_3D_*.nc')
# ORA5_merivel = ORA5_merivel.to_netcdf('vomecrty_combined.nc')

# ORA_t_Depth = ORA5_temp.deptht; ORA_t_Temp2 = ORA5_temp.votemper; ORA_t_time2 = ORA5_temp.time_counter
# ORA_t_lat = ORA5_temp.nav_lat; ORA_t_lon = ORA5_temp.nav_lon;

# print('ORA5 time', ORAtime2.values)
print('ORA5 dataset', ORA5_zonvel.variables.items())
print('ORA5 dataset', ORA5_merivel.variables.items())
# print('ORAtemp3D items', ORA5_temp.variables.items())
# print('ORAtemp3D latitude', ORAlat[:,0:10], ORAlat.shape)
# print('ORAtemp3D longitude', ORAlon[:,0:10].values, ORAlon.shape)
# print('latitude', np.min(ORA_t_lat), np.max(ORA_t_lat))
# print('longitude', np.min(ORA_t_lon), np.max(ORA_t_lon)


#------------------------------------------------#
##### extract data for selected co-ordinates #####
#------------------------------------------------#

#ZONAL
def extract(ORA5_zonvel, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_zonvel['nav_lon'].values.ravel(), ORA5_zonvel['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_zonvel['nav_lon'].shape))
    return ORA5_zonvel.isel(x=j0, y=i0).squeeze()
selectORAZonalO = extract(ORA5_zonvel, -16.25, 48.75);# print('extracts',selectORATemp)
selectORAZonal1 = extract(ORA5_zonvel, -16.25, 48.5);#North
selectORAZonal2 = extract(ORA5_zonvel, -16.0, 48.75);#East
selectORAZonal3 = extract(ORA5_zonvel, -16.0, 48.5);#West
# selectORAZonal4 = extract(ORA5_zonvel, -16.25, 46.75);#South


#MERIDONAL
def extract(ORA5_merivel, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_merivel['nav_lon'].values.ravel(), ORA5_merivel['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_merivel['nav_lon'].shape))
    return ORA5_merivel.isel(x=j0, y=i0).squeeze()
selectORAMeridO = extract(ORA5_merivel, -16.25, 48.75);# print('extracts',selectORASalt)
selectORAMerid1 = extract(ORA5_merivel, -16.25, 48.5);#North
selectORAMerid2 = extract(ORA5_merivel, -16.0, 48.75);#East
selectORAMerid3 = extract(ORA5_merivel, -16.0, 48.5);#West
# selectORAMerid4 = extract(ORA5_merivel, -16.25, 46.75);#South


#WATER FLUX
def extract(ORA5_wtrflx, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_wtrflx['nav_lon'].values.ravel(), ORA5_wtrflx['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_wtrflx['nav_lon'].shape))
    return ORA5_wtrflx.isel(x=j0, y=i0).squeeze()
selectORAwFluxO = extract(ORA5_wtrflx, -16.25, 48.75);# print('extracts',selectORASalt)
selectORAwFlux1 = extract(ORA5_wtrflx, -16.25, 50.75);#North
selectORAwFlux2 = extract(ORA5_wtrflx, -14.25, 48.75);#East
selectORAwFlux3 = extract(ORA5_wtrflx, -18.25, 48.75);#West
selectORAwFlux4 = extract(ORA5_wtrflx, -16.25, 46.75);#South

# ## extract depth information for specific time counters
# selectORAZonalO = selectORAZonalO.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAZonal1 = selectORAZonal1.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAZonal2 = selectORAZonal2.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAZonal3 = selectORAZonal3.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAZonal4 = selectORAZonal4.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)


# selectORAMeridO = selectORAMeridO.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAMerid1 = selectORAMerid1.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAMerid2 = selectORAMerid2.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAMerid3 = selectORAMerid3.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORAMerid4 = selectORAMerid4.sel(time_counter=0)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)

#--------------------------#
##### plotting figures #####
#--------------------------#
#ZONAL and MERIDONAL
month = np.array([datetime(2012,9,1),datetime(2012,10,1),datetime(2012,11,1),datetime(2012,12,1),
    datetime(2013,1,1),datetime(2013,2,1),datetime(2013,3,1),datetime(2013,4,1),datetime(2013,5,1),
    datetime(2013,6,1),datetime(2013,7,1),datetime(2013,8,1)])

# fig,ax = plt.subplots(1,1)
# # # plt.plot(selectORAwFluxO.time_counter,selectORAwFluxO.sowaflup, label = r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax.plot(month,selectORAZonalO.vozocrtx, label = r'ORA5 Zonal (48.75 $\degree$N, 16.25 $\degree$W)')
# # ax.plot(month,selectORAMeridO.vomecrty, label = r'GTM Meridonal (48.75 $\degree$N, 16.25 $\degree$W)')
# # ax.plot(month,selectORAwFlux1.sowaflup, linestyle = 'dashed', label = r'ORA5 (48.5 $\degree$N, 16.25 $\degree$W)')
# # ax.plot(month,selectORAwFlux2.sowaflup, linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.00 $\degree$W)')
# # ax.plot(month,selectORAwFlux3.sowaflup, linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.00 $\degree$W)')
# ax.legend(loc='lower right')
# ax.set_xlim(datetime(2012,9,1),datetime(2013, 8, 1))
# plt.ylabel(r'Velocity - $m s^{-1}$', fontsize=12); plt.xlabel(r'Time -days since 2012-09-16 00:00:00', fontsize=12                                                                                      )
# # # plt.ylim(-1000,0); plt.xlim(35.1,35.8)
# plt.grid()
# plt.show()


# #ZONAL
# plt.plot(selectORAZonalO.vozocrtx,np.negative(selectORAZonalO.depthu), label = r'ORA5 (48.7 $\degree$N, 16.2 $\degree$W) 2012-09')
# plt.plot(selectORAZonal1.vozocrtx,np.negative(selectORAZonal1.depthu), linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W) 2012-09')
# plt.plot(selectORAZonal2.vozocrtx,np.negative(selectORAZonal2.depthu), linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W) 2012-09')
# plt.plot(selectORAZonal3.vozocrtx,np.negative(selectORAZonal3.depthu), linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W) 2012-09')
# plt.plot(selectORAZonal4.vozocrtx,np.negative(selectORAZonal4.depthu), linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W) 2012-09')
# plt.plot(Sep12_u,GTMdepth, label = 'GOTM 2012-09')
# plt.legend(loc='lower right')
# plt.ylabel('Depth -m'); plt.xlabel(r'Zonal Velocity - $m s^{-1}$')
# plt.ylim(-1000,0);# plt.xlim(7,19)
# plt.grid()
# plt.show()


# fig,ax = plt.subplots(1,3,sharey=True,sharex=True)

# ax[0].plot(selectORAZonalO.vozocrtx.sel(time_counter=0),np.negative(selectORAZonalO.depthu), color = 'blue', label = r'ORA5 (48.75$\degree$N, 16.25$\degree$W)')
# # ax[0].plot(selectORAZonal1.vozocrtx.sel(time_counter=91),np.negative(selectORAZonal1.depthu), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# # ax[0].plot(selectORAZonal2.vozocrtx.sel(time_counter=91),np.negative(selectORAZonal2.depthu), linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.20 $\degree$W)')
# # ax[0].plot(selectORAZonal3.vozocrtx.sel(time_counter=91),np.negative(selectORAZonal3.depthu), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.20 $\degree$W)')
# ax[0].plot(kppSep12_u,GTMdepth, label = 'GOTM KPP scheme', color = 'red')
# ax[0].plot(kepsSep12_u,GTMdepth, label = r'GOTM k-$\epsilon$ scheme', color = 'green')
# ax[0].grid();
# ax[0].set_ylim(-350,0);# ax[0].set_xlim(7,19);
# ax[0].tick_params(axis='both', which='major', labelsize=12)

# ax[1].plot(selectORAZonalO.vozocrtx.sel(time_counter=30),np.negative(selectORAZonalO.depthu), color = 'blue', label = r'ORA5 (48.75$\degree$N, 16.25$\degree$W)')
# # ax[1].plot(selectORAZonal1.vozocrtx.sel(time_counter=122),np.negative(selectORAZonal1.depthu), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# # ax[1].plot(selectORAZonal2.vozocrtx.sel(time_counter=30),np.negative(selectORAZonal2.depthu), linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.20 $\degree$W)')
# # ax[1].plot(selectORAZonal3.vozocrtx.sel(time_counter=30),np.negative(selectORAZonal3.depthu), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.20 $\degree$W)')
# ax[1].plot(kppOct12_u,GTMdepth, label = 'GOTM KPP scheme', color = 'red')
# ax[1].plot(kepsOct12_u,GTMdepth, label = r'GOTM k-$\epsilon$ scheme', color = 'green')
# ax[1].grid(); 
# ax[1].set_ylim(-350,0);# ax[1].set_xlim(7,19);
# ax[1].tick_params(axis='both', which='major', labelsize=12)

# ax[2].plot(selectORAZonalO.vozocrtx.sel(time_counter=61),np.negative(selectORAZonalO.depthu), color = 'blue', label = r'ORA5 (48.75$\degree$N, 16.25$\degree$W)')
# # ax[2].plot(selectORAZonal1.vozocrtx.sel(time_counter=61),np.negative(selectORAZonal1.depthu), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# # ax[2].plot(selectORAZonal2.vozocrtx.sel(time_counter=61),np.negative(selectORAZonal2.depthu), linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.20 $\degree$W)')
# # ax[2].plot(selectORAZonal3.vozocrtx.sel(time_counter=61),np.negative(selectORAZonal3.depthu), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.20 $\degree$W)')
# ax[2].plot(kppNov12_u,GTMdepth, label = 'GOTM KPP scheme', color = 'red')
# ax[2].plot(kepsNov12_u,GTMdepth, label = r'GOTM k-$\epsilon$ scheme', color = 'green')
# ax[2].grid(); ax[2].legend(bbox_to_anchor=(1, 0.90),loc='center left', fontsize='small')
# ax[2].set_ylim(-350,0);# ax[2].set_xlim(7,19);
# ax[2].tick_params(axis='both', which='major', labelsize=12)


# fig.text(0.5, 0.01, r'Zonal Velocity - $m s^{-1}$', ha='center', va='center', fontsize=12)
# fig.text(0.008, 0.5, r'Ocean Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')

# ax[0].text(0.025, 0.05, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
# ax[0].text(0.9, 0.045, 'SEP', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[0].transAxes)

# ax[1].text(0.025, 0.05, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[1].text(0.9, 0.045, 'OCT', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[1].transAxes)

# ax[2].text(0.025, 0.05, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
# ax[2].text(0.9, 0.045, 'NOV', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[2].transAxes)

# plt.show()


# # #MERIDONAL
# # plt.plot(selectORAMeridO.vomecrty,np.negative(selectORAMeridO.depthv), label = r'ORA5 (48.7 $\degree$N, 16.2 $\degree$W) 2012-09')
# # plt.plot(selectORAMerid1.vomecrty,np.negative(selectORAMerid1.depthv), linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W) 2012-09')
# # plt.plot(selectORAMerid2.vomecrty,np.negative(selectORAMerid2.depthv), linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W) 2012-09')
# # plt.plot(selectORAMerid3.vomecrty,np.negative(selectORAMerid3.depthv), linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W) 2012-09')
# # plt.plot(selectORAMerid4.vomecrty,np.negative(selectORAMerid4.depthv), linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W) 2012-09')
# # plt.plot(Sep12_v,GTMdepth, label = 'GOTM 2012-09')
# # plt.legend(loc='lower right')
# # plt.ylabel('Depth -m'); plt.xlabel(r'Meridonal Velocity - $m s^{-1}$')
# # plt.ylim(-1000,0);# plt.xlim(35.1,35.8)
# # plt.grid()
# # plt.show()

# fig,ax = plt.subplots(1,3,sharey=True,sharex=True)

# ax[0].plot(selectORAMeridO.vomecrty.sel(time_counter=0),np.negative(selectORAMeridO.depthv), label = r'ORA5 (48.75$\degree$N, 16.25$\degree$W)')
# # ax[0].plot(selectORAMerid1.vomecrty.sel(time_counter=0),np.negative(selectORAMerid1.depthv), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# # ax[0].plot(selectORAMerid2.vomecrty.sel(time_counter=0),np.negative(selectORAMerid2.depthv), linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.00 $\degree$W)')
# # ax[0].plot(selectORAMerid3.vomecrty.sel(time_counter=0),np.negative(selectORAMerid3.depthv), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.00 $\degree$W)')
# ax[0].plot(kppSep12_v,GTMdepth, label = 'GOTM KPP scheme')
# ax[0].plot(kepsSep12_v,GTMdepth, label = r'GOTM k-$\epsilon$ scheme')
# ax[0].grid();
# ax[0].set_ylim(-350,0);# ax[0].set_xlim(7,19);
# ax[0].tick_params(axis='both', which='major', labelsize=12)

# ax[1].plot(selectORAMeridO.vomecrty.sel(time_counter=30),np.negative(selectORAMeridO.depthv), label = r'ORA5 (48.75$\degree$N, 16.25$\degree$W)')
# # ax[1].plot(selectORAMerid1.vomecrty.sel(time_counter=30),np.negative(selectORAMerid1.depthv), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# # ax[1].plot(selectORAMerid2.vomecrty.sel(time_counter=30),np.negative(selectORAMerid2.depthv), linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.20 $\degree$W)')
# # ax[1].plot(selectORAMerid3.vomecrty.sel(time_counter=30),np.negative(selectORAMerid3.depthv), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.20 $\degree$W)')
# ax[1].plot(kppOct12_v,GTMdepth, label = 'GOTM KPP scheme')
# ax[1].plot(kepsOct12_v,GTMdepth, label = r'GOTM k-$\epsilon$ scheme')
# ax[1].grid(); 
# ax[1].set_ylim(-350,0);# ax[1].set_xlim(7,19);
# ax[1].tick_params(axis='both', which='major', labelsize=12)

# ax[2].plot(selectORAMeridO.vomecrty.sel(time_counter=61),np.negative(selectORAMeridO.depthv), label = r'ORA5 (48.75$\degree$N, 16.25$\degree$W)')
# # ax[2].plot(selectORAMerid1.vomecrty.sel(time_counter=61),np.negative(selectORAMerid1.depthv), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# # ax[2].plot(selectORAMerid2.vomecrty.sel(time_counter=61),np.negative(selectORAMerid2.depthv), linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.20 $\degree$W)')
# # ax[2].plot(selectORAMerid3.vomecrty.sel(time_counter=61),np.negative(selectORAMerid3.depthv), linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.20 $\degree$W)')
# ax[2].plot(kppNov12_v,GTMdepth, label = 'GOTM KPP scheme')
# ax[2].plot(kepsNov12_v,GTMdepth, label = r'GOTM k-$\epsilon$ scheme')
# ax[2].grid(); ax[2].legend(bbox_to_anchor=(1, 0.90),loc='center left', fontsize='small')
# ax[2].set_ylim(-350,0);# ax[2].set_xlim(7,19);
# ax[2].tick_params(axis='both', which='major', labelsize=12)


# fig.text(0.5, 0.01, r'Meridonal Velocity - $m s^{-1}$', ha='center', va='center', fontsize=12)
# fig.text(0.008, 0.5, r'Ocean Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')

# ax[0].text(0.025, 0.05, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
# ax[0].text(0.9, 0.045, 'SEP', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[0].transAxes)

# ax[1].text(0.025, 0.05, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[1].text(0.9, 0.045, 'OCT', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[1].transAxes)

# ax[2].text(0.025, 0.05, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
# ax[2].text(0.9, 0.045, 'NOV', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[2].transAxes)

# plt.show()

#WATERFLUX
plt.plot(selectORAwFluxO.time_counter,selectORAwFluxO.sowaflup, label = r'ORA5 (48.7$\degree$N, 16.2$\degree$W)')
# plt.plot(selectORAwFlux1.time_counter,selectORAwFlux1.sowaflup, linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W)')
# plt.plot(selectORAwFlux2.time_counter, selectORAwFlux2.sowaflup, linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W)')
# plt.plot(selectORAwFlux3.time_counter, selectORAwFlux3.sowaflup, linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W)')
# plt.plot(selectORAwFlux4.time_counter, selectORAwFlux4.sowaflup, linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W)')
plt.legend(bbox_to_anchor=(1, 0.90), loc='lower right', fontsize='small')
plt.ylabel(r'Water flux - $Kg \ m^{-2} \ s^{-1}$', fontsize=12); plt.xlabel(r'Time -days since 2012-09-16 00:00:00', fontsize=12                                                                                      )
# plt.ylim(-1000,0); plt.xlim(35.1,35.8)
plt.grid()
plt.show()
exit()


# #  ################################  # #
# 3-MONTHLY PANELS FOR VARIABLE PROFILES #
# #  ################################  # #
# !! disable lines 167-178 for plot the following figures #
