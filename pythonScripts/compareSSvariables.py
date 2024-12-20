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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, 
    LogLocator, ScalarFormatter,LogFormatter)

# plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

#############################################################################################
### -------------------------------------- OSMOSIS -------------------------------------- ###
#############################################################################################
glider = xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; GldrDens = glider.pot_den; 
Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)
# 
new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days
glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})
# print('time', new_time[0:100].values)

#-------------------------------#
##### averaged monthly data #####
#-------------------------------#

## September 2012
gliderSep12 = np.where((new_time > 244) & (new_time <= 274)); 
Sep12_time = new_time[gliderSep12]; Sep12_time = Sep12_time.mean(axis=0); 
Sep12_temp = GldrTemp[0:276,2]; Sep12_temp = Sep12_temp.mean(axis=0)
Sep12_sal = GldrSalt[0:276,2]; Sep12_sal = Sep12_sal.mean(axis=0)
## October 2012
gliderOct12 = np.where((new_time > 274) & (new_time <= 305)); 
Oct12_time = new_time[gliderOct12]; Oct12_time = Oct12_time.mean(axis=0); 
Oct12_temp = GldrTemp[277:588,2]; Oct12_temp = Oct12_temp.mean(axis=0)
Oct12_sal = GldrSalt[277:588,2]; Oct12_sal = Oct12_sal.mean(axis=0)
## November 2012
gliderNov12 = np.where((new_time > 305) & (new_time <= 335)); 
Nov12_time = new_time[gliderNov12]; Nov12_time = Nov12_time.mean(axis=0); 
Nov12_temp = GldrTemp[589:894,2]; Nov12_temp = Nov12_temp.mean(axis=0)
Nov12_sal = GldrSalt[589:894,2]; Nov12_sal = Nov12_sal.mean(axis=0)
## December 2012
gliderDec12 = np.where((new_time > 335) & (new_time <= 366)); 
Dec12_time = new_time[gliderDec12]; Dec12_time = Dec12_time.mean(axis=0); 
Dec12_temp = GldrTemp[895:1217,2]; Dec12_temp = Dec12_temp.mean(axis=0)
Dec12_sal = GldrSalt[895:1217,2]; Dec12_sal = Dec12_sal.mean(axis=0)
## January 2013
gliderJan13 = np.where((new_time > 366) & (new_time <= 397)); 
Jan13_time = new_time[gliderJan13]; Jan13_time = Jan13_time.mean(axis=0); 
Jan13_temp = GldrTemp[1218:1552,2]; Jan13_temp = Jan13_temp.mean(axis=0)
Jan13_sal = GldrSalt[1218:1552,2]; Jan13_sal = Jan13_sal.mean(axis=0)
## February 2013
gliderFeb13 = np.where((new_time > 397) & (new_time <= 456)); 
Feb13_time = new_time[gliderFeb13]; Feb13_time = Feb13_time.mean(axis=0); 
Feb13_temp = GldrTemp[1553:2248,2]; Feb13_temp = Feb13_temp.mean(axis=0)
Feb13_sal = GldrSalt[1553:2248,2]; Feb13_sal = Feb13_sal.mean(axis=0)
## March 2013
gliderMar13 = np.where((new_time > 456) & (new_time <= 487)); 
Mar13_time = new_time[gliderMar13]; Mar13_time = Mar13_time.mean(axis=0); 
Mar13_temp = GldrTemp[2249:2632,2]; Mar13_temp = Mar13_temp.mean(axis=0)
Mar13_sal = GldrSalt[2249:2632,2]; Mar13_sal = Mar13_sal.mean(axis=0)
## April 2013
gliderApr13 = np.where((new_time > 487) & (new_time <= 517)); 
Apr13_time = new_time[gliderApr13]; Apr13_time = Apr13_time.mean(axis=0); 
Apr13_temp = GldrTemp[2633:2992,2]; Apr13_temp = Apr13_temp.mean(axis=0)
Apr13_sal = GldrSalt[2633:2992,2]; Apr13_sal = Apr13_sal.mean(axis=0)
## May 2013
gliderMay13 = np.where((new_time > 517) & (new_time <= 548)); 
May13_time = new_time[gliderMay13]; May13_time = May13_time.mean(axis=0); 
May13_temp = GldrTemp[2993:3356,2]; May13_temp = May13_temp.mean(axis=0)
May13_sal = GldrSalt[2993:3356,2]; May13_sal = May13_sal.mean(axis=0)
## June 2013
gliderJun13 = np.where((new_time > 548) & (new_time <= 578)); 
Jun13_time = new_time[gliderJun13]; Jun13_time = Jun13_time.mean(axis=0); 
Jun13_temp = GldrTemp[3357:3682,2]; Jun13_temp = Jun13_temp.mean(axis=0)
Jun13_sal = GldrSalt[3357:3682,2]; Jun13_sal = Jun13_sal.mean(axis=0)
## July 2013
gliderJul13 = np.where((new_time > 578) & (new_time <= 609)); 
Jul13_time = new_time[gliderJul13]; Jul13_time = Jul13_time.mean(axis=0); 
Jul13_temp = GldrTemp[3683:4000,2]; Jul13_temp = Jul13_temp.mean(axis=0)
Jul13_sal = GldrSalt[3683:4000,2]; Jul13_sal = Jul13_sal.mean(axis=0)
## August 2013
gliderAug13 = np.where((new_time > 609) & (new_time <= 640)); 
Aug13_time = new_time[gliderAug13]; Aug13_time = Aug13_time.mean(axis=0); 
Aug13_temp = GldrTemp[4001:4095,2]; Aug13_temp = Aug13_temp.mean(axis=0)
Aug13_sal = GldrSalt[4001:4095,2]; Aug13_sal = Aug13_sal.mean(axis=0)

gliderOS_SST = np.c_[Sep12_temp,Oct12_temp,Nov12_temp,Dec12_temp,
    Jan13_temp,Feb13_temp,Mar13_temp,Apr13_temp,May13_temp,Jun13_temp,Jul13_temp,Aug13_temp]
gliderOS_SSS = np.c_[Sep12_sal,Oct12_sal,Nov12_sal,Dec12_sal,
    Jan13_sal,Feb13_sal,Mar13_sal,Apr13_sal,May13_sal,Jun13_sal,Jul13_sal,Aug13_sal]
# gliderOS_time = np.c_[Sep12_time,Oct12_time,Nov12_time,Dec12_time,
    # Jan13_time,Feb13_time,Mar13_time,Apr13_time,May13_time,Jun13_time,Jul13_time,Aug13_time]
gliderOS_time = np.c_[0,30,61,91,122,152,181,212,242,273,303,334]

gliderOS_SST = np.array(gliderOS_SST[0])
gliderOS_SSS = np.array(gliderOS_SSS[0])
gliderOS_time = np.asarray(gliderOS_time)
# print('glider time', gliderOS_time)
# print('glider SST', gliderOS_SST)
# exit()

# ### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
gotmkeps = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_keps.nc", decode_times=False)

# KPP variables
GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMsalt = gotmkpp.salt; #
GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMu = gotmkpp.u; GTMv = gotmkpp.v; 
GTMsst = gotmkpp.sst ; GTMsss = gotmkpp.sss;
time = gotmkpp.time


# k-eps variables
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkesalt = gotmkeps.salt; 
GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; 
GTMkesst = gotmkeps.sst ; GTMkesss = gotmkeps.sss;
ketime = gotmkeps.time


## Convert time
time2 = 261 + (time[:]/86400);
ketime2 = 261 + (ketime[:]/86400)

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})

## September 2012
kppSep12 = np.where((time2 > 244) & (time2 <= 274)); kppSep12_time = time2[kppSep12]; 
kppSep12_temp = GTMtemp.isel(lat=0,lon=0)[0:1871,0]; kppSep12_temp = kppSep12_temp.mean(axis=0);
kppSep12_sal = GTMsalt.isel(lat=0,lon=0)[0:1871,0]; kppSep12_sal = kppSep12_sal.mean(axis=0);

kepsSep12 = np.where((ketime2 > 244) & (ketime2 <= 274)); kepsSep12_time = ketime2[kepsSep12]; 
kepsSep12_temp = GTMketemp.isel(lat=0,lon=0)[0:1871,0]; kepsSep12_temp = kepsSep12_temp.mean(axis=0);
kepsSep12_sal = GTMkesalt.isel(lat=0,lon=0)[0:1871,0]; kepsSep12_sal = kepsSep12_sal.mean(axis=0);

## October 2012
kppOct12 = np.where((time2 > 274) & (time2 <= 305)); kppOct12_time = time2[kppOct12]; 
kppOct12_temp = GTMtemp.isel(lat=0,lon=0)[1872:6335,0]; kppOct12_temp = kppOct12_temp.mean(axis=0);
kppOct12_sal = GTMsalt.isel(lat=0,lon=0)[1872:6335,0]; kppOct12_sal = kppOct12_sal.mean(axis=0); 

kepsOct12 = np.where((ketime2 > 274) & (ketime2 <= 305)); kepsOct12_time = ketime2[kepsOct12]; 
kepsOct12_temp = GTMketemp.isel(lat=0,lon=0)[1872:6335,0]; kepsOct12_temp = kepsOct12_temp.mean(axis=0);
kepsOct12_sal = GTMkesalt.isel(lat=0,lon=0)[1872:6335,0]; kepsOct12_sal = kepsOct12_sal.mean(axis=0);

## November 2012
kppNov12 = np.where((time2 > 305) & (time2 <= 335)); kppNov12_time = time2[kppNov12]; 
kppNov12_temp = GTMtemp.isel(lat=0,lon=0)[6336:10655,0]; kppNov12_temp = kppNov12_temp.mean(axis=0);
kppNov12_sal = GTMsalt.isel(lat=0,lon=0)[6336:10655,0]; kppNov12_sal = kppNov12_sal.mean(axis=0);

kepsNov12 = np.where((ketime2 > 305) & (ketime2 <= 335)); kepsNov12_time = ketime2[kepsNov12]; 
kepsNov12_temp = GTMketemp.isel(lat=0,lon=0)[6336:10655,0]; kepsNov12_temp = kepsNov12_temp.mean(axis=0);
kepsNov12_sal = GTMkesalt.isel(lat=0,lon=0)[6336:10655,0]; kepsNov12_sal = kepsNov12_sal.mean(axis=0);

## December 2012
kppDec12 = np.where((time2 > 335) & (time2 <= 366)); kppDec12_time = time2[kppDec12]; 
kppDec12_temp = GTMtemp.isel(lat=0,lon=0)[10656:15119,0]; kppDec12_temp = kppDec12_temp.mean(axis=0);
kppDec12_sal = GTMsalt.isel(lat=0,lon=0)[10656:15119,0]; kppDec12_sal = kppDec12_sal.mean(axis=0);

kepsDec12 = np.where((ketime2 > 335) & (ketime2 <= 366)); kepsDec12_time = ketime2[kepsDec12]; 
kepsDec12_temp = GTMketemp.isel(lat=0,lon=0)[10656:15119,0]; kepsDec12_temp = kepsDec12_temp.mean(axis=0);
kepsDec12_sal = GTMkesalt.isel(lat=0,lon=0)[10656:15119,0]; kepsDec12_sal = kepsDec12_sal.mean(axis=0);

## January 2013
kppJan13 = np.where((time2 > 366) & (time2 <= 397)); kppJan13_time = time2[kppJan13]; 
kppJan13_temp = GTMtemp.isel(lat=0,lon=0)[15120:19583,0]; kppJan13_temp = kppJan13_temp.mean(axis=0);
kppJan13_sal = GTMsalt.isel(lat=0,lon=0)[15120:19583,0]; kppJan13_sal = kppJan13_sal.mean(axis=0);

kepsJan13 = np.where((ketime2 > 366) & (ketime2 <= 397)); kepsJan13_time = ketime2[kepsJan13]; 
kepsJan13_temp = GTMketemp.isel(lat=0,lon=0)[15120:19583,0]; kepsJan13_temp = kepsJan13_temp.mean(axis=0);
kepsJan13_sal = GTMkesalt.isel(lat=0,lon=0)[15120:19583,0]; kepsJan13_sal = kepsJan13_sal.mean(axis=0);

## February 2013
kppFeb13 = np.where((time2 > 397) & (time2 <= 456)); kppFeb13_time = time2[kppFeb13]; 
kppFeb13_temp = GTMtemp.isel(lat=0,lon=0)[19584:28079,0]; kppFeb13_temp = kppFeb13_temp.mean(axis=0);
kppFeb13_sal = GTMsalt.isel(lat=0,lon=0)[19584:28079,0]; kppFeb13_sal = kppFeb13_sal.mean(axis=0);

kepsFeb13 = np.where((ketime2 > 397) & (ketime2 <= 456)); kepsFeb13_time = ketime2[kepsFeb13]; 
kepsFeb13_temp = GTMketemp.isel(lat=0,lon=0)[19584:28079,0]; kepsFeb13_temp = kepsFeb13_temp.mean(axis=0);
kepsFeb13_sal = GTMkesalt.isel(lat=0,lon=0)[19584:28079,0]; kepsFeb13_sal = kepsFeb13_sal.mean(axis=0)

## March 2013
kppMar13 = np.where((time2 > 456) & (time2 <= 487)); kppMar13_time = time2[kppMar13]; 
kppMar13_temp = GTMtemp.isel(lat=0,lon=0)[28080:32543,0]; kppMar13_temp = kppMar13_temp.mean(axis=0);
kppMar13_sal = GTMsalt.isel(lat=0,lon=0)[28080:32543,0]; kppMar13_sal = kppMar13_sal.mean(axis=0);

kepsMar13 = np.where((ketime2 > 456) & (ketime2 <= 487)); kepsMar13_time = ketime2[kepsMar13]; 
kepsMar13_temp = GTMketemp.isel(lat=0,lon=0)[28080:32543,0]; kepsMar13_temp = kepsMar13_temp.mean(axis=0);
kepsMar13_sal = GTMkesalt.isel(lat=0,lon=0)[28080:32543,0]; kepsMar13_sal = kepsMar13_sal.mean(axis=0);

## April 2013
kppApr13 = np.where((time2 > 487) & (time2 <= 517)); kppApr13_time = time2[kppApr13]; 
kppApr13_temp = GTMtemp.isel(lat=0,lon=0)[32544:36863,0]; kppApr13_temp = kppApr13_temp.mean(axis=0);
kppApr13_sal = GTMsalt.isel(lat=0,lon=0)[32544:36863,0]; kppApr13_sal = kppApr13_sal.mean(axis=0);

kepsApr13 = np.where((ketime2 > 487) & (ketime2 <= 517)); kepsApr13_time = ketime2[kepsApr13]; 
kepsApr13_temp = GTMketemp.isel(lat=0,lon=0)[32544:36863,0]; kepsApr13_temp = kepsApr13_temp.mean(axis=0);
kepsApr13_sal = GTMkesalt.isel(lat=0,lon=0)[32544:36863,0]; kepsApr13_sal = kepsApr13_sal.mean(axis=0);

## May 2013
kppMay13 = np.where((time2 > 517) & (time2 <= 548)); kppMay13_time = time2[kppMay13]; 
kppMay13_temp = GTMtemp.isel(lat=0,lon=0)[36864:41327,0]; kppMay13_temp = kppMay13_temp.mean(axis=0);
kppMay13_sal = GTMsalt.isel(lat=0,lon=0)[36864:41327,0]; kppMay13_sal = kppMay13_sal.mean(axis=0);

kepsMay13 = np.where((ketime2 > 517) & (ketime2 <= 548)); kepsMay13_time = ketime2[kepsMay13]; 
kepsMay13_temp = GTMketemp.isel(lat=0,lon=0)[36864:41327,0]; kepsMay13_temp = kepsMay13_temp.mean(axis=0);
kepsMay13_sal = GTMkesalt.isel(lat=0,lon=0)[36864:41327,0]; kepsMay13_sal = kepsMay13_sal.mean(axis=0);

## June 2013
kppJun13 = np.where((time2 > 548) & (time2 <= 578)); kppJun13_time = time2[kppJun13]; 
kppJun13_temp = GTMtemp.isel(lat=0,lon=0)[41328:45647,0]; kppJun13_temp = kppJun13_temp.mean(axis=0);
kppJun13_sal = GTMsalt.isel(lat=0,lon=0)[41328:45647,0]; kppJun13_sal = kppJun13_sal.mean(axis=0);

kepsJun13 = np.where((ketime2 > 548) & (ketime2 <= 578)); kepsJun13_time = ketime2[kepsJun13]; 
kepsJun13_temp = GTMketemp.isel(lat=0,lon=0)[41328:45647,0]; kepsJun13_temp = kepsJun13_temp.mean(axis=0);
kepsJun13_sal = GTMkesalt.isel(lat=0,lon=0)[41328:45647,0]; kepsJun13_sal = kepsJun13_sal.mean(axis=0);

## July 2013
kppJul13 = np.where((time2 > 578) & (time2 <= 609)); kppJul13_time = time2[kppJul13]; 
kppJul13_temp = GTMtemp.isel(lat=0,lon=0)[45648:50111,0]; kppJul13_temp = kppJul13_temp.mean(axis=0);
kppJul13_sal = GTMsalt.isel(lat=0,lon=0)[45648:50111,0]; kppJul13_sal = kppJul13_sal.mean(axis=0);

kepsJul13 = np.where((ketime2 > 578) & (ketime2 <= 609)); kepsJul13_time = ketime2[kepsJul13]; 
kepsJul13_temp = GTMketemp.isel(lat=0,lon=0)[45648:50111,0]; kepsJul13_temp = kepsJul13_temp.mean(axis=0);
kepsJul13_sal = GTMkesalt.isel(lat=0,lon=0)[45648:50111,0]; kepsJul13_sal = kepsJul13_sal.mean(axis=0);

## August 2013
kppAug13 = np.where((time2 > 609) & (time2 <= 640)); kppAug13_time = time2[kppAug13]; 
kppAug13_temp = GTMtemp.isel(lat=0,lon=0)[50112:52595,0]; kppAug13_temp = kppAug13_temp.mean(axis=0);
kppAug13_sal = GTMsalt.isel(lat=0,lon=0)[50112:52595,0]; kppAug13_sal = kppAug13_sal.mean(axis=0);

kepsAug13 = np.where((ketime2 > 609) & (ketime2 <= 640)); kepsAug13_time = ketime2[kepsAug13]; 
kepsAug13_temp = GTMketemp.isel(lat=0,lon=0)[50112:52595,0]; kepsAug13_temp = kepsAug13_temp.mean(axis=0);
kepsAug13_sal = GTMkesalt.isel(lat=0,lon=0)[50112:52595,0]; kepsAug13_sal = kepsAug13_sal.mean(axis=0);

KPP_SST = np.c_[kppSep12_temp,kppOct12_temp,kppNov12_temp,kppDec12_temp,kppJan13_temp,kppFeb13_temp,
    kppMar13_temp,kppApr13_temp,kppMay13_temp,kppJun13_temp,kppJul13_temp,kppAug13_temp]
keps_SST = np.c_[kepsSep12_temp,kepsOct12_temp,kepsNov12_temp,kepsDec12_temp,kepsJan13_temp,kepsFeb13_temp,
    kepsMar13_temp,kepsApr13_temp,kepsMay13_temp,kepsJun13_temp,kepsJul13_temp,kepsAug13_temp]

KPP_SSS = np.c_[kppSep12_sal,kppOct12_sal,kppNov12_sal,kppDec12_sal,kppJan13_sal,kppFeb13_sal,
    kppMar13_sal,kppApr13_sal,kppMay13_sal,kppJun13_sal,kppJul13_sal,kppAug13_sal]
keps_SSS = np.c_[kepsSep12_sal,kepsOct12_sal,kepsNov12_sal,kepsDec12_sal,kepsJan13_sal,kepsFeb13_sal,
    kepsMar13_sal,kepsApr13_sal,kepsMay13_sal,kepsJun13_sal,kepsJul13_sal,kepsAug13_sal]
KPP_SST = np.array(KPP_SST[0])
keps_SST = np.array(keps_SST[0])
KPP_SSS = np.array(KPP_SSS[0])
keps_SSS = np.array(keps_SSS[0])
print('kpp temp', KPP_SST)
print('keps temp', keps_SST)
print('kpp salt', KPP_SSS)
print('keps salt', keps_SSS)

#############################################################################################
### ----------------------------- Ocean Reanalysis System 5 ----------------------------- ###
#############################################################################################

ORA5_SST = xr.open_dataset('sosstsst_combined.nc', decode_times=False);
ORA5_SSS = xr.open_dataset('sosaline_combined.nc', decode_times=False);

# # combine netcdf4 files with same predixes
# ORA5_SST = xr.open_mfdataset('sosstsst_control_monthly_highres_2D_*.nc')
# ORA5_SST = ORA5_SST.to_netcdf('sosstsst_combined.nc')
# ORA5_SSS = xr.open_mfdataset('sosaline_control_monthly_highres_2D_*.nc')
# ORA5_SSS = ORA5_SSS.to_netcdf('sosaline_combined.nc')

# print('ORA5 SST dataset', ORA5_SST.variables.items())
# print('ORA5 SSS dataset', ORA5_SSS.variables.items())


#------------------------------------------------#
##### extract data for selected co-ordinates #####
#------------------------------------------------#
#SEA SURFACE TEMPERATURE
def extract(ORA5_SST, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_SST['nav_lon'].values.ravel(), 
        ORA5_SST['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_SST['nav_lon'].shape))
    return ORA5_SST.isel(x=j0, y=i0).squeeze()
selectORASST = extract(ORA5_SST, -16.25, 48.75); print('extracts',selectORASST)

#SEA SURFACE SALINITY
def extract(ORA5_SSS, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_SSS['nav_lon'].values.ravel(), 
        ORA5_SSS['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_SSS['nav_lon'].shape))
    return ORA5_SSS.isel(x=j0, y=i0).squeeze()
selectORASSS = extract(ORA5_SSS, -16.25, 48.75); print('extracts',selectORASSS)


# selectORASST = selectORASST.sel(time_counter=334)
# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORASSS = selectORASSS.sel(time_counter=334)
# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)

# print('glider SST :', gliderOS_SST.values)
# exit()

month = np.array([datetime(2012,9,1),datetime(2012,10,1),datetime(2012,11,1),datetime(2012,12,1),
    datetime(2013,1,1),datetime(2013,2,1),datetime(2013,3,1),datetime(2013,4,1),datetime(2013,5,1),
    datetime(2013,6,1),datetime(2013,7,1),datetime(2013,8,1)])

#--------------------------#
##### plotting figures #####
#--------------------------#
# # SEA SURFACE TEMPERATURE

# fig,ax = plt.subplots(1,1)
# ax.plot(month, gliderOS_SST, label = r'O5MOSIS (48.7 $\degree$N, 16.2 $\degree$W)')
# ax.plot(month, selectORASST.sosstsst, label = r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax.plot(month, KPP_SST, label = 'GOTM KPP scheme')
# ax.plot(month, keps_SST, label = r'GOTM k-$\epsilon$ scheme')
# ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
# # ax.legend(loc='center right')
# # # ax.xaxis.set_major_locator(mdates.YearLocator())
# # # ax.xaxis.set_major_locator(mdates.DateLocator())
# ax.set_xlim(datetime(2012,9,1),datetime(2013, 8, 1))
# plt.ylabel(r'Sea Surface Temperature - $\degree$ C'); plt.xlabel('Time')
# plt.grid()
# plt.show()

# #SALINITY
# fig,ax = plt.subplots(1,1)
# ax.plot(month,gliderOS_SSS, label = r'O5MOSIS (48.7 $\degree$N, 16.2 $\degree$W)')
# ax.plot(month,selectORASSS.sosaline, label = r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax.plot(month, KPP_SSS, label = 'GOTM KPP scheme')
# ax.plot(month, keps_SSS, label = r'GOTM k-$\epsilon$ scheme')
# ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
# # ax.legend(loc='lower right')
# ax.set_xlim(datetime(2012,9,1),datetime(2013, 8, 1))
# plt.ylabel(r'Sea Surface Salinity - PSU'); plt.xlabel('Time')
# # plt.ylim(-1000,0); plt.xlim(35.1,35.8)
# plt.grid()
# plt.show()

###########################
# CORRELATION COEFFICIENT #
###########################

SST_CC = np.corrcoef(gliderOS_SST,selectORASST.sosstsst)
print('SST correlaion OS vs ORA5', SST_CC)

SSS_CC = np.corrcoef(gliderOS_SSS,selectORASSS.sosaline)
print('SSS correlaion OS vs ORA5', SSS_CC)

exit()
# plot options
selectORATemp.plot(selectORATemp.votemper,np.negative(selectORATemp.deptht), 
    label = 'ORA5 time=%d'%selectORATemp.time_counter)