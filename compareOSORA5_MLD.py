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
print('time', new_time[0:100].values)

#########################################
### SECTION OSMOSIS MIXED LAYER DEPTH ###
#########################################

### ----------------------------- mixed layer depth .csv file ----------------------------- ###
def read_lines():
    with open('mld.csv', 'rU') as mldFile:
        reader = csv.reader(mldFile)
        for row in reader:
            yield [ float(i) for i in row ]
for i in read_lines():
    print('mld :',i)

# to get a list, instead of a generator, use
# xy = list(read_lines())

def read__lines():
    with open('dats.csv', 'rU') as datFile:
        reader_ = csv.reader(datFile)
        for row in reader_:
            yield [ float(j) for j in row ]
for j in read__lines():
    print('date :',j)
j = np.array(j)
j = j + 241

## ------------------------------- moving average MLDs--------------------------------- ##

def moving_average(i, n=3) :
    ret = np.cumsum(i, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# i = np.arange(10); print('set :',i); example = moving_average(i); 
# print('moving average example :', example)

df_mld = moving_average(i)
df_mld5 = moving_average(i,n=5)
df_mld7 = moving_average(i,n=7)

df_dats = moving_average(j)
df_dats5 = moving_average(j,n=5)
df_dats7 = moving_average(j,n=7)

## September 2012
gliderSep12 = np.where((df_dats5 > 244) & (df_dats5 <= 274)); Sep12_time = df_dats5[gliderSep12]; 
Sep12_mld = df_mld5[0:285]; Sep12_mld = Sep12_mld.mean(axis=0)
## October 2012
gliderOct12 = np.where((df_dats5 > 274) & (df_dats5 <= 305)); Oct12_time = df_dats5[gliderOct12]; 
Oct12_mld = df_mld5[286:597]; Oct12_mld = Oct12_mld.mean(axis=0)
## November 2012
gliderNov12 = np.where((df_dats5 > 305) & (df_dats5 <= 335)); Nov12_time = df_dats5[gliderNov12]; 
Nov12_mld = df_mld5[598:903]; Nov12_mld = Nov12_mld.mean(axis=0)
## December 2012
gliderDec12 = np.where((df_dats5 > 335) & (df_dats5 <= 366)); Dec12_time = df_dats5[gliderDec12]; 
Dec12_mld = df_mld5[904:1226]; Dec12_mld = Dec12_mld.mean(axis=0)
## January 2013
gliderJan13 = np.where((df_dats5 > 366) & (df_dats5 <= 397)); Jan13_time = df_dats5[gliderJan13]; 
Jan13_mld = df_mld5[1227:1562]; Jan13_mld = Jan13_mld.mean(axis=0)
## February 2013
gliderFeb13 = np.where((df_dats5 > 397) & (df_dats5 <= 456)); Feb13_time = df_dats5[gliderFeb13]; 
Feb13_mld = df_mld5[1563:2258]; Feb13_mld = Feb13_mld.mean(axis=0)
## March 2013
gliderMar13 = np.where((df_dats5 > 456) & (df_dats5 <= 487)); Mar13_time = df_dats5[gliderMar13]; 
Mar13_mld = df_mld5[2259:2644]; Mar13_mld = Mar13_mld.mean(axis=0)
## April 2013
gliderApr13 = np.where((df_dats5 > 487) & (df_dats5 <= 517)); Apr13_time = df_dats5[gliderApr13]; 
Apr13_mld = df_mld5[2645:3003]; Apr13_mld = Apr13_mld.mean(axis=0)
## May 2013
gliderMay13 = np.where((df_dats5 > 517) & (df_dats5 <= 548)); May13_time = df_dats5[gliderMay13]; 
May13_mld = df_mld5[3004:3366]; May13_mld = May13_mld.mean(axis=0)
## June 2013
gliderJun13 = np.where((df_dats5 > 548) & (df_dats5 <= 578)); Jun13_time = df_dats5[gliderJun13]; 
Jun13_mld = df_mld5[3367:3690]; Jun13_mld = Jun13_mld.mean(axis=0)
## July 2013
gliderJul13 = np.where((df_dats5 > 578) & (df_dats5 <= 609)); Jul13_time = df_dats5[gliderJul13]; 
Jul13_mld = df_mld5[3691:4008]; Jul13_mld = Jul13_mld.mean(axis=0)
## August 2013
gliderAug13 = np.where((df_dats5 > 609) & (df_dats5 <= 640)); Aug13_time = df_dats5[gliderAug13]; 
Aug13_mld = df_mld5[4009:4091]; Aug13_mld = Aug13_mld.mean(axis=0)


gliderOS_MLD = np.c_[Sep12_mld,Oct12_mld,Nov12_mld,Dec12_mld,Jan13_mld,Feb13_mld,Mar13_mld,Apr13_mld,May13_mld,Jun13_mld,Jul13_mld,Aug13_mld]
# gliderOS_time = np.c_[Sep12_time,Oct12_time,Nov12_time,Dec12_time,Jan13_time,Feb13_time,Mar13_time,Apr13_time,May13_time,Jun13_time,Jul13_time,Aug13_time]
gliderOS_time = np.c_[0,30,61,91,122,152,181,212,242,273,303,334]

gliderOS_MLD = np.array(gliderOS_MLD[0])
gliderOS_time = np.asarray(gliderOS_time)
# print('glider time', gliderOS_time)
# print('glider MLD', gliderOS_MLD)
# exit()

# ########################################################################
# ### ----------------------------- GOTM ----------------------------- ###
# ########################################################################
gotmkpp = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
gotmkeps = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_keps.nc", decode_times=False)

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMsalt = gotmkpp.salt; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMu = gotmkpp.u; GTMv = gotmkpp.v; time = gotmkpp.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkesalt = gotmkeps.salt; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time

## Convert time
time2 = 261 + (time[:]/86400);
ketime2 = 261 + (ketime[:]/86400)

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})

## September 2012
kppSep12 = np.where((time2 > 244) & (time2 <= 274)); kppSep12_time = time2[kppSep12]; 
kppSep12_temp = GTMtemp[0:1871,:]; kppSep12_temp = kppSep12_temp.mean(axis=0); kppSep12_temp = kppSep12_temp.isel(lat=0,lon=0)
kppSep12_sal = GTMsalt[0:1871,:]; kppSep12_sal = kppSep12_sal.mean(axis=0); kppSep12_sal = kppSep12_sal.isel(lat=0,lon=0)

kepsSep12 = np.where((ketime2 > 244) & (ketime2 <= 274)); kepsSep12_time = ketime2[kepsSep12]; 
kepsSep12_temp = GTMketemp[0:1871,:]; kepsSep12_temp = kepsSep12_temp.mean(axis=0); kepsSep12_temp = kepsSep12_temp.isel(lat=0,lon=0)
kepsSep12_sal = GTMkesalt[0:1871,:]; kepsSep12_sal = kepsSep12_sal.mean(axis=0); kepsSep12_sal = kepsSep12_sal.isel(lat=0, lon=0)

## October 2012
kppOct12 = np.where((time2 > 274) & (time2 <= 305)); kppOct12_time = time2[kppOct12]; 
kppOct12_temp = GTMtemp[1872:6335,:]; kppOct12_temp = kppOct12_temp.mean(axis=0); kppOct12_temp = kppOct12_temp.isel(lat=0, lon=0)
kppOct12_sal = GTMsalt[1872:6335,:]; kppOct12_sal = kppOct12_sal.mean(axis=0); kppOct12_sal = kppOct12_sal.isel(lat=0, lon=0)

kepsOct12 = np.where((ketime2 > 274) & (ketime2 <= 305)); kepsOct12_time = ketime2[kepsOct12]; 
kepsOct12_temp = GTMketemp[1872:6335,:]; kepsOct12_temp = kepsOct12_temp.mean(axis=0); kepsOct12_temp = kepsOct12_temp.isel(lat=0, lon=0)
kepsOct12_sal = GTMkesalt[1872:6335,:]; kepsOct12_sal = kepsOct12_sal.mean(axis=0); kepsOct12_sal = kepsOct12_sal.isel(lat=0, lon=0)

## November 2012
kppNov12 = np.where((time2 > 305) & (time2 <= 335)); kppNov12_time = time2[kppNov12]; 
kppNov12_temp = GTMtemp[6336:10655,:]; kppNov12_temp = kppNov12_temp.mean(axis=0); kppNov12_temp = kppNov12_temp.isel(lat=0, lon=0)
kppNov12_sal = GTMsalt[6336:10655,:]; kppNov12_sal = kppNov12_sal.mean(axis=0); kppNov12_sal = kppNov12_sal.isel(lat=0, lon=0)

kepsNov12 = np.where((ketime2 > 305) & (ketime2 <= 335)); kepsNov12_time = ketime2[kepsNov12]; 
kepsNov12_temp = GTMketemp[6336:10655,:]; kepsNov12_temp = kepsNov12_temp.mean(axis=0); kepsNov12_temp = kepsNov12_temp.isel(lat=0, lon=0)
kepsNov12_sal = GTMkesalt[6336:10655,:]; kepsNov12_sal = kepsNov12_sal.mean(axis=0); kepsNov12_sal = kepsNov12_sal.isel(lat=0, lon=0)

## December 2012
kppDec12 = np.where((time2 > 335) & (time2 <= 366)); kppDec12_time = time2[kppDec12]; 
kppDec12_temp = GTMtemp[10656:15119,:]; kppDec12_temp = kppDec12_temp.mean(axis=0); kppDec12_temp = kppDec12_temp.isel(lat=0, lon=0)
kppDec12_sal = GTMsalt[10656:15119,:]; kppDec12_sal = kppDec12_sal.mean(axis=0); kppDec12_sal = kppDec12_sal.isel(lat=0, lon=0)

kepsDec12 = np.where((ketime2 > 335) & (ketime2 <= 366)); kepsDec12_time = ketime2[kepsDec12]; 
kepsDec12_temp = GTMtemp[10656:15119,:]; kepsDec12_temp = kepsDec12_temp.mean(axis=0); kepsDec12_temp = kepsDec12_temp.isel(lat=0, lon=0)
kepsDec12_sal = GTMkesalt[10656:15119,:]; kepsDec12_sal = kepsDec12_sal.mean(axis=0); kepsDec12_sal = kepsDec12_sal.isel(lat=0, lon=0)


## January 2013
kppJan13 = np.where((time2 > 366) & (time2 <= 397)); kppJan13_time = time2[kppJan13]; 
kppJan13_temp = GTMtemp[15120:19583,:]; kppJan13_temp = kppJan13_temp.mean(axis=0); kppJan13_temp = kppJan13_temp.isel(lat=0, lon=0)
kppJan13_sal = GTMsalt[15120:19583,:]; kppJan13_sal = kppJan13_sal.mean(axis=0); kppJan13_sal = kppJan13_sal.isel(lat=0, lon=0)

kepsJan13 = np.where((ketime2 > 366) & (ketime2 <= 397)); kepsJan13_time = ketime2[kepsJan13]; 
kepsJan13_temp = GTMketemp[15120:19583,:]; kepsJan13_temp = kepsJan13_temp.mean(axis=0); kepsJan13_temp = kepsJan13_temp.isel(lat=0, lon=0)
kepsJan13_sal = GTMkesalt[15120:19583,:]; kepsJan13_sal = kepsJan13_sal.mean(axis=0); kepsJan13_sal = kepsJan13_sal.isel(lat=0, lon=0)

## February 2013
kppFeb13 = np.where((time2 > 397) & (time2 <= 456)); kppFeb13_time = time2[kppFeb13]; 
kppFeb13_temp = GTMtemp[19584:28079,:]; kppFeb13_temp = kppFeb13_temp.mean(axis=0); kppFeb13_temp = kppFeb13_temp.isel(lat=0, lon=0)
kppFeb13_sal = GTMsalt[19584:28079,:]; kppFeb13_sal = kppFeb13_sal.mean(axis=0); kppFeb13_sal = kppFeb13_sal.isel(lat=0, lon=0)

kepsFeb13 = np.where((ketime2 > 397) & (ketime2 <= 456)); kepsFeb13_time = ketime2[kepsFeb13]; 
kepsFeb13_temp = GTMketemp[19584:28079,:]; kepsFeb13_temp = kepsFeb13_temp.mean(axis=0); kepsFeb13_temp = kepsFeb13_temp.isel(lat=0, lon=0)
kepsFeb13_sal = GTMkesalt[19584:28079,:]; kepsFeb13_sal = kepsFeb13_sal.mean(axis=0)

## March 2013
kppMar13 = np.where((time2 > 456) & (time2 <= 487)); kppMar13_time = time2[kppMar13]; 
kppMar13_temp = GTMtemp[28080:32543,:]; kppMar13_temp = kppMar13_temp.mean(axis=0); kppMar13_temp = kppMar13_temp.isel(lat=0, lon=0)
kppMar13_sal = GTMsalt[28080:32543,:]; kppMar13_sal = kppMar13_sal.mean(axis=0); kppMar13_sal = kppMar13_sal.isel(lat=0, lon=0)

kepsMar13 = np.where((ketime2 > 456) & (ketime2 <= 487)); kepsMar13_time = ketime2[kepsMar13]; 
kepsMar13_temp = GTMketemp[28080:32543,:]; kepsMar13_temp = kepsMar13_temp.mean(axis=0); kepsMar13_temp = kepsMar13_temp.isel(lat=0, lon=0)
kepsMar13_sal = GTMkesalt[28080:32543,:]; kepsMar13_sal = kepsMar13_sal.mean(axis=0); kepsMar13_sal = kepsMar13_sal.isel(lat=0, lon=0)

## April 2013
kppApr13 = np.where((time2 > 487) & (time2 <= 517)); kppApr13_time = time2[kppApr13]; 
kppApr13_temp = GTMtemp[32544:36863,:]; kppApr13_temp = kppApr13_temp.mean(axis=0); kppApr13_temp = kppApr13_temp.isel(lat=0, lon=0)
kppApr13_sal = GTMsalt[32544:36863,:]; kppApr13_sal = kppApr13_sal.mean(axis=0); kppApr13_sal = kppApr13_sal.isel(lat=0, lon=0)

kepsApr13 = np.where((ketime2 > 487) & (ketime2 <= 517)); kepsApr13_time = ketime2[kepsApr13]; 
kepsApr13_temp = GTMketemp[32544:36863,:]; kepsApr13_temp = kepsApr13_temp.mean(axis=0); kepsApr13_temp = kepsApr13_temp.isel(lat=0, lon=0)
kepsApr13_sal = GTMkesalt[32544:36863,:]; kepsApr13_sal = kepsApr13_sal.mean(axis=0); kepsApr13_sal = kepsApr13_sal.isel(lat=0,lon=0)

## May 2013
kppMay13 = np.where((time2 > 517) & (time2 <= 548)); kppMay13_time = time2[kppMay13]; 
kppMay13_temp = GTMtemp[36864:41327,:]; kppMay13_temp = kppMay13_temp.mean(axis=0); kppMay13_temp = kppMay13_temp.isel(lat=0, lon=0)
kppMay13_sal = GTMsalt[36864:41327,:]; kppMay13_sal = kppMay13_sal.mean(axis=0); kppMay13_sal = kppMay13_sal.isel(lat=0, lon=0)

kepsMay13 = np.where((ketime2 > 517) & (ketime2 <= 548)); kepsMay13_time = ketime2[kepsMay13]; 
kepsMay13_temp = GTMketemp[36864:41327,:]; kepsMay13_temp = kepsMay13_temp.mean(axis=0); kepsMay13_temp = kepsMay13_temp.isel(lat=0, lon=0)
kepsMay13_sal = GTMkesalt[36864:41327,:]; kepsMay13_sal = kepsMay13_sal.mean(axis=0); kepsMay13_sal = kepsMay13_sal.isel(lat=0, lon=0)

## June 2013
kppJun13 = np.where((time2 > 548) & (time2 <= 578)); kppJun13_time = time2[kppJun13]; 
kppJun13_temp = GTMtemp[41328:45647,:]; kppJun13_temp = kppJun13_temp.mean(axis=0); kppJun13_temp = kppJun13_temp.isel(lat=0, lon=0)
kppJun13_sal = GTMsalt[41328:45647,:]; kppJun13_sal = kppJun13_sal.mean(axis=0); kppJun13_sal = kppJun13_sal.isel(lat=0, lon=0)

kepsJun13 = np.where((ketime2 > 548) & (ketime2 <= 578)); kepsJun13_time = ketime2[kepsJun13]; 
kepsJun13_temp = GTMketemp[41328:45647,:]; kepsJun13_temp = kepsJun13_temp.mean(axis=0); kepsJun13_temp = kepsJun13_temp.isel(lat=0, lon=0)
kepsJun13_sal = GTMkesalt[41328:45647,:]; kepsJun13_sal = kepsJun13_sal.mean(axis=0); kepsJun13_sal = kepsJun13_sal.isel(lat=0, lon=0)

## July 2013
kppJul13 = np.where((time2 > 578) & (time2 <= 609)); kppJul13_time = time2[kppJul13]; 
kppJul13_temp = GTMtemp[45648:50111,:]; kppJul13_temp = kppJul13_temp.mean(axis=0); kppJul13_temp = kppJul13_temp.isel(lat=0, lon=0)
kppJul13_sal = GTMsalt[45648:50111,:]; kppJul13_sal = kppJul13_sal.mean(axis=0); kppJul13_sal = kppJul13_sal.isel(lat=0, lon=0)

kepsJul13 = np.where((ketime2 > 578) & (ketime2 <= 609)); kepsJul13_time = ketime2[kepsJul13]; 
kepsJul13_temp = GTMketemp[45648:50111,:]; kepsJul13_temp = kepsJul13_temp.mean(axis=0); kepsJul13_temp = kepsJul13_temp.isel(lat=0, lon=0)
kepsJul13_sal = GTMkesalt[45648:50111,:]; kepsJul13_sal = kepsJul13_sal.mean(axis=0); kepsJul13_sal = kepsJul13_sal.isel(lat=0, lon=0)

## August 2013
kppAug13 = np.where((time2 > 609) & (time2 <= 640)); kppAug13_time = time2[kppAug13]; 
kppAug13_temp = GTMtemp[50112:52595,:]; kppAug13_temp = kppAug13_temp.mean(axis=0); kppAug13_temp = kppAug13_temp.isel(lat=0, lon=0)
kppAug13_sal = GTMsalt[50112:52595,:]; kppAug13_sal = kppAug13_sal.mean(axis=0); kppAug13_sal = kppAug13_sal.isel(lat=0, lon=0)

kepsAug13 = np.where((ketime2 > 609) & (ketime2 <= 640)); kepsAug13_time = ketime2[kepsAug13]; 
kepsAug13_temp = GTMketemp[50112:52595,:]; kepsAug13_temp = kepsAug13_temp.mean(axis=0); kepsAug13_temp = kepsAug13_temp.isel(lat=0, lon=0)
kepsAug13_sal = GTMkesalt[50112:52595,:]; kepsAug13_sal = kepsAug13_sal.mean(axis=0); kepsAug13_sal = kepsAug13_sal.isel(lat=0, lon=0)


#  ----------------------------- calculate MIXED LAYER DEPTH #1 -----------------------------  #

## Psedocode 
# 1. extract temperature values for depth z_ref
# 2. new_temp = temp(z_ref) - delta_t
# 3. depth at new_temp

########################
## original KPP model ##
mld_tempKPP = GTMtemp.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMtemp-mld_tempKPP).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
mldKPP = GTMtemp.z[z_indexes]; 
mldKPP = mldKPP.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

##########################
## original k-eps model ##
mld_tempkeps = GTMketemp.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMketemp-mld_tempkeps).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
mldkeps = GTMketemp.z[z_indexes]; 
mldkeps = mldkeps.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## September 2012
kppSep12MLD = np.where((time2 > 244) & (time2 <= 274)); kppSep12_time = time2[kppSep12MLD]; 
kppSep12_mld = mldKPP[0:1871]; kppSep12_mld = kppSep12_mld.mean(axis=0);

kepsSep12MLD = np.where((ketime2 > 244) & (ketime2 <= 274)); kepsSep12_time = ketime2[kepsSep12MLD]; 
kepsSep12_mld = mldkeps[0:1871]; kepsSep12_mld = kepsSep12_mld.mean(axis=0);

## October 2012
kppOct12MLD = np.where((time2 > 274) & (time2 <= 305)); kppOct12_time = time2[kppOct12MLD]; 
kppOct12_mld = mldKPP[1872:6335]; kppOct12_mld = kppOct12_mld.mean(axis=0);

kepsOct12MLD = np.where((ketime2 > 274) & (ketime2 <= 305)); kepsOct12_time = ketime2[kepsOct12MLD]; 
kepsOct12_mld = mldkeps[1872:6335]; kepsOct12_mld = kepsOct12_mld.mean(axis=0);

## November 2012
kppNov12MLD = np.where((time2 > 305) & (time2 <= 335)); kppNov12_time = time2[kppNov12MLD]; 
kppNov12_mld = mldKPP[6336:10655]; kppNov12_mld = kppNov12_mld.mean(axis=0);

kepsNov12MLD = np.where((ketime2 > 305) & (ketime2 <= 335)); kepsNov12_time = ketime2[kepsNov12MLD]; 
kepsNov12_mld = mldkeps[6336:10655]; kepsNov12_mld = kepsNov12_mld.mean(axis=0);

## December 2012
kppDec12MLD = np.where((time2 > 335) & (time2 <= 366)); kppDec12_time = time2[kppDec12]; 
kppDec12_mld = mldKPP[10656:15119]; kppDec12_mld = kppDec12_mld.mean(axis=0);

kepsDec12MLD = np.where((ketime2 > 335) & (ketime2 <= 366)); kepsDec12_time = ketime2[kepsDec12]; 
kepsDec12_mld = mldkeps[10656:15119]; kepsDec12_mld = kepsDec12_mld.mean(axis=0);

## January 2013
kppJan13MLD = np.where((time2 > 366) & (time2 <= 397)); kppJan13_time = time2[kppJan13MLD]; 
kppJan13_mld = mldKPP[15120:19583]; kppJan13_mld = kppJan13_mld.mean(axis=0);

kepsJan13MLD = np.where((ketime2 > 366) & (ketime2 <= 397)); kepsJan13_time = ketime2[kepsJan13MLD]; 
kepsJan13_mld = mldkeps[15120:19583]; kepsJan13_mld = kepsJan13_mld.mean(axis=0);

## February 2013
kppFeb13MLD = np.where((time2 > 397) & (time2 <= 456)); kppFeb13_time = time2[kppFeb13MLD]; 
kppFeb13_mld = mldKPP[19584:28079]; kppFeb13_mld = kppFeb13_mld.mean(axis=0);

kepsFeb13MLD = np.where((ketime2 > 397) & (ketime2 <= 456)); kepsFeb13_time = ketime2[kepsFeb13MLD]; 
kepsFeb13_mld = mldkeps[19584:28079]; kepsFeb13_mld = kepsFeb13_mld.mean(axis=0);

## March 2013
kppMar13MLD = np.where((time2 > 456) & (time2 <= 487)); kppMar13_time = time2[kppMar13MLD]; 
kppMar13_mld = mldKPP[28080:32543]; kppMar13_mld = kppMar13_mld.mean(axis=0);

kepsMar13MLD = np.where((ketime2 > 456) & (ketime2 <= 487)); kepsMar13_time = ketime2[kepsMar13MLD]; 
kepsMar13_mld = mldkeps[28080:32543]; kepsMar13_mld = kepsMar13_mld.mean(axis=0);

## April 2013
kppApr13MLD = np.where((time2 > 487) & (time2 <= 517)); kppApr13_time = time2[kppApr13MLD]; 
kppApr13_mld = mldKPP[32544:36863]; kppApr13_mld = kppApr13_mld.mean(axis=0);

kepsApr13MLD = np.where((ketime2 > 487) & (ketime2 <= 517)); kepsApr13_time = ketime2[kepsApr13MLD]; 
kepsApr13_mld = mldkeps[32544:36863]; kepsApr13_mld = kepsApr13_mld.mean(axis=0);

## May 2013
kppMay13MLD = np.where((time2 > 517) & (time2 <= 548)); kppMay13_time = time2[kppMay13MLD]; 
kppMay13_mld = mldKPP[36864:41327]; kppMay13_mld = kppMay13_mld.mean(axis=0);

kepsMay13MLD = np.where((ketime2 > 517) & (ketime2 <= 548)); kepsMay13_time = ketime2[kepsMay13MLD]; 
kepsMay13_mld = mldkeps[36864:41327]; kepsMay13_mld = kepsMay13_mld.mean(axis=0);

## June 2013
kppJun13MLD = np.where((time2 > 548) & (time2 <= 578)); kppJun13_time = time2[kppJun13MLD]; 
kppJun13_mld = mldKPP[41328:45647]; kppJun13_mld = kppJun13_mld.mean(axis=0);

kepsJun13MLD = np.where((ketime2 > 548) & (ketime2 <= 578)); kepsJun13_time = ketime2[kepsJun13MLD]; 
kepsJun13_mld = mldkeps[41328:45647]; kepsJun13_mld = kepsJun13_mld.mean(axis=0);

## July 2013
kppJul13MLD = np.where((time2 > 578) & (time2 <= 609)); kppJul13_time = time2[kppJul13MLD]; 
kppJul13_mld = mldKPP[45648:50111]; kppJul13_mld = kppJul13_mld.mean(axis=0);

kepsJul13MLD = np.where((ketime2 > 578) & (ketime2 <= 609)); kepsJul13_time = ketime2[kepsJul13MLD]; 
kepsJul13_mld = mldkeps[45648:50111]; kepsJul13_mld = kepsJul13_mld.mean(axis=0);

## August 2013
kppAug13MLD = np.where((time2 > 609) & (time2 <= 640)); kppAug13_time = time2[kppAug13MLD]; 
kppAug13_mld = mldKPP[50112:52595]; kppAug13_mld = kppAug13_mld.mean(axis=0);

kepsAug13MLD = np.where((ketime2 > 609) & (ketime2 <= 640)); kepsAug13_time = ketime2[kepsAug13MLD]; 
kepsAug13_mld = mldkeps[50112:52595]; kepsAug13_mld = kepsAug13_mld.mean(axis=0);


KPP_MLD = np.c_[kppSep12_mld,kppOct12_mld,kppNov12_mld,kppDec12_mld,kppJan13_mld,kppFeb13_mld,kppMar13_mld,
    kppApr13_mld,kppMay13_mld,kppJun13_mld,kppJul13_mld,kppAug13_mld]
keps_MLD = np.c_[kepsSep12_mld,kepsOct12_mld,kepsNov12_mld,kepsDec12_mld,kepsJan13_mld,kepsFeb13_mld,kepsMar13_mld,
    kepsApr13_mld,kepsMay13_mld,kepsJun13_mld,kepsJul13_mld,kepsAug13_mld]

KPP_MLD = np.array(KPP_MLD[0])
keps_MLD = np.array(keps_MLD[0])



# #  ----------------------------- calculate MIXED LAYER DEPTH #2 -----------------------------  #
# # # new diagnosis to avoid the argmin() component which results in MLD collapses

# ## Psuedocode
# # 1. loop through all levels of depth in the temperature variables from bottom to the top
# # 2. flag where the sign changes from negative to positive
# # 3. establish the position of the element to be MLD 


# MLDkppindx = []; MLDkeindx = []; 
# ########################
# ## use orig KPP model ##
# mld_tempKPP = GTMtemp.sel(z=-11.0, method='nearest') - 0.2;
# mld_tempKPP = (GTMtemp-mld_tempKPP).isel(lat=0,lon=0).values

# for time in mld_tempKPP:
#     for depth in range(len(time)):
#         if time[depth] > 0:
#             MLDkppindx.append(depth)
#             break
# MLDkppindx = np.asarray(MLDkppindx); mldKPP = GTMtemp.z[MLDkppindx];

# ########################
# ## use orig k-eps model ##
# mld_tempkeps = GTMketemp.sel(z=-11.0, method='nearest') - 0.2;
# mld_tempkeps = (GTMketemp-mld_tempkeps).isel(lat=0,lon=0).values

# for time in mld_tempkeps:
#     for depth in range(len(time)):
#         if time[depth] > 0:
#             MLDkeindx.append(depth)
#             break
# MLDkeindx = np.asarray(MLDkeindx); mldkeps = GTMketemp.z[MLDkeindx];

# ## September 2012
# kppSep12MLD = np.where((time2 > 244) & (time2 <= 274)); kppSep12_time = time2[kppSep12MLD]; 
# kppSep12_mld = mldKPP[0:1871]; kppSep12_mld = kppSep12_mld.mean(axis=0);

# kepsSep12MLD = np.where((ketime2 > 244) & (ketime2 <= 274)); kepsSep12_time = ketime2[kepsSep12MLD]; 
# kepsSep12_mld = mldkeps[0:1871]; kepsSep12_mld = kepsSep12_mld.mean(axis=0);

# ## October 2012
# kppOct12MLD = np.where((time2 > 274) & (time2 <= 305)); kppOct12_time = time2[kppOct12MLD]; 
# kppOct12_mld = mldKPP[1872:6335]; kppOct12_mld = kppOct12_mld.mean(axis=0);

# kepsOct12MLD = np.where((ketime2 > 274) & (ketime2 <= 305)); kepsOct12_time = ketime2[kepsOct12MLD]; 
# kepsOct12_mld = mldkeps[1872:6335]; kepsOct12_mld = kepsOct12_mld.mean(axis=0);

# ## November 2012
# kppNov12MLD = np.where((time2 > 305) & (time2 <= 335)); kppNov12_time = time2[kppNov12MLD]; 
# kppNov12_mld = mldKPP[6336:10655]; kppNov12_mld = kppNov12_mld.mean(axis=0);

# kepsNov12MLD = np.where((ketime2 > 305) & (ketime2 <= 335)); kepsNov12_time = ketime2[kepsNov12MLD]; 
# kepsNov12_mld = mldkeps[6336:10655]; kepsNov12_mld = kepsNov12_mld.mean(axis=0);

# ## December 2012
# kppDec12MLD = np.where((time2 > 335) & (time2 <= 366)); kppDec12_time = time2[kppDec12]; 
# kppDec12_mld = mldKPP[10656:15119]; kppDec12_mld = kppDec12_mld.mean(axis=0);

# kepsDec12MLD = np.where((ketime2 > 335) & (ketime2 <= 366)); kepsDec12_time = ketime2[kepsDec12]; 
# kepsDec12_mld = mldkeps[10656:15119]; kepsDec12_mld = kepsDec12_mld.mean(axis=0);

# ## January 2013
# kppJan13MLD = np.where((time2 > 366) & (time2 <= 397)); kppJan13_time = time2[kppJan13MLD]; 
# kppJan13_mld = mldKPP[15120:19583]; kppJan13_mld = kppJan13_mld.mean(axis=0);

# kepsJan13MLD = np.where((ketime2 > 366) & (ketime2 <= 397)); kepsJan13_time = ketime2[kepsJan13MLD]; 
# kepsJan13_mld = mldkeps[15120:19583]; kepsJan13_mld = kepsJan13_mld.mean(axis=0);

# ## February 2013
# kppFeb13MLD = np.where((time2 > 397) & (time2 <= 456)); kppFeb13_time = time2[kppFeb13MLD]; 
# kppFeb13_mld = mldKPP[19584:28079]; kppFeb13_mld = kppFeb13_mld.mean(axis=0);

# kepsFeb13MLD = np.where((ketime2 > 397) & (ketime2 <= 456)); kepsFeb13_time = ketime2[kepsFeb13MLD]; 
# kepsFeb13_mld = mldkeps[19584:28079]; kepsFeb13_mld = kepsFeb13_mld.mean(axis=0);

# ## March 2013
# kppMar13MLD = np.where((time2 > 456) & (time2 <= 487)); kppMar13_time = time2[kppMar13MLD]; 
# kppMar13_mld = mldKPP[28080:32543]; kppMar13_mld = kppMar13_mld.mean(axis=0);

# kepsMar13MLD = np.where((ketime2 > 456) & (ketime2 <= 487)); kepsMar13_time = ketime2[kepsMar13MLD]; 
# kepsMar13_mld = mldkeps[28080:32543]; kepsMar13_mld = kepsMar13_mld.mean(axis=0);

# ## April 2013
# kppApr13MLD = np.where((time2 > 487) & (time2 <= 517)); kppApr13_time = time2[kppApr13MLD]; 
# kppApr13_mld = mldKPP[32544:36863]; kppApr13_mld = kppApr13_mld.mean(axis=0);

# kepsApr13MLD = np.where((ketime2 > 487) & (ketime2 <= 517)); kepsApr13_time = ketime2[kepsApr13MLD]; 
# kepsApr13_mld = mldkeps[32544:36863]; kepsApr13_mld = kepsApr13_mld.mean(axis=0);

# ## May 2013
# kppMay13MLD = np.where((time2 > 517) & (time2 <= 548)); kppMay13_time = time2[kppMay13MLD]; 
# kppMay13_mld = mldKPP[36864:41327]; kppMay13_mld = kppMay13_mld.mean(axis=0);

# kepsMay13MLD = np.where((ketime2 > 517) & (ketime2 <= 548)); kepsMay13_time = ketime2[kepsMay13MLD]; 
# kepsMay13_mld = mldkeps[36864:41327]; kepsMay13_mld = kepsMay13_mld.mean(axis=0);

# ## June 2013
# kppJun13MLD = np.where((time2 > 548) & (time2 <= 578)); kppJun13_time = time2[kppJun13MLD]; 
# kppJun13_mld = mldKPP[41328:45647]; kppJun13_mld = kppJun13_mld.mean(axis=0);

# kepsJun13MLD = np.where((ketime2 > 548) & (ketime2 <= 578)); kepsJun13_time = ketime2[kepsJun13MLD]; 
# kepsJun13_mld = mldkeps[41328:45647]; kepsJun13_mld = kepsJun13_mld.mean(axis=0);

# ## July 2013
# kppJul13MLD = np.where((time2 > 578) & (time2 <= 609)); kppJul13_time = time2[kppJul13MLD]; 
# kppJul13_mld = mldKPP[45648:50111]; kppJul13_mld = kppJul13_mld.mean(axis=0);

# kepsJul13MLD = np.where((ketime2 > 578) & (ketime2 <= 609)); kepsJul13_time = ketime2[kepsJul13MLD]; 
# kepsJul13_mld = mldkeps[45648:50111]; kepsJul13_mld = kepsJul13_mld.mean(axis=0);

# ## August 2013
# kppAug13MLD = np.where((time2 > 609) & (time2 <= 640)); kppAug13_time = time2[kppAug13MLD]; 
# kppAug13_mld = mldKPP[50112:52595]; kppAug13_mld = kppAug13_mld.mean(axis=0);

# kepsAug13MLD = np.where((ketime2 > 609) & (ketime2 <= 640)); kepsAug13_time = ketime2[kepsAug13MLD]; 
# kepsAug13_mld = mldkeps[50112:52595]; kepsAug13_mld = kepsAug13_mld.mean(axis=0);


KPP_MLD = np.c_[kppSep12_mld,kppOct12_mld,kppNov12_mld,kppDec12_mld,kppJan13_mld,kppFeb13_mld,kppMar13_mld,
    kppApr13_mld,kppMay13_mld,kppJun13_mld,kppJul13_mld,kppAug13_mld]
keps_MLD = np.c_[kepsSep12_mld,kepsOct12_mld,kepsNov12_mld,kepsDec12_mld,kepsJan13_mld,kepsFeb13_mld,kepsMar13_mld,
    kepsApr13_mld,kepsMay13_mld,kepsJun13_mld,kepsJul13_mld,kepsAug13_mld]

KPP_MLD = np.array(KPP_MLD[0])
keps_MLD = np.array(keps_MLD[0])

#############################################################################################
### ----------------------------- Ocean Reanalysis System 5 ----------------------------- ###
#############################################################################################

ORA5_temp = xr.open_dataset('votemper_combined.nc', decode_times=False);
ORA5_salt = xr.open_dataset('vosaline_combined.nc', decode_times=False);

## combine netcdf4 files with same predixes
# ORA5_salt = xr.open_mfdataset('vosaline_control_monthly_highres_3D_*.nc')
# ORA5_salt = ORA5_salt.to_netcdf('vosaline_combined.nc')

# ORA_t_Depth = ORA5_temp.deptht; ORA_t_Temp2 = ORA5_temp.votemper; ORA_t_time2 = ORA5_temp.time_counter
# ORA_t_lat = ORA5_temp.nav_lat; ORA_t_lon = ORA5_temp.nav_lon;

# print('ORA5 time', ORAtime2.values)
# print('ORA5 dataset', ORA5_temp.variables.items())
print('ORA5 dataset', ORA5_salt.deptht.values)
# print('ORAtemp3D items', ORA5_temp.variables.items())
# print('ORAtemp3D latitude', ORAlat[:,0:10], ORAlat.shape)
# print('ORAtemp3D longitude', ORAlon[:,0:10].values, ORAlon.shape)
# print('latitude', np.min(ORA_t_lat), np.max(ORA_t_lat))
# print('longitude', np.min(ORA_t_lon), np.max(ORA_t_lon)
#------------------------------------------------#
##### extract data for selected co-ordinates #####
#------------------------------------------------#
#TEMPERATURE
def extract(ORA5_temp, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_temp['nav_lon'].values.ravel(), ORA5_temp['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_temp['nav_lon'].shape))
    return ORA5_temp.isel(x=j0, y=i0).squeeze()
selectORATemp = extract(ORA5_temp, -16.25, 48.75);# print('extracts',selectORATemp)

#SALINITY
def extract(ORA5_salt, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_salt['nav_lon'].values.ravel(), ORA5_salt['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_salt['nav_lon'].shape))
    return ORA5_salt.isel(x=j0, y=i0).squeeze()
selectORASalt = extract(ORA5_salt, -16.25, 48.75);# print('extracts',selectORASalt)


# selectORATemp = selectORATemp.sel(time_counter=334)
# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORASalt = selectORASalt.sel(time_counter=334)
# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORATemp = selectORATemp.votemper

### OLD MLD DIAGNOSIS ###
mld_ORATemp = selectORATemp.sel(deptht=11.0, method='nearest') - 0.2;
z_indexes = abs(selectORATemp-mld_ORATemp).argmin('deptht'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
mldORA = selectORATemp.deptht[z_indexes]; 
# mldORA = mldORA.isel(nav_lat=0,nav_lon=0);  
print('mixed layer depth :',mldORA, mldORA.shape)

###########################################################################
### ----------------------------- FIGURES ----------------------------- ###
###########################################################################

month1 = np.array([datetime(2012,9,1),datetime(2012,10,1),datetime(2012,11,1),datetime(2012,12,1),
    datetime(2013,1,1),datetime(2013,2,1),datetime(2013,3,1),datetime(2013,4,1),datetime(2013,5,1),
    datetime(2013,6,1),datetime(2013,7,1),datetime(2013,8,1)])
month2 = np.array([datetime(2012,9,1),datetime(2012,10,1),datetime(2012,11,1),datetime(2012,12,1),
    datetime(2013,1,1),datetime(2013,2,1),datetime(2013,3,1),datetime(2013,4,1),datetime(2013,5,1),
    datetime(2013,6,1),datetime(2013,7,1),datetime(2013,8,1),datetime(2013,9,1)])

fig,ax = plt.subplots(1,1)

ax.plot(month1,np.negative(gliderOS_MLD),color='blue',label=r'OSMOSIS (48.7 $\degree$N, 16.2 $\degree$W)')
ax.plot(month2,np.negative(mldORA),color='orangered',label=r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
ax.plot(month1,KPP_MLD,color='burgundy',label='GOTM KPP scheme')
ax.plot(month1,keps_MLD,color='green',label=r'GOTM k-$\epsilon$ scheme')
ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
# ax.xaxis.set_major_Kcator(mdates.DateLocator())
ax.set_xlim(datetime(2012,9,1),datetime(2013, 8, 1))
plt.ylabel(r'Mixed Layer Depth -$m$'); plt.xlabel(r'Time')
plt.grid()
plt.show()