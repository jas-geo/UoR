## ------------- Written code to extract initial conditions from Glider ------------- ##
import numpy as np
import math
import sys
import csv
import os
import xarray as xr
import pandas as pd
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format})

### ----------------------------- GLIDER ----------------------------- ##

gliderdata= Dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc'); 
glider = xr.open_dataset("../gotm-4.0.0-kpp/simulations/glider_timeseries.nc", decode_times=False)

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; Gldrdays = glider.dats; GldrSalt = glider.prac_sal;

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days

glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})

### ----------------------------- GOTM ----------------------------- ##
gotmkppLaRibhLPveP3000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_placebo1.nc",decode_times=False)
# gotmkppLaRibhLPveP3000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP3000.nc",decode_times=False)
# gotmkppLaRibhLPveP3000Annual = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP3000.nc",decode_times=False)

## variables
GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; 
GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMsalt = gotmkpp.salt; 

# GTMtempLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.temp; GTMdepthLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.z; 
# GTMnumLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.num; GTMnuhLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.nuh; 
# GTMuLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.u; GTMvLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.v; 
# GTMsaltLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.salt; 

# GTMtempLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.temp; GTMdepthLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.z; 
# GTMnumLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.num; GTMnuhLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.nuh; 
# GTMuLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.u; GTMvLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.v; 

## time
time = gotmkpp.time; time2 = 470.0 + (time[:]/86400); 
gotmkpp = gotmkpp.assign(time2=("time", time2)); 
gotmkpp = gotmkpp.swap_dims({"time" : "time2"})

# timeLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.time; time2LaRibhLPveP3000 = 470.0 + (timeLaRibhLPveP3000[:]/86400); 
# gotmkppLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.assign(time2LaRibhLPveP3000=("time", time2LaRibhLPveP3000)); 
# gotmkppLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.swap_dims({"time" : "time2LaRibhLPveP3000"})

# timeLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.time; time2LaRibhLPveP3000Annual = 267.75 + (timeLaRibhLPveP3000Annual[:]/86400); 
# gotmkppLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.assign(time2LaRibhLPveP3000Annual=("time", time2LaRibhLPveP3000Annual)); 
# gotmkppLaRibhLPveP3000Annual = gotmkppLaRibhLPveP3000Annual.swap_dims({"time" : "time2LaRibhLPveP3000Annual"})

# print(gotmkppLaRibhLPveP3000Annual.variables.items())

##########################################################################################
##  ---------------------------  obtain initial conditions  --------------------------- ##
##########################################################################################
# print('Glider depth', GldrDepth.values)
# t=464 : 2340 # t=470 : 2410 # t=480 : 2543 # 
print('Glider fixed time :', new_time[2605:])
exit()
## temperature ##
# print('Glider temp for time=464', GldrTemp[2340,:].values, len(GldrTemp[2340,:]))
Gldrtemp470 = GldrTemp[2410,:].values
t_profile470 = np.c_[GldrDepth,GldrTemp[2410,:].values]
print('Glider temp for time=470', t_profile470)

## salinity ##
# print('Glider salt for time=464', GldrSalt[2340,:].values, len(GldrSalt[2340,:]))
GldrSalt470 = GldrSalt[2410,:].values
s_profile470 = np.c_[GldrDepth,GldrSalt[2410,:].values]
print('Glider salt for time=470', s_profile470)

##########################################################################################
### --------------------------- GOTM Mixed Layer Depths --------------------------- ###
##########################################################################################

## TEMPERATURE ##
# GTMTemp_TimeP3000 = np.c_[time2LaRibhLPveP3000, GTMtempLaRibhLPveP3000]
# print('GTM time&temp h/L>3000 values', GTMtempLaRibhLPveP3000[0:10]) 
## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP3000 = GTMtempLaRibhLPveP3000.sel(z=-5.0, method='nearest') - 0.2; 
# print('initial MLD',mld_tempKPPLaRibhLPveP3000)
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
difftempMLD = abs(GTMtempLaRibhLPveP3000-mld_tempKPPLaRibhLPveP3000)
# print('diff in temp & MLD', difftempMLD)
z_indexes = difftempMLD.argmin('z'); 
# print('z indexes', z_indexes)
# print('z values', GTMtempLaRibhLPveP3000.z.values)
mldKPPLaRibhLPveP3000 = GTMtempLaRibhLPveP3000.z[z_indexes];
# print('MLD with lat and lon', mldKPPLaRibhLPveP3000)
mldKPPLaRibhLPveP3000 = mldKPPLaRibhLPveP3000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)
# print('GOTM final MLD h/L > 3000', mldKPPLaRibhLPveP3000)

## SALINTIY ##
## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_saltKPPLaRibhLPveP3000 = GTMsaltLaRibhLPveP3000.sel(z=-5.0, method='nearest') - 0.03; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMsaltLaRibhLPveP3000-mld_saltKPPLaRibhLPveP3000).argmin('z'); # print(z_indexes)
mldSKPPLaRibhLPveP3000 = GTMsaltLaRibhLPveP3000.z[z_indexes]; 
mldSKPPLaRibhLPveP3000 = mldSKPPLaRibhLPveP3000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)
# print('GOTM MLD h/L > 3000', mldKPPLaRibhLPveP3000Annual)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP3000Annual = GTMtempLaRibhLPveP3000Annual.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP3000Annual-mld_tempKPPLaRibhLPveP3000Annual).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP3000Annual = GTMtempLaRibhLPveP3000Annual.z[z_indexes]; 
mldKPPLaRibhLPveP3000Annual = mldKPPLaRibhLPveP3000Annual.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)
# print('GOTM MLD h/L > 3000', mldKPPLaRibhLPveP3000Annual)

###############################################################################################
### ----------------------------- mixed layer depth .csv file ----------------------------- ###
###############################################################################################
def read_lines():
    with open('mld.csv', 'rU') as mldFile:
        reader = csv.reader(mldFile)
        for row in reader:
            yield [ float(i) for i in row ]
for i in read_lines():
    print('mld :',i)

def read_lines():
    with open('dats.csv', 'rU') as datFile:
        reader_ = csv.reader(datFile)
        for row in reader_:
            yield [ float(j) for j in row ]

for j in read_lines():
    print('date :',j)

j = np.array(j)
j = j + 241

print(j)

i_j = np.c_[j,i]
print(i_j)


##########################################################################################
## ------------------------------- moving average MLDs--------------------------------- ##
##########################################################################################

def moving_average(i, n=3) :
    ret = np.cumsum(i, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

df_mld = moving_average(i)
df_mld5 = moving_average(i,n=5)
df_mld7 = moving_average(i,n=7)

df_dats = moving_average(j)
df_dats5 = moving_average(j,n=5)
df_dats7 = moving_average(j,n=7)

print('manualy moving average :', df_mld)

############################################################
## Comparison of Annual and 480-633MLDs using hL > [3000] ##
############################################################

fig, ax = plt.subplots(2,1)

# ax[0].plot(df_dats,np.negative(df_mld), color='rosybrown',linestyle='solid')#,label = 'OSMOSIS moving average n=5')
ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2LaRibhLPveP3000,mldKPPLaRibhLPveP3000,color='maroon',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(465,520);# ax[0].set_ylim(-400,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2LaRibhLPveP3000,mldSKPPLaRibhLPveP3000,color='maroon',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(465,520); ax[1].set_ylim(-400,0)

# plt.legend(loc='best',fontsize='small')

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs using hL thresholds w NON-LOCAL mixing [3000,3500,4000]')
plt.show()

exit()
##########################################################################################
### --------------------------- OSMOSIS Mixed Layer Depths --------------------------- ###
##########################################################################################
mld_tempOS = GldrTemp.sel(pressure=11) - 0.2;
print('Glider MLD temp at reference',mld_tempOS[2543:2600])
print('Glider temp',GldrTemp[2543:2600,0:15])
z_indexes = abs(GldrTemp[2543:2600,0:15]-mld_tempOS[2543:2600]).argmin('pressure'); 
print(z_indexes)
mldOS = GldrTemp.pressure[z_indexes]; 
# mldOS = mldOS.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)
print(mldOS)

exit()

# ##  ---------------------------  obtain initial conditions  --------------------------- ##
# # print('Glider depth', GldrDepth.values)
print('Glider fixed time :', new_time[2340:2434])
# # temperature
# # print('Glider temp for time=464', GldrTemp[2340,:].values, len(GldrTemp[2340,:]))
# Gldrtemp464 = GldrTemp[2340,:].values
# t_profile464 = np.c_[GldrDepth,GldrTemp[2340,:].values]
# print('Glider temp for time=464', t_profile464)
# # salinity
# # print('Glider salt for time=464', GldrSalt[2340,:].values, len(GldrSalt[2340,:]))
# GldrSalt464 = GldrSalt[2340,:].values
# s_profile464 = np.c_[GldrDepth,GldrSalt[2340,:].values]
# print('Glider salt for time=464', s_profile464)