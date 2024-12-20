# -*- coding: utf-8 -*-
"""
Created on Tue March 15 10:46:56 2023
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
import netCDF4
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

# owsPAPA = Dataset('ows_papa.nc')
PAPAtemp = TempPAPAtemp.T_20; PAPAdepth = np.negative(TempPAPAtemp.depth);  # access a variable in the file 
# SPAPAtemp = STempPAPAtemp.T_20; SPAPAdepth = STempPAPAtemp.depth;  # access a variable in the file 
# print('depth', PAPAdepth)

# # convert 3D to 2D
# PAPAtemp = PAPAtemp.isel(lat=0,lon=0)[:,0]
# # print('temperature', PAPAtemp)


# time array
PAPAtime = TempPAPAtemp.time; #2012-09-04 00:00:00 8856
PAPAtime2 = 248 + (PAPAtime[:]/1440); # 2012-09-04 00:00:00 -- 2013-09-07 00:00:00
# PAPAtime = TempPAPAtemp.time; #2012-09-01 00:00:00

# print('time', PAPAtime[0:100])

# Create the starting date as a `datetime` object.
start = dt(2012, 9, 4, 0, 0, 0)
# List initialiser.
result = [start]

# Build a list of datetime objects for each hour of the year.
# for i in range(len(PAPAtime)):
for i in range(1, 8856):
    start += td(seconds=3600)
    result.append(start)

# Initialise a DataFrame data structure.
df = pd.DataFrame({'dates': result})
# Add integer days and time
df['time_delta'] = df['dates'] - pd.datetime(2012,1,1,0,0,0)
# Add floating days 
df['total_days_td'] = df['time_delta'].dt.total_seconds() / (24 * 60 * 60)
# Remove the datetime object columns.
df.drop(['dates'], inplace=True, axis=1); df.drop(['time_delta'], inplace=True, axis=1) # drop the columns
# create array and convert to 1D
df = np.array(df); df = df.flatten(); # print(df, df.shape)
# # plot
# plt.plot(df,PAPAtemp)
# plt.grid()
# plt.show()
# ##############
# ## PROFILES ##
# ##############

#  OWSP_sel = [0,10,20];# [0,13515]
#  # TEMPERATURE
#  for i in OWSP_sel:
#      PAPAtmp = PAPAtemp.isel(lon=0,lat=0)[i,:]; PAPAdpth = PAPAdepth[:]
#      PAPAtme = PAPAtime2[i]
#      plt.plot(PAPAtmp,PAPAdpth,label = r'OWS Papa %0.2f day'%PAPAtme.values);
#      plt.legend(loc='best', fontsize='small')
#      plt.ylabel('Depth -m');plt.xlabel(r'Temperature -\deg C'); 
#      plt.ylim(-200,10)
#      plt.grid(); 


# PAPAtemp = PAPAtemp.isel(lat=0,lon=0)
# # plt.plot(PAPAtemp,PAPAdepth.T)

# plt.plot(PAPAtemp[0,:],PAPAdepth,label = r'OWS Papa %0.2f day'%PAPAtime2[0].values)
# plt.plot(PAPAtemp[90,:],PAPAdepth,label = r'OWS Papa %0.2f day'%PAPAtime2[90].values)
# plt.plot(PAPAtemp[900,:],PAPAdepth,label = r'OWS Papa %0.2f day'%PAPAtime2[900].values)
# plt.plot(PAPAtemp[4000,:],PAPAdepth,label = r'OWS Papa %0.2f day'%PAPAtime2[4000].values)
# plt.plot(PAPAtemp[9460,:],PAPAdepth,label = r'OWS Papa %0.2f day'%PAPAtime2[9640].values)
# plt.plot(PAPAtemp[50000,:],PAPAdepth,label = r'OWS Papa %0.2f day'%PAPAtime2[50000].values)
# plt.ylabel('Depth -m');plt.xlabel(r'Temperature - $\degree$ C'); 
# plt.legend(loc='best', fontsize='small')
# plt.grid()
# plt.show()


##############################
### SECTION GOTM VARIABLES ###
##############################

### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
# gotmkppKm = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_2Km.nc", decode_times=False)
# gotmkppKh = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_2Kh.nc", decode_times=False)
# ## LaT addition ##
gotmkppLaRib = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat.nc",decode_times=False)
# ## h/L correction ##
# gotmkppLaRibhLPveP100 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP100.nc",decode_times=False)
gotmkppLaRibhLPveP1500 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP1500.nc",decode_times=False)
# gotmkppLaRibhLPveP10000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP10000.nc",decode_times=False)


#######################
## collect variables ##

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMsal = gotmkpp.salt; GTMo2 = gotmkpp.o2_obs; GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMxflx = gotmkpp.gamu; GTMyflx = gotmkpp.gamv; time = gotmkpp.time
# GTMtempKm = gotmkppKm.temp; GTMdepthKm = gotmkppKm.z; GTMnumKm = gotmkppKm.num; GTMnuhKm = gotmkppKm.nuh; GTMsalKm = gotmkppKm.salt; GTMo2Km = gotmkppKm.o2_obs; GTMuKm = gotmkppKm.u; GTMvKm = gotmkppKm.v; GTMxflxKm = gotmkppKm.gamu; GTMyflxKm = gotmkppKm.gamv; timeKm = gotmkppKm.time
# GTMtempKh = gotmkppKh.temp; GTMdepthKh = gotmkppKh.z; GTMnumKh = gotmkppKh.num; GTMnuhKh = gotmkppKh.nuh; GTMsalKh = gotmkppKh.salt; GTMo2Kh = gotmkppKh.o2_obs; GTMuKh = gotmkppKh.u; GTMvKh = gotmkppKh.v; GTMxflxKh = gotmkppKh.gamu; GTMyflxKh = gotmkppKh.gamv; timeKh = gotmkppKh.time
GTMtempLaRib = gotmkppLaRib.temp; GTMdepthLaRib = gotmkppLaRib.z; GTMnumLaRib = gotmkppLaRib.num; GTMnuhLaRib = gotmkppLaRib.nuh; GTMuLaRib = gotmkppLaRib.u; GTMvLaRib = gotmkppLaRib.v; timeLaRib = gotmkppLaRib.time

# GTMtempLaRibhLPveP100 = gotmkppLaRibhLPveP100.temp; GTMdepthLaRibhLPveP100 = gotmkppLaRibhLPveP100.z; GTMnumLaRibhLPveP100 = gotmkppLaRibhLPveP100.num; GTMnuhLaRibhLPveP100 = gotmkppLaRibhLPveP100.nuh; GTMuLaRibhLPveP100 = gotmkppLaRibhLPveP100.u; GTMvLaRibhLPveP100 = gotmkppLaRibhLPveP100.v; timeLaRibhLPveP100 = gotmkppLaRibhLPveP100.time
GTMtempLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.temp; GTMdepthLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.z; GTMnumLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.num; GTMnuhLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.nuh; GTMuLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u; GTMvLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.v; timeLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.time
# GTMtempLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.temp; GTMdepthLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.z; GTMnumLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.num; GTMnuhLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.nuh; GTMuLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u; GTMvLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.v; timeLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.time


##################
## convert time ##

time2 = 267.75 + (time[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
# time2Km = 267.75 + (timeKm[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
# time2Kh = 267.75 + (timeKh[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
time2LaRib = 267.75 + (timeLaRib[:]/86400); 

# time2LaRibhLPveP100 = 267.75 + (timeLaRibhLPveP100[:]/86400); 
time2LaRibhLPveP1500 = 267.75 + (timeLaRibhLPveP1500[:]/86400); 
# time2LaRibhLPveP10000 = 267.75 + (timeLaRibhLPveP10000[:]/86400); 

###################
# time allocation #

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
# gotmkppKm = gotmkppKm.assign(time2Km=("time", time2Km)); gotmkppKm = gotmkppKm.swap_dims({"time" : "time2Km"})
# gotmkppKh = gotmkppKh.assign(time2Kh=("time", time2Kh)); gotmkppKh = gotmkppKh.swap_dims({"time" : "time2Kh"})
gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})

# gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.assign(time2LaRibhLPveP100=("time", time2LaRibhLPveP100)); gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.swap_dims({"time" : "time2LaRibhLPveP100"})
gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.assign(time2LaRibhLPveP1500=("time", time2LaRibhLPveP1500)); gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.swap_dims({"time" : "time2LaRibhLPveP1500"})
# gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.assign(time2LaRibhLPveP10000=("time", time2LaRibhLPveP10000)); gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.swap_dims({"time" : "time2LaRibhLPveP10000"})


# ######################
# ## SEA SURFACE PLTS ##
# ######################
# ## convert 3D to 2D
# GTMtemp = GTMtemp.isel(lat=0,lon=0)
# ## sea surface temperature
# GTMtemp = GTMtemp.sel(z=-1, method='nearest')

# # Sea surface temperature plots
# fig,ax = plt.subplots(1,1,sharex=True)

# ax.plot(df,PAPAtemp,label='PMEL PAPA')
# ax.plot(time2,GTMtemp,label='GOTM KPP')
# ax.grid(); ax.legend(loc='lower right')
# ax.tick_params(axis='both', which='major', labelsize=12)


# fig.text(0.5, 0.01, r'Time -Julian days', ha='center', va='center', fontsize=12)
# fig.text(0.008, 0.5, r'Mixed Layer Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')


# plt.show()



##########################################
### SECTION GOTM and PAPA DIAGNOSE MLD ###
##########################################

#  ----------------------------- calculate MIXED LAYER DEPTH -----------------------------  #
# # new diagnosis to avoid the argmin() component which results in MLD collapses

## Psuedocode
# 1. loop through all levels of depth in the temperature variables from bottom to the top
# 2. flag where the sign changes from negative to positive
# 3. establish the position of the element to be MLD 

MLDval = []; 
MLDindx = []; MLDKmindx = []; MLDKhindx = []; MLDKPPLaRibindx = []; 
MLDKPPLaRibhLPveP100indx = []; MLDKPPLaRibhLPveP1500indx = []; MLDKPPLaRibhLPveP10000indx = [];
########################
## use orig KPP model ##
mld_tempKPP = GTMtemp.sel(z=-11.0, method='nearest') - 0.2;
mld_tempKPP = (GTMtemp-mld_tempKPP).isel(lat=0,lon=0).values
# print('the difference',(mld_tempKPP)[0:50,469:499])
# print('the difference',(mld_tempKPP)[0:50,469:499].shape)

for time in mld_tempKPP:
    for depth in range(len(time)):
        if time[depth] > 0:
            MLDindx.append(depth)
            break
MLDindx = np.asarray(MLDindx); mldKPP = GTMtemp.z[MLDindx];
# print('index for MLD :', MLDindx, MLDindx.shape)
# print('mldKPP:', mldKPP, mldKPP.shape)

##############################
## use Lat and Rib original ##
mld_tempKPPLaRib = GTMtempLaRib.sel(z=-11.0, method='nearest') - 0.2; 
mld_tempKPPLaRib = (GTMtempLaRib-mld_tempKPPLaRib).isel(lat=0,lon=0).values

for time in mld_tempKPPLaRib:
    for depth in range(len(time)):
        if time[depth] > 0:
            MLDKPPLaRibindx.append(depth)
            break
MLDKPPLaRibindx = np.asarray(MLDKPPLaRibindx); mldKPPLaRib = GTMtempLaRib.z[MLDKPPLaRibindx];


# ##################################################################
# ## use Lat and Rib original with hL condition (Ch=1 for hL>100) ##
# mld_tempKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
# mld_tempKPPLaRibhLPveP100 = (GTMtempLaRibhLPveP100-mld_tempKPPLaRibhLPveP100).isel(lat=0,lon=0).values

# for time in mld_tempKPPLaRibhLPveP100:
#     for depth in range(len(time)):
#         if time[depth] > 0:
#             MLDKPPLaRibhLPveP100indx.append(depth)
#             break
# MLDKPPLaRibhLPveP100indx = np.asarray(MLDKPPLaRibhLPveP100indx); 
# mldKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.z[MLDKPPLaRibhLPveP100indx];


#############################################################################
## use Lat and Rib original with hL condition (Ch=1 for hL>1500) for gamma ##
mld_tempKPPLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.sel(z=-11.0, method='nearest') - 0.2;
mld_tempKPPLaRibhLPveP1500 = (GTMtempLaRibhLPveP1500-mld_tempKPPLaRibhLPveP1500).isel(lat=0,lon=0).values

for time in mld_tempKPPLaRibhLPveP1500:
    for depth in range(len(time)):
        if time[depth] > 0:
            MLDKPPLaRibhLPveP1500indx.append(depth)
            break
MLDKPPLaRibhLPveP1500indx = np.asarray(MLDKPPLaRibhLPveP1500indx); 
mldKPPLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.z[MLDKPPLaRibhLPveP1500indx];


# ##############################################################################################
# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>10000) for gamma ##
# mld_tempKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.sel(z=-11.0, method='nearest') - 0.2; 
# mld_tempKPPLaRibhLPveP10000 = (GTMtempLaRibhLPveP10000-mld_tempKPPLaRibhLPveP10000).isel(lat=0,lon=0).values

# for time in mld_tempKPPLaRibhLPveP10000:
#     for depth in range(len(time)):
#         if time[depth] > 0:
#             MLDKPPLaRibhLPveP10000indx.append(depth)
#             break
# MLDKPPLaRibhLPveP10000indx = np.asarray(MLDKPPLaRibhLPveP10000indx); 
# mldKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.z[MLDKPPLaRibhLPveP10000indx];


# ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# ## -------------- PAPA --------------  ##
# ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# # print('temp dataset',TempPAPAtemp); 

# print('temp',PAPAtemp)

# profPAPAtemp = PAPAtemp.sel(time=10,lon=0,lat=0).values
#  ##PROFILE

# plt.plot(profPAPAtemp, PAPAdepth.values)
# plt.ylabel('Mixed layer depth -m')
# plt.xlabel('Time from 04-Sept-2012 until 30-Sep-2013')
# plt.title('OWS Papa Mixed layer depth evolution')
# plt.show()


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

# exit()
# PAPAtemp_imp = []

# for i, ele in enumerate(PAPAtemp):
     
#     # replace if greater than K
#     if ele > 100 :
#         previous = PAPAtemp[i-1]
#         PAPAtemp_imp.append(previous)
#     else :
#         PAPAtemp_imp.append(ele)

# print('improved dataset of OWS Papa temp :', PAPAtemp[9640,0:10].values)

#for values in temp:
    #where it takes 1.e+35f, caluclate the average between the value time[-1] and time[+1]
    # update the temp array with these values

PAPAtemp = pd.DataFrame(TempPAPAtemp.T_20)
PAPAtemp = PAPAtemp.fillna(PAPAtemp.mean())

mld_tempPAPA = PAPAtemp.sel(depth=10.0, method='nearest') - 0.2; 
z_indexes = abs(PAPAtemp-mld_tempPAPA).argmin('depth'); # print(z_indexes)
mldPAPA = PAPAtemp.depth[z_indexes]; 
mldPAPA = mldPAPA.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# print('index for MLD :', MLDindx, MLDindx.shape)
# print('temp:', PAPAtemp[9640,:].values)#, mldPAPA.shape)
# print('PAPA time:', PAPAtime2[53000:].values, PAPAtime.shape)#, PAPAtime.shape)
# exit()
# print('depth:', PAPAdepth[53000:].values)#, mldPAPA.shape)
# print('mldKPP:', mldPAPA[0:10].values, mldPAPA.shape)#, mldPAPA.shape)

print('temp:', PAPAtemp[9640,0:10].values)#, mldPAPA.shape)
print(np.where(PAPAtemp[9640,0:10]  == 1.00000004e+35))

# print(np.where(PAPAtemp >1000))


plt.plot(PAPAtime2, mldPAPA)
plt.ylabel('Mixed layer depth -m')
plt.xlabel('Time from 04-Sept-2012 until 30-Sep-2013')
plt.title('OWS Papa Mixed layer depth evolution')
plt.show()

exit()


fig, ax = plt.subplots(1,1,sharex=True)

# ax=ax.twiny(); 
# ax.set_axisbelow(True)
# # ax1.set_xlim(mdate[0],mdate[-1])
# xa = ax0.get_xaxis(); 
# xa.axis_date(); 
# xa.set_minor_locator(mdates.DayLocator(interval=1)); 
# xa.set_major_locator(mdates.MonthLocator()); 
# xa.set_major_formatter(mdates.DateFormatter('%b-%y'))
# ax.plot(yearmonth,fric_u,color='white',linestyle='dotted')


ax.plot(PAPAtime.values,mldPAPA, color='indigo',linestyle='dashed',label = 'PAPA')
# ax.plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
# ax.plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'_nolegend_')
# ax.plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta = 1$; $(h/L_L)_{MAX}=1500$')
ax.grid(); ax.legend(loc='lower right')
# ax[1].set_xlim(265,635);# ax[1].set_ylim(-200,0)
ax.set_xlim(480,550); ax.set_ylim(-150,0)
ax.tick_params(axis='both', which='major', labelsize=12)


# plt.legend(loc='lower right',fontsize='small')

fig.text(0.5, 0.01, r'Time -Julian days', ha='center', va='center', fontsize=12)
fig.text(0.008, 0.5, r'Mixed Layer Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')
# # # ax.text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# plt.suptitle(r'Comparison of mixed layer depths with $La_T$ and $h/L$ thresholds')
plt.show()

exit()



fig, ax = plt.subplots(1,1,sharex=True)

# ax=ax.twiny(); 
# ax.set_axisbelow(True)
# # ax1.set_xlim(mdate[0],mdate[-1])
# xa = ax0.get_xaxis(); 
# xa.axis_date(); 
# xa.set_minor_locator(mdates.DayLocator(interval=1)); 
# xa.set_major_locator(mdates.MonthLocator()); 
# xa.set_major_formatter(mdates.DateFormatter('%b-%y'))
# ax.plot(yearmonth,fric_u,color='white',linestyle='dotted')


ax.plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = '_nolegend_')
ax.plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
ax.plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'_nolegend_')
ax.plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta = 1$; $(h/L_L)_{MAX}=1500$')
ax.grid(); ax.legend(loc='lower right')
# ax[1].set_xlim(265,635);# ax[1].set_ylim(-200,0)
ax.set_xlim(480,550); ax[1].set_ylim(-150,0)
ax.tick_params(axis='both', which='major', labelsize=12)


# plt.legend(loc='lower right',fontsize='small')

fig.text(0.5, 0.01, r'Time -Julian days', ha='center', va='center', fontsize=12)
fig.text(0.008, 0.5, r'Mixed Layer Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')
# # # ax.text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
ax[0].text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.02, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[2].text(0.02, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# plt.suptitle(r'Comparison of mixed layer depths with $La_T$ and $h/L$ thresholds')
plt.show()










######################
exit()
######################


plt.plot(PAPAtemp,PAPAdepth,label='PAPA')
plt.grid()
plt.show()


exit()

##!! The temperature is set up in such a way that it si set to zero in GOTM model
#temperature = data.variables['temp']
#Temp = temperature[:]

## Extract the depth data
Depth = data.variables['z']
Z = Depth[:]
iDepth = data.variables['zi']
Zi = iDepth[:]

#data0 = cdms.open('wave_breaking.nc')
# print(data.variables.keys())
