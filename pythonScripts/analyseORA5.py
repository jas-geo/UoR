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
# ### ------------------------------- Argo Floats ----------------------------- ##
# ARGO_1 = xr.open_dataset('../GL_PR_PF_6901228.nc', decode_times=False);

# # print(ARGO_1.variables.items())

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

#-------------------------------#
##### averaged monthly data #####
#-------------------------------#

## September 2012
gliderSep12 = np.where((new_time > 244) & (new_time <= 274)); Sep12_time = new_time[gliderSep12]; 
Sep12_temp = GldrTemp[0:276,:]; Sep12_temp = Sep12_temp.mean(axis=0)
Sep12_sal = GldrSalt[0:276,:]; Sep12_sal = Sep12_sal.mean(axis=0)
## October 2012
gliderOct12 = np.where((new_time > 274) & (new_time <= 305)); Oct12_time = new_time[gliderOct12]; 
Oct12_temp = GldrTemp[277:588,:]; Oct12_temp = Oct12_temp.mean(axis=0)
Oct12_sal = GldrSalt[277:588,:]; Oct12_sal = Oct12_sal.mean(axis=0)
## November 2012
gliderNov12 = np.where((new_time > 305) & (new_time <= 335)); Nov12_time = new_time[gliderNov12]; 
Nov12_temp = GldrTemp[589:894,:]; Nov12_temp = Nov12_temp.mean(axis=0)
Nov12_sal = GldrSalt[589:894,:]; Nov12_sal = Nov12_sal.mean(axis=0)
## December 2012
gliderDec12 = np.where((new_time > 335) & (new_time <= 366)); Dec12_time = new_time[gliderDec12]; 
Dec12_temp = GldrTemp[895:1217,:]; Dec12_temp = Dec12_temp.mean(axis=0)
Dec12_sal = GldrSalt[895:1217,:]; Dec12_sal = Dec12_sal.mean(axis=0)
## January 2013
gliderJan13 = np.where((new_time > 366) & (new_time <= 397)); Jan13_time = new_time[gliderJan13]; 
Jan13_temp = GldrTemp[1218:1552,:]; Jan13_temp = Jan13_temp.mean(axis=0)
Jan13_sal = GldrSalt[1218:1552,:]; Jan13_sal = Jan13_sal.mean(axis=0)
## February 2013
gliderFeb13 = np.where((new_time > 397) & (new_time <= 456)); Feb13_time = new_time[gliderFeb13]; 
Feb13_temp = GldrTemp[1553:2248,:]; Feb13_temp = Feb13_temp.mean(axis=0)
Feb13_sal = GldrSalt[1553:2248,:]; Feb13_sal = Feb13_sal.mean(axis=0)
## March 2013
gliderMar13 = np.where((new_time > 456) & (new_time <= 487)); Mar13_time = new_time[gliderMar13]; 
Mar13_temp = GldrTemp[2249:2632,:]; Mar13_temp = Mar13_temp.mean(axis=0)
Mar13_sal = GldrSalt[2249:2632,:]; Mar13_sal = Mar13_sal.mean(axis=0)
## April 2013
gliderApr13 = np.where((new_time > 487) & (new_time <= 517)); Apr13_time = new_time[gliderApr13]; 
Apr13_temp = GldrTemp[2633:2992,:]; Apr13_temp = Apr13_temp.mean(axis=0)
Apr13_sal = GldrSalt[2633:2992,:]; Apr13_sal = Apr13_sal.mean(axis=0)
## May 2013
gliderMay13 = np.where((new_time > 517) & (new_time <= 548)); May13_time = new_time[gliderMay13]; 
May13_temp = GldrTemp[2993:3356,:]; May13_temp = May13_temp.mean(axis=0)
May13_sal = GldrSalt[2993:3356,:]; May13_sal = May13_sal.mean(axis=0)
## June 2013
gliderJun13 = np.where((new_time > 548) & (new_time <= 578)); Jun13_time = new_time[gliderJun13]; 
Jun13_temp = GldrTemp[3357:3682,:]; Jun13_temp = Jun13_temp.mean(axis=0)
Jun13_sal = GldrSalt[3357:3682,:]; Jun13_sal = Jun13_sal.mean(axis=0)
## July 2013
gliderJul13 = np.where((new_time > 578) & (new_time <= 609)); Jul13_time = new_time[gliderJul13]; 
Jul13_temp = GldrTemp[3683:4000,:]; Jul13_temp = Jul13_temp.mean(axis=0)
Jul13_sal = GldrSalt[3683:4000,:]; Jul13_sal = Jul13_sal.mean(axis=0)
## August 2013
gliderAug13 = np.where((new_time > 609) & (new_time <= 640)); Aug13_time = new_time[gliderAug13]; 
Aug13_temp = GldrTemp[4001:4095,:]; Aug13_temp = Aug13_temp.mean(axis=0)
Aug13_sal = GldrSalt[4001:4095,:]; Aug13_sal = Aug13_sal.mean(axis=0)

#############################################################################################
### ----------------------------- Ocean Reanalysis System 5 ----------------------------- ###
#############################################################################################

ORA5_temp = xr.open_dataset('votemper_combined.nc', decode_times=False);
ORA5_salt = xr.open_dataset('vosaline_combined.nc', decode_times=False);
ORA5_wtrflx = xr.open_dataset('sowaflup_combined.nc', decode_times=False)

# # combine netcdf4 files with same predixes
# ORA5_wtrflx = xr.open_mfdataset('sowaflup_control_monthly_highres_2D_*.nc')
# ORA5_wtrflx = ORA5_wtrflx.to_netcdf('sowaflup_combined.nc')

ORA_t_Depth = ORA5_temp.deptht; ORA_t_Temp2 = ORA5_temp.votemper; ORA_t_time2 = ORA5_temp.time_counter
ORA_t_lat = ORA5_temp.nav_lat; ORA_t_lon = ORA5_temp.nav_lon;

# print('ORA5 time', ORAtime2.values)
# print('ORA5 dataset', ORA5_temp.variables.items())
# print('ORA5 dataset', ORA5_wtrflx.variables.items())
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
selectORATempO = extract(ORA5_temp, -16.25, 48.75);# print('extracts',selectORATemp)
selectORATemp1 = extract(ORA5_temp, -16.25, 48.50);#North
selectORATemp2 = extract(ORA5_temp, -16.00, 48.75);#East
selectORATemp3 = extract(ORA5_temp, -16.00, 48.50);#West
# selectORATemp4 = extract(ORA5_temp, -16.25, 46.75);#South


#SALINITY
def extract(ORA5_salt, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_salt['nav_lon'].values.ravel(), ORA5_salt['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_salt['nav_lon'].shape))
    return ORA5_salt.isel(x=j0, y=i0).squeeze()
selectORASaltO = extract(ORA5_salt, -16.25, 48.75);# print('extracts',selectORASalt)
selectORASalt1 = extract(ORA5_salt, -16.25, 48.50);#North
selectORASalt2 = extract(ORA5_salt, -16.00, 48.75);#East
selectORASalt3 = extract(ORA5_salt, -16.00, 48.50);#West
# selectORASalt4 = extract(ORA5_salt, -16.25, 46.75);#South


#WATER FLUX
def extract(ORA5_wtrflx, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_wtrflx['nav_lon'].values.ravel(), ORA5_wtrflx['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_wtrflx['nav_lon'].shape))
    return ORA5_wtrflx.isel(x=j0, y=i0).squeeze()
selectORAwFluxO = extract(ORA5_wtrflx, -16.25, 48.75);# print('extracts',selectORASalt)
selectORAwFlux1 = extract(ORA5_wtrflx, -16.25, 48.5);#North
selectORAwFlux2 = extract(ORA5_wtrflx, -16.0, 48.75);#East
selectORAwFlux3 = extract(ORA5_wtrflx, -16.0, 48.5);#West
# selectORAwFlux4 = extract(ORA5_wtrflx, -16.25, 46.75);#South

# print(selectORATempO.info())


# ## extract depth information for specific time counters
# selectORATempO = selectORATempO.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORATemp1 = selectORATemp1.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORATemp2 = selectORATemp2.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORATemp3 = selectORATemp3.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORATemp4 = selectORATemp4.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)


# selectORASaltO = selectORASaltO.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORASalt1 = selectORASalt1.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORASalt2 = selectORASalt2.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORASalt3 = selectORASalt3.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
# selectORASalt4 = selectORASalt4.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)

#--------------------------#
##### plotting figures #####
#--------------------------#
# #TEMPERATURE
# plt.plot(selectORATempO.votemper,np.negative(selectORATempO.deptht), label = r'ORA5 (48.7 $\degree$N, 16.2 $\degree$W) 2013-04')
# plt.plot(selectORATemp1.votemper,np.negative(selectORATemp1.deptht), linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W) 2013-04')
# plt.plot(selectORATemp2.votemper,np.negative(selectORATemp2.deptht), linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W) 2013-04')
# plt.plot(selectORATemp3.votemper,np.negative(selectORATemp3.deptht), linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W) 2013-04')
# plt.plot(selectORATemp4.votemper,np.negative(selectORATemp4.deptht), linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W) 2013-04')
# plt.plot(Apr13_temp,np.negative(glider.pres_grid), label = 'O5MOSIS 2013-04')
# plt.legend(loc='lower right')
# plt.ylabel('Depth -m'); plt.xlabel(r'Temperature - $\degree$ C')
# plt.ylim(-1000,0); plt.xlim(7,19)
# plt.grid()
# plt.show()

# #SALINITY
# plt.plot(selectORASaltO.vosaline,np.negative(selectORASaltO.deptht), label = r'ORA5 (48.7 $\degree$N, 16.2 $\degree$W) 2013-04')
# plt.plot(selectORASalt1.vosaline,np.negative(selectORASalt1.deptht), linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W) 2013-04')
# plt.plot(selectORASalt2.vosaline,np.negative(selectORASalt2.deptht), linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W) 2013-04')
# plt.plot(selectORASalt3.vosaline,np.negative(selectORASalt3.deptht), linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W) 2013-04')
# plt.plot(selectORASalt4.vosaline,np.negative(selectORASalt4.deptht), linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W) 2013-04')
# plt.plot(Apr13_sal,np.negative(glider.pres_grid), label = 'O5MOSIS 2013-04')
# plt.legend(loc='lower right')
# plt.ylabel('Depth -m'); plt.xlabel(r'Salinity - PSU')
# plt.ylim(-1000,0); plt.xlim(35.1,35.8)
# plt.grid()
# plt.show()
month = np.array([datetime(2012,9,1),datetime(2012,10,1),datetime(2012,11,1),datetime(2012,12,1),
    datetime(2013,1,1),datetime(2013,2,1),datetime(2013,3,1),datetime(2013,4,1),datetime(2013,5,1),
    datetime(2013,6,1),datetime(2013,7,1),datetime(2013,8,1)])

# #WATERFLUX

# fig,ax = plt.subplots(1,1)
# # plt.plot(selectORAwFluxO.time_counter,selectORAwFluxO.sowaflup, label = r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax.plot(month,selectORAwFluxO.sowaflup, label = r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax.plot(month,selectORAwFlux1.sowaflup, linestyle = 'dashed', label = r'ORA5 (48.5 $\degree$N, 16.25 $\degree$W)')
# ax.plot(month,selectORAwFlux2.sowaflup, linestyle = 'dotted', label = r'ORA5 (48.75 $\degree$N, 16.00 $\degree$W)')
# ax.plot(month,selectORAwFlux3.sowaflup, linestyle = 'dashed', label = r'ORA5 (48.50 $\degree$N, 16.00 $\degree$W)')
# ax.legend(loc='lower right')
# # ax.xaxis.set_major_locator(mdates.YearLocator())
# # ax.xaxis.set_major_locator(mdates.DateLocator())
# ax.set_xlim(datetime(2012,9,1),datetime(2013, 8, 1))
# plt.ylabel(r'Water flux - $Kg \ m^{-2} \ s^{-1}$', fontsize=12);# plt.xlabel(r'Time -days since 2012-09-16 00:00:00', fontsize=12                                                                                      )
# # plt.ylim(-1000,0); plt.xlim(35.1,35.8)
# plt.grid()
# plt.show()



# #  ################################  # #
# 3-MONTHLY PANELS FOR VARIABLE PROFILES #
# #  ################################  # #
# !! disable lines 167-178 for plot the following figures #

#TEMPERATURE
## S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
fig, ax = plt.subplots(1,3,sharey=True,sharex=True)

ax[0].plot(selectORATempO.votemper.sel(time_counter=0), np.negative(selectORATempO.deptht), label=r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax[0].plot(selectORATemp1.votemper.sel(time_counter=91), np.negative(selectORATemp1.deptht), linestyle = 'dashed', label=r'ORA5 (48.5 $\degree$N, 16.25 $\degree$W)')
# ax[0].plot(selectORATemp2.votemper.sel(time_counter=91), np.negative(selectORATemp2.deptht), linestyle = 'dotted', label=r'ORA5 (48.75 $\degree$N, 16.0 $\degree$W)')
# ax[0].plot(selectORATemp3.votemper.sel(time_counter=91), np.negative(selectORATemp3.deptht), linestyle = 'dashed', label=r'ORA5 (48.5 $\degree$N, 16.0 $\degree$W)')
ax[0].plot(Sep12_temp, np.negative(glider.pres_grid), label=r'OSMOSIS')
ax[0].grid(); ax[0].legend(loc='center right',fontsize='small')
ax[0].set_xlim(10,19); ax[0].set_ylim(-500,0)
# ax[0].set_xlim(10,14); ax[0].set_ylim(-600,0)
ax[0].tick_params(axis='both', which='major', labelsize=12)


ax[1].plot(selectORATempO.votemper.sel(time_counter=30), np.negative(selectORATempO.deptht), label=r'_nolegend_')
# ax[1].plot(selectORATemp1.votemper.sel(time_counter=122), np.negative(selectORATemp1.deptht), linestyle = 'dashed', label=r'_nolegend_')
# ax[1].plot(selectORATemp2.votemper.sel(time_counter=122), np.negative(selectORATemp2.deptht), linestyle = 'dotted', label=r'_nolegend_')
# ax[1].plot(selectORATemp3.votemper.sel(time_counter=122), np.negative(selectORATemp3.deptht), linestyle = 'dashed', label=r'_nolegend_')
ax[1].plot(Oct12_temp, np.negative(glider.pres_grid), label=r'_nolegend_')
ax[1].grid();
ax[1].set_xlim(10,19); ax[1].set_ylim(-500,0)
# ax[1].set_xlim(10,14); ax[0].set_ylim(-600,0)
ax[1].tick_params(axis='both', which='major', labelsize=12)


ax[2].plot(selectORATempO.votemper.sel(time_counter=61), np.negative(selectORATempO.deptht), label=r'_nolegend_')
# ax[2].plot(selectORATemp1.votemper.sel(time_counter=152), np.negative(selectORATemp1.deptht), linestyle = 'dashed', label=r'_nolegend_')
# ax[2].plot(selectORATemp2.votemper.sel(time_counter=152), np.negative(selectORATemp2.deptht), linestyle = 'dotted', label=r'_nolegend_')
# ax[2].plot(selectORATemp3.votemper.sel(time_counter=152), np.negative(selectORATemp3.deptht), linestyle = 'dashed', label=r'_nolegend_')
ax[2].plot(Nov12_temp, np.negative(glider.pres_grid), label=r'_nolegend_')
ax[2].grid();
ax[2].set_xlim(10,19); ax[2].set_ylim(-500,0)
# ax[1].set_xlim(10,14); ax[0].set_ylim(-600,0)
ax[2].tick_params(axis='both', which='major', labelsize=12)


fig.text(0.5, 0.01, r'Temperature -$\degree C$', ha='center', va='center', fontsize=12)
fig.text(0.008, 0.5, r'Ocean Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')

ax[0].text(0.025, 0.05, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.9, 0.045, 'SEP', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[0].transAxes)

ax[1].text(0.025, 0.05, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.9, 0.045, 'OCT', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[1].transAxes)

ax[2].text(0.025, 0.05, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].text(0.9, 0.045, 'NOV', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[2].transAxes)

plt.show()

# SALINITY 
## S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
fig, ax = plt.subplots(1,3,sharey=True,sharex=True)

ax[0].plot(selectORASaltO.vosaline.sel(time_counter=0), np.negative(selectORASaltO.deptht), label=r'_nolegend_')
# ax[0].plot(selectORASalt1.vosaline.sel(time_counter=91), np.negative(selectORASalt1.deptht), linestyle = 'dashed', label=r'_nolegend_')
# ax[0].plot(selectORASalt2.vosaline.sel(time_counter=91), np.negative(selectORASalt2.deptht), linestyle = 'dotted', label=r'_nolegend_')
# ax[0].plot(selectORASalt3.vosaline.sel(time_counter=91), np.negative(selectORASalt3.deptht), linestyle = 'dashed', label=r'_nolegend_')
ax[0].plot(Sep12_sal, np.negative(glider.pres_grid), label=r'_nolegend_')
ax[0].grid();# ax[0].legend(loc='upper left') 
ax[0].set_xlim(35.1,35.8); ax[0].set_ylim(-500,0)
ax[0].tick_params(axis='both', which='major', labelsize=12)


ax[1].plot(selectORASaltO.vosaline.sel(time_counter=30), np.negative(selectORASaltO.deptht), label=r'_nolegend_')
# ax[1].plot(selectORASalt1.vosaline.sel(time_counter=122), np.negative(selectORASalt1.deptht), linestyle = 'dashed', label=r'_nolegend_')
# ax[1].plot(selectORASalt2.vosaline.sel(time_counter=122), np.negative(selectORASalt2.deptht), linestyle = 'dotted', label=r'_nolegend_')
# ax[1].plot(selectORASalt3.vosaline.sel(time_counter=122), np.negative(selectORASalt3.deptht), linestyle = 'dashed', label=r'_nolegend_')
ax[1].plot(Oct12_sal, np.negative(glider.pres_grid), label=r'_nolegend_')
ax[1].grid();
ax[1].set_xlim(35.1,35.8); ax[1].set_ylim(-500,0)
ax[1].tick_params(axis='both', which='major', labelsize=12)


ax[2].plot(selectORASaltO.vosaline.sel(time_counter=61), np.negative(selectORASaltO.deptht), label=r'ORA5 (48.75 $\degree$N, 16.25 $\degree$W)')
# ax[2].plot(selectORASalt1.vosaline.sel(time_counter=152), np.negative(selectORASalt1.deptht), linestyle = 'dashed', label=r'ORA5 (48.50 $\degree$N, 16.25 $\degree$W)')
# ax[2].plot(selectORASalt2.vosaline.sel(time_counter=152), np.negative(selectORASalt2.deptht), linestyle = 'dotted', label=r'ORA5 (48.75 $\degree$N, 16.00 $\degree$W)')
# ax[2].plot(selectORASalt3.vosaline.sel(time_counter=152), np.negative(selectORASalt3.deptht), linestyle = 'dashed', label=r'ORA5 (48.50 $\degree$N, 16.00 $\degree$W)')
ax[2].plot(Nov12_sal, np.negative(glider.pres_grid), label=r'OSMOSIS (48.7 $\degree$N, 16.2 $\degree$W)')
ax[2].grid(); ax[2].legend(loc='center left',fontsize='small') 
ax[2].set_xlim(35.1,35.8); ax[2].set_ylim(-500,0)
ax[2].tick_params(axis='both', which='major', labelsize=12)

fig.text(0.5, 0.01, r'Salinity -$PSU$', ha='center', va='center', fontsize=12)
fig.text(0.008, 0.5, r'Ocean Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')

ax[0].text(0.025, 0.05, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.9, 0.045, 'SEP', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[0].transAxes)

ax[1].text(0.025, 0.05, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.9, 0.045, 'OCT', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[1].transAxes)

ax[2].text(0.025, 0.05, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[2].text(0.9, 0.045, 'NOV', horizontalalignment='center', verticalalignment='center', fontsize = 12, bbox = dict(facecolor = 'wheat', alpha = 0.5), transform=ax[2].transAxes)

plt.show()


exit()

# plot options
selectORATemp.plot(selectORATemp.votemper,np.negative(selectORATemp.deptht), label = 'ORA5 time=%d'%selectORATemp.time_counter)


