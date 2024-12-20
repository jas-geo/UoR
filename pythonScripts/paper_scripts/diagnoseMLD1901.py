## Sensible heat flux and solar radiation ##
import numpy as np
import math
import sys
import csv
from scipy.stats.stats import pearsonr
from scipy import interpolate
import os
import xarray as xr
import pandas as pd
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
import matplotlib.dates as mdates
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, ScalarFormatter,LogFormatter)
from bisect import bisect_right

plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

####################################
### SECTION COLLECT OBSERVATIONS ###
####################################

### ----------------------------- GLIDER ----------------------------- ##
# file1 = 'C:/home/users/mc837749/Documents/gotm-4.0.0/simulations/annualOSMOSISforcing_2502/glider_timeseries.nc'
glider = xr.open_dataset('~/Documents/gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; GldrDens = glider.pot_den; 
Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen



## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days


new_time[2538]=479.780199
new_time[2539]=479.796973
new_time[2540]=479.813747
# print(new_time[2535:2544].values)

new_time2 = []

def julian_to_date():
	initial_date = "2012-01-01 00:00:00"
	global new_time2
	for i in new_time.values:
		date = pd.to_datetime(initial_date) + pd.DateOffset(days=i)
		new_time2.append(date)
		# print(date)
	return print('Time format is in Date-Time == SUCCESS')

julian_to_date()

# print(new_time2)

glider = glider.assign(new_time2=("time", new_time2))
glider.new_time2.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time2', inplace=True); glider = glider.swap_dims({"time" : "new_time2"})
##############################
### SECTION GOTM VARIABLES ###
##############################

### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
gotmkppKm = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_2Km.nc", decode_times=False)
gotmkppKh = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_2Kh.nc", decode_times=False)
## LaT addition ##
gotmkppLaRib = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat.nc",decode_times=False)
## h/L correction ##
gotmkppLaRibhLPveP100 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP100.nc",decode_times=False)
gotmkppLaRibhLPveP1500 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP1500.nc",decode_times=False)
gotmkppLaRibhLPveP10000 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP10000.nc",decode_times=False)


#######################
## collect variables ##

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMsal = gotmkpp.salt; GTMo2 = gotmkpp.o2_obs; GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMxflx = gotmkpp.gamu; GTMyflx = gotmkpp.gamv; time = gotmkpp.time
GTMtempKm = gotmkppKm.temp; GTMdepthKm = gotmkppKm.z; GTMnumKm = gotmkppKm.num; GTMnuhKm = gotmkppKm.nuh; GTMsalKm = gotmkppKm.salt; GTMo2Km = gotmkppKm.o2_obs; GTMuKm = gotmkppKm.u; GTMvKm = gotmkppKm.v; GTMxflxKm = gotmkppKm.gamu; GTMyflxKm = gotmkppKm.gamv; timeKm = gotmkppKm.time
GTMtempKh = gotmkppKh.temp; GTMdepthKh = gotmkppKh.z; GTMnumKh = gotmkppKh.num; GTMnuhKh = gotmkppKh.nuh; GTMsalKh = gotmkppKh.salt; GTMo2Kh = gotmkppKh.o2_obs; GTMuKh = gotmkppKh.u; GTMvKh = gotmkppKh.v; GTMxflxKh = gotmkppKh.gamu; GTMyflxKh = gotmkppKh.gamv; timeKh = gotmkppKh.time
GTMtempLaRib = gotmkppLaRib.temp; GTMdepthLaRib = gotmkppLaRib.z; GTMnumLaRib = gotmkppLaRib.num; GTMnuhLaRib = gotmkppLaRib.nuh; GTMuLaRib = gotmkppLaRib.u; GTMvLaRib = gotmkppLaRib.v; timeLaRib = gotmkppLaRib.time

GTMtempLaRibhLPveP100 = gotmkppLaRibhLPveP100.temp; GTMdepthLaRibhLPveP100 = gotmkppLaRibhLPveP100.z; GTMnumLaRibhLPveP100 = gotmkppLaRibhLPveP100.num; GTMnuhLaRibhLPveP100 = gotmkppLaRibhLPveP100.nuh; GTMuLaRibhLPveP100 = gotmkppLaRibhLPveP100.u; GTMvLaRibhLPveP100 = gotmkppLaRibhLPveP100.v; timeLaRibhLPveP100 = gotmkppLaRibhLPveP100.time
GTMtempLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.temp; GTMdepthLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.z; GTMnumLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.num; GTMnuhLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.nuh; GTMuLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u; GTMvLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.v; timeLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.time
GTMtempLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.temp; GTMdepthLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.z; GTMnumLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.num; GTMnuhLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.nuh; GTMuLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u; GTMvLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.v; timeLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.time


##################
## convert time ##

time2 = 267.75 + (time[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
time2Km = 267.75 + (timeKm[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
time2Kh = 267.75 + (timeKh[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
time2LaRib = 267.75 + (timeLaRib[:]/86400); 

time2LaRibhLPveP100 = 267.75 + (timeLaRibhLPveP100[:]/86400); 
time2LaRibhLPveP1500 = 267.75 + (timeLaRibhLPveP1500[:]/86400); 
time2LaRibhLPveP10000 = 267.75 + (timeLaRibhLPveP10000[:]/86400); 

# time_list = [time2,time2Km,time2Kh,time2LaRib,time2LaRibhLPveP100,time2LaRibhLPveP1500,time2LaRibhLPveP10000]

# def julian_to_date_GOTM():
# 	initial_date = "2012-01-01 00:00:00"
# 	for i in time_list:
# 		for j in i.values:
# 			GTMdate = pd.to_datetime(initial_date) + pd.DateOffset(days=j)
# 			print(date)
# 	return print('GOTM Time format is in Date-Time == SUCCESS')

# julian_to_date_GOTM()
# print()
###################
# time allocation #

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkppKm = gotmkppKm.assign(time2Km=("time", time2Km)); gotmkppKm = gotmkppKm.swap_dims({"time" : "time2Km"})
gotmkppKh = gotmkppKh.assign(time2Kh=("time", time2Kh)); gotmkppKh = gotmkppKh.swap_dims({"time" : "time2Kh"})
gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})

gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.assign(time2LaRibhLPveP100=("time", time2LaRibhLPveP100)); gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.swap_dims({"time" : "time2LaRibhLPveP100"})
gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.assign(time2LaRibhLPveP1500=("time", time2LaRibhLPveP1500)); gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.swap_dims({"time" : "time2LaRibhLPveP1500"})
gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.assign(time2LaRibhLPveP10000=("time", time2LaRibhLPveP10000)); gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.swap_dims({"time" : "time2LaRibhLPveP10000"})


#################################
### SECTION GOTM DIAGNOSE MLD ###
#################################

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
## use multiple of 2 for Km ##
mld_tempKPPKm = GTMtempKm.sel(z=-11.0, method='nearest') - 0.2;
mld_tempKPPKm = (GTMtempKm-mld_tempKPPKm).isel(lat=0,lon=0).values

for time in mld_tempKPPKm:
	for depth in range(len(time)):
		if time[depth] > 0:
			MLDKmindx.append(depth)
			break
MLDKmindx = np.asarray(MLDKmindx); mldKPPKm = GTMtempKm.z[MLDKmindx];

##############################
## use multiple of 2 for Kh ##
mld_tempKPPKh = GTMtempKh.sel(z=-11.0, method='nearest') - 0.2;
mld_tempKPPKh = (GTMtempKh-mld_tempKPPKh).isel(lat=0,lon=0).values

for time in mld_tempKPPKh:
	for depth in range(len(time)):
		if time[depth] > 0:
			MLDKhindx.append(depth)
			break
MLDKhindx = np.asarray(MLDKhindx); mldKPPKh = GTMtempKh.z[MLDKhindx];

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


##################################################################
## use Lat and Rib original with hL condition (Ch=1 for hL>100) ##
mld_tempKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
mld_tempKPPLaRibhLPveP100 = (GTMtempLaRibhLPveP100-mld_tempKPPLaRibhLPveP100).isel(lat=0,lon=0).values

for time in mld_tempKPPLaRibhLPveP100:
	for depth in range(len(time)):
		if time[depth] > 0:
			MLDKPPLaRibhLPveP100indx.append(depth)
			break
MLDKPPLaRibhLPveP100indx = np.asarray(MLDKPPLaRibhLPveP100indx); 
mldKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.z[MLDKPPLaRibhLPveP100indx];


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


##############################################################################################
## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma ##
mld_tempKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.sel(z=-11.0, method='nearest') - 0.2; 
mld_tempKPPLaRibhLPveP10000 = (GTMtempLaRibhLPveP10000-mld_tempKPPLaRibhLPveP10000).isel(lat=0,lon=0).values

for time in mld_tempKPPLaRibhLPveP10000:
	for depth in range(len(time)):
		if time[depth] > 0:
			MLDKPPLaRibhLPveP10000indx.append(depth)
			break
MLDKPPLaRibhLPveP10000indx = np.asarray(MLDKPPLaRibhLPveP10000indx); 
mldKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.z[MLDKPPLaRibhLPveP10000indx];



# ## potential definition
# def diagnoseMLD():
# 	for z_level,time in range(GTMtemp.shape[1]):
# 		if GTMtemp[z_level,time] > 0:
# 			MLDindx = z_level
# 			break
# 		MLD[time] = np.where(MLDindx)
# 	print(MLD)

# diagnoseMLD()

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

# print(j)
# to get a list, instead of a generator, use
# x_y = list(read__lines())
## Plot of MLD from OSMOSIS Observations
# plt.plot(j,np.negative(i),color='olivedrab',linestyle='dotted'); 
# plt.xlabel('Time -Julian day'); plt.ylabel('Mixed Layer Depth -m')
# plt.grid()
# plt.axis([260,341,-400,0]); 
# plt.title('OSMOSIS Mixed Layer Depth')
# plt.savefig('OSMOSIS_Mixed_Layer_Depth_.png'); plt.show()

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


## ----------------------------------------------------------------------------------------- ##
###################################################################
## ------------------- momentum flux times --------------------- ##
###################################################################

concat_func = lambda x,y: str(x) + " " + str(y)
mflux = np.genfromtxt('../../gotm-4.0.0-edit/simulations/annual_edited/momentum_flux.dat')

mdate = []; mtime = []; xTau = []; yTau = []; xStokes = []; yStokes = []; depthStokes = []

with open('../../gotm-4.0.0-edit/simulations/annual_edited/momentum_flux.dat', 'r') as f:
	for row in f:
		r, s, t, u, v, w, z = row.split()
		mdate.append(r)
		mtime.append(s)
		xTau.append(t)
		yTau.append(u)
		xStokes.append(v)
		yStokes.append(w)
		depthStokes.append(z)
depthStokes = np.asarray(depthStokes); depthStokes = depthStokes.astype('float64')

# print(mdate); print(mtime)
DT_m = list(map(concat_func,mdate,mtime)) # list the map function
DT_m = np.asarray(DT_m)
# print('DT_m',DT_m)
mdate = np.asarray(mdate);# mdate = mdate.astype('str')

mtime = np.asarray(mtime)

# print(DT.shape)
yearDaysM = np.empty((1465,0), int)
yearmonth = np.empty((1465,0), int)
# yearM = np.empty((1465,0), int)
# monthM = np.empty((1465,0), int)
# dayM = np.empty((1465,0), int)
# datess = np.empty((1465,0))

for i in range(len(DT_m)):
	yearDaysM= np.append(yearDaysM, datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_yday)
	yr = datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_year
	mth = datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_mon
	dys = datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_mday
	hrs = datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_hour
	mins = datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_min
	sec = datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_sec
	yearmonth = np.append(yearmonth , str(yr) + "-" + str(mth) + "-" + str(dys) + ' ' + str(hrs) + ":" + str(mins) + ":" + str(sec))
	# datess= np.append(datess, datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetup

la_t = np.zeros(1465)

#---------------------------------------------------------------------------------------------#
#-----------------------------------------------------#
#          Surface sea temperature - degreeC          #

# OSMOSIS
GldrTemp = GldrTemp[:,2] # [2]=[-5 m]
# print('OSMOSIS Temp at -5 m : ',GldrTemp)
# print('OSMOSIS depth : ',GldrDepth.values)

# GldrTemp = GldrTemp.sel(pressure=-1)
#  Original KPP scheme
GTMtemp = GTMtemp.isel(lon=0, lat=0)
GTMtemp = GTMtemp.sel(z=-5, method='nearest')
#  Original KPP Lat Rib
GTMtempLaRib = GTMtempLaRib.isel(lon=0, lat=0)
GTMtempLaRib = GTMtempLaRib.sel(z=-5, method='nearest')
# #  KPP nuh*2 Rib*2
# GTMtempRib2nuh2 = GTMtempRib2nuh2.isel(lon=0, lat=0)
# GTMtempRib2nuh2 = GTMtempRib2nuh2.sel(z=-5, method='nearest')
# KPP with La_T Ri_B hL>100
GTMtempLaRibhLPveP100 = GTMtempLaRibhLPveP100.isel(lon=0, lat=0)
GTMtempLaRibhLPveP100 = GTMtempLaRibhLPveP100.sel(z=-5, method='nearest')
# KPP with La_T Ri_B hL>1500
GTMtempLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.isel(lon=0, lat=0)
GTMtempLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.sel(z=-5, method='nearest')
# KPP with La_T Ri_B hL>10000
GTMtempLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.isel(lon=0, lat=0)
GTMtempLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.sel(z=-5, method='nearest')

#-----------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#-----------------------------------------------------#
#          Surface sea temperature - degreeC          #
#-----------------------------------------------------#

GTM_hL_LaRib = gotmkppLaRib.hL; GTM_hLR_LaRib = gotmkppLaRib.hLR; GTM_Lat_LaRib = gotmkppLaRib.La_t; GTM_tdotus_LaRib = gotmkppLaRib.tdotus; 
GTM_hL_LaRibhLPveP100 = gotmkppLaRibhLPveP100.hL; GTM_hLR_LaRibhLPveP100 = gotmkppLaRibhLPveP100.hLR; GTM_Lat_LaRibhLPveP100 = gotmkppLaRibhLPveP100.La_t; GTM_tdotus_LaRibhLPveP100 = gotmkppLaRibhLPveP100.tdotus; 
# GTM_hL_LaRibhLPveN100 = gotmkppLaRibhLPveN100.hL; GTM_hLR_LaRibhLPveN100 = gotmkppLaRibhLPveN100.hLR; GTM_Lat_LaRibhLPveN100 = gotmkppLaRibhLPveN100.La_t; GTM_tdotus_LaRibhLPveN100 = gotmkppLaRibhLPveN100.tdotus; 
GTM_hL_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.hL; GTM_hLR_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.hLR; GTM_Lat_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.La_t; GTM_tdotus_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.tdotus; GTM_Theatflux_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.total; GTM_fric_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u_taus; GTM_us_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u_s; 
# GTM_hL_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.hL; GTM_hLR_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.hLR; GTM_Lat_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.La_t; GTM_tdotus_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.tdotus; GTM_Theatflux_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.total; GTM_fric_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u_taus;


################################
# extract data for hL analysis #
################################

# GTM_hLR_LaRibhL2 = GTM_hLR_LaRibhL2.isel(lon=0,lat=0)
GTM_hLR_LaRib = GTM_hLR_LaRib.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP100 = GTM_hLR_LaRibhLPveP100.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP1500 = GTM_hLR_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP10000 = GTM_hLR_LaRibhLPveP10000.isel(lon=0,lat=0)



############################
### SECTION VISUAL PLOTS ###
############################
fig, ax = plt.subplots(1,1,sharex=True,figsize=(40,18))

# # ##########################
# # # month-year on top-axis #

ax0=ax.twiny(); 
ax0.tick_params(axis='both', which='major', labelsize=18)
ax0.set_axisbelow(True)
# ax0.set_xlim(time2[0],time2[-1])
xa0 = ax0.get_xaxis(); 
xa0.axis_date(); 
xa0.set_minor_locator(mdates.DayLocator(interval=1)); 
xa0.set_major_locator(mdates.MonthLocator()); 
xa0.set_major_formatter(mdates.DateFormatter('%b-%y'))
# ax0.plot(new_time2,GldrTemp,color='white',linestyle='dotted')
ax0.plot(yearmonth,la_t,color='white')
# ax0.set_xlim(480,540)

# # compare Km and Kh multiples #

# ax.plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = r'OSMOSIS')
# ax.plot(time2,mldKPP,color='black',linestyle='solid',label = r'GOTM KPP: $K_m$, $K_h$')
# ax.plot(time2Km,mldKPPKm,color='blue',linestyle='solid',label = r'GOTM KPP: $2*K_m$, $K_h$')
# ax.plot(time2Kh,mldKPPKh,color='red',linestyle='solid',label = r'GOTM KPP: $K_m$; $2*K_h$')
# ax.grid(); ax.legend(loc='lower right',fontsize=25)
# ax.set_xlim(485,540); ax.set_ylim(-300,0)
# ax.tick_params(axis='both', which='major', labelsize=23)
# ax.set_ylabel(r'Mixed Layer Depth -$m$',fontsize=30); ax.set_xlabel(r'Time -Julian days',fontsize=30)

# plt.show()
# plt.savefig("compare_Kh_Km_Fig.2_paper.png")

### compare h/L thresholds ###
##############################

# SST
ax.plot(new_time,GldrTemp,color='indigo',linestyle='dashed',label ='OSMOSIS')
ax.plot(time2,GTMtemp,color='black',linestyle='solid',label =r'KPP $La_{T}$ OFF')
ax.plot(time2,GTMtempLaRib,color='goldenrod',linestyle='solid',label =r'KPP $La_{T}$ $\gamma = \delta = 1$')
ax.plot(time2LaRibhLPveP100,GTMtempLaRibhLPveP100,color='limegreen',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ (h/L_L)_{MAX}=6.1$')
ax.plot(time2LaRibhLPveP1500,GTMtempLaRibhLPveP1500,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ (h/L_L)_{MAX}=90.9$')
ax.plot(time2LaRibhLPveP10000,GTMtempLaRibhLPveP10000,color='dodgerblue',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ (h/L_L)_{MAX}=606.2$')
ax.grid(); ax.legend(loc='upper left',fontsize=25)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlim(265,635); ax.set_ylim(10,24)
# ax[1].set_ylim(11,23)
ax.set_ylabel(r'SST -$\degree C$',fontsize=30)
ax.set_xlabel(r'Time -Julian days',fontsize=30)

plt.show()
plt.savefig("compare_SST_Fig.8_paper.png")

# # MLD

# ax[0].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = r'OSMOSIS')
# ax[0].plot(time2,mldKPP,color='black',linestyle='solid',label = r'KPP $La_{T}$ OFF')
# ax[0].plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta = 1$')
# ax[0].plot(time2LaRibhLPveP100,mldKPPLaRibhLPveP100,color='limegreen',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta=1$; $(h/L_L)_{MAX}=6.1$')
# ax[0].grid(); ax[0].legend(loc='lower right',fontsize=21)
# # ax[0].set_xlim(265,635);# ax[0].set_ylim(-200,0)
# ax[0].set_xlim(480,550); ax[0].set_ylim(-150,0)
# ax[0].tick_params(axis='both', which='major', labelsize=20)
# # ax[0].set_ylabel(r'Mixed Layer Depth -$m$')


# ax[1].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = '_nolegend_')
# ax[1].plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
# ax[1].plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'_nolegend_')
# ax[1].plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta = 1$; $(h/L_L)_{MAX}=90.9$')
# ax[1].grid(); ax[1].legend(loc='lower right',fontsize=21)
# # ax[1].set_xlim(265,635);# ax[1].set_ylim(-200,0)
# ax[1].set_xlim(480,550); ax[1].set_ylim(-150,0)
# ax[1].set_ylabel('Mixed Layer Depth -m',fontsize=25)
# ax[1].tick_params(axis='both', which='major', labelsize=20)


# ax[2].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = '_nolegend_')
# ax[2].plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
# ax[2].plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'_nolegend_')
# ax[2].plot(time2LaRibhLPveP10000,mldKPPLaRibhLPveP10000,color='dodgerblue',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta = 1$; $(h/L_L)_{MAX}=606.2$')
# ax[2].grid(); ax[2].legend(loc='lower right',fontsize=21)
# # ax[2].set_xlim(265,635);# ax[2].set_ylim(-200,0)
# ax[2].set_xlim(480,550); ax[2].set_ylim(-150,0)
# ax[2].set_xlabel('Time -Julian days',fontsize=25)
# ax[2].tick_params(axis='both', which='major', labelsize=20)

# # plt.legend(loc='lower right',fontsize='small')

# # fig.text(0.5, 0.01, r'Time -Julian days', ha='center', va='center', fontsize=12)
# # fig.text(0.008, 0.5, r'Mixed Layer Depth -$m$', ha='center', va='center', fontsize=12, rotation='vertical')
# # # # ax.text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
# ax[0].text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes,fontsize=20)
# ax[1].text(0.02, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes,fontsize=20)
# ax[2].text(0.02, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes,fontsize=20)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# # plt.suptitle(r'Comparison of mixed layer depths with $La_T$ and $h/L$ thresholds')
# plt.show()
# # plt.savefig('compare_MLDs_Fig.6_paper.png')
# plt.savefig('compare_MLDs_Fig.7_paper.png')

exit()





### OLDER VERSION
###############################################################################################
### ----------------------------- calculate mixed layer depth ----------------------------- ###
###############################################################################################

## Psuedocode 
# 1. extract temperature values for depth z_ref
# 2. new_temp = temp(z_ref) - delta_t
# 3. depth at new_temp

# MLD_DT02 = depth where (T = T(z_ref) ± delta_t)
# z_ref = -10 -- reference depth;   T(z_ref) -- temperature at reference depth          
# delta_t = 0.2 -- ;    


## using KPP model
mld_tempKPP = GTMtemp.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes1 = abs(GTMtemp-mld_tempKPP).argmin('z');#  print(z_indexes.values)
mldKPP = GTMtemp.z[z_indexes1];
mldKPP = mldKPP.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original
mld_tempKPPLaRib = GTMtempLaRib.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes2 = abs(GTMtempLaRib-mld_tempKPPLaRib).argmin('z'); # print(z_indexes)
mldKPPLaRib = GTMtempLaRib.z[z_indexes2]; 
mldKPPLaRib = mldKPPLaRib.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100)
mld_tempKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes3 = abs(GTMtempLaRibhLPveP100-mld_tempKPPLaRibhLPveP100).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes3.isel(lat=0,lon=0)[55])
mldKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.z[z_indexes3]; 
mldKPPLaRibhLPveP100 = mldKPPLaRibhLPveP100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>1500) for gamma
mld_tempKPPLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes4 = abs(GTMtempLaRibhLPveP1500-mld_tempKPPLaRibhLPveP1500).argmin('z'); # print(z_indexes)
print((GTMtempLaRibhLPveP1500-mld_tempKPPLaRibhLPveP1500).isel(lat=0,lon=0)[55,469:499])
mldKPPLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.z[z_indexes4]; 
mldKPPLaRibhLPveP1500 = mldKPPLaRibhLPveP1500.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes5 = abs(GTMtempLaRibhLPveP10000-mld_tempKPPLaRibhLPveP10000).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.z[z_indexes5]; 
mldKPPLaRibhLPveP10000 = mldKPPLaRibhLPveP10000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# testing


timeMLD_KPP = np.c_[time2[55], mldKPP[55], GTMtemp.isel(lon=0,lat=0)[0,55]]
timeMLD_P100 = np.c_[time2LaRibhLPveP100[55], mldKPPLaRibhLPveP100[55], GTMtempLaRibhLPveP100.isel(lon=0,lat=0)[0,55]]
timeMLD_P1500 = np.c_[time2LaRibhLPveP1500[55], mldKPPLaRibhLPveP1500[55], GTMtempLaRibhLPveP1500.isel(lon=0,lat=0)[0,55]]
timeMLD_P10000 = np.c_[time2LaRibhLPveP10000[55], mldKPPLaRibhLPveP10000[55], GTMtempLaRibhLPveP10000.isel(lon=0,lat=0)[0,55]]


print('MLD original KPP :', timeMLD_KPP)
print(r'MLD KPP $(h/L > 100)$ :', timeMLD_P100)
print(r'MLD KPP $(h/L > 1500)$ :', timeMLD_P1500)
print(r'MLD KPP $(h/L > 10000)$ :', timeMLD_P10000)
