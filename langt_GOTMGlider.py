## ------------- Written code to examine the periods with high and low La_t ------------- ##
import numpy as np
import math
import os
import xarray as xr
import pandas as pd
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt

### ----------------------------- GLIDER ----------------------------- ##

gliderdata= Dataset('glider_timeseries.nc'); glider = xr.open_dataset("glider_timeseries.nc", decode_times=False)

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; Gldrdays = glider.dats

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days

glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})
# print('Glider time from 2012-09-22',new_time[194:938].values)

### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("OSMOSIS_cruise_winter.nc", decode_times=False)
gotmkeps = xr.open_dataset("../201220_keps/OSMOSIS_cruise_winter.nc", decode_times=False)

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; time = gotmkpp.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; ketime = gotmkeps.time

## Convert time
time2 = 261 + (time[:]/86400); ketime2 = 261 + (ketime[:]/86400)
gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})

# print('GOTM time from 2012-09-17',time2.shape)
# print('GOTM time from 2012-09-17',ketime2.shape)
# print('Glider time from 2012-09-04',new_time[143:950])

# -------------- PLOTS -------------- #
# rvL = [0:45]
### plot
for i in rvL:
	# kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	# kpptime = time2[i]; 
	kepstime = ketime2[i]
	# plt.plot(kpptemp1d,kppdepth,label = 'GOTM_kpp day %s'%kpptime.values); plt.legend(loc=(0.45,0.35))
	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc=(0.45,0.35))
plt.xlabel('Temp -C'); plt.ylabel('Depth -m'); plt.title('GOTM_keps_test Low La_t values')
plt.axis([11,18,-400,0])
plt.savefig('GOTM_keps_Temp_lowLat_test.png')
plt.show()

exit()


## select desired rows for the loop to plot specific days
# Gldr_sel = [949] # range [143:950]
# GTM_sel = [5696] # range [0:5696]

# ## Compared plot GLider and GOTM_KPP
# for j in Gldr_sel:
# 	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
# 	Gldrtime = new_time[j] 
# 	for i in GTM_sel:
# 		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
# 		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
# 		kpptime = time2[i]
# 		kepstime = ketime2[i]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_KPP Day %s'%kpptime.values); plt.legend(loc=(0.45,0.35))
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps Day %s'%kepstime.values); plt.legend(loc=(0.45,0.35))
# 	plt.plot(Gldrtemp1d,Gldrdepth,label='Glider Day %s'%Gldrtime.values); plt.legend(loc=(0.45,0.35))
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	plt.title('Compared Temperature Profiles : Glider and GOTM w KPP,k-eps')
# 	plt.axis([11,18,-400,0])
# plt.savefig('Gldr_GOTM_kppkeps_Temp_profile_select_2412test.png')
# # plt.show()
exit()

### ----------------------------- .dat file data ----------------------------- ###

hdate = []; htime = []

with open('heat_flux.dat', 'r') as f:
	for row in f:
		a, b, c, d = row.split()
		hdate.append(a)
		htime.append(b)
# print(hdate); print(htime)
concat_func = lambda x,y: str(x) + " " + str(y)

DT = list(map(concat_func,hdate,htime)) # list the map function
DT = np.asarray(DT)
hdate = np.asarray(hdate); htime = np.asarray(htime)

yearDays = np.empty((634,0), int)
for i in range(len(DT)):
	yearDays= np.append(yearDays, datetime.strptime(DT[i],'%Y-%m-%d %H:%M:%S').timetuple().tm_yday)
# print(yearDays)
hours = np.empty((634,0), int)
for tim in range(len(htime)):
	timeDays = datetime.strptime(htime[tim], '%H:%M:%S')#.split(':')
	hours = np.append(hours, timeDays.hour)

hourstodays = hours/24
# print(hourstodays)
convertedDays = yearDays + hourstodays
print('Converted Days: ', convertedDays)

## CALCULATIONS

hflux = np.genfromtxt('heat_flux.dat')
mflux = np.genfromtxt('momentum_flux.dat')
print('Heat Flux Data: ', hflux); print(hflux.shape)


## add the converted days as a new column to the array
heatFlux = np.c_[convertedDays,hflux]; # print(heatFlux)

momFlux = np.c_[convertedDays,mflux]

## calculate u_s and u_*
# u_s = math.sqrt((momFlux[0,3])**2 + (momFlux[0,4])**2)
# print(momFlux[0][3],momFlux[0,4]); print(u_s);

u_s = np.empty((634,0), float)
fric_u = np.empty((634,0), float)
## Density of air == 1.225 k/m**3
for time in range(len(momFlux)):
	u_s = np.append(u_s, math.sqrt((momFlux[time][5])**2 + (momFlux[time][6])**2))
	fric_u = np.append(fric_u, math.sqrt(((momFlux[time][3])**2 + (momFlux[time][4])**2)/1.225))
print('Stoke drift veloctiy :',u_s); print(u_s.shape)
print('friction velocity :',fric_u); print(fric_u.shape)

## calculate la_t

la_t = np.empty((634,0), float)
for i in range(len(u_s)):
	la_t = np.append(la_t, math.sqrt(fric_u[i]/u_s[i]))
print('turbulent Langmuir number :',la_t); print(la_t.shape)

## mean value of la_t

la_tMean = np.mean(la_t)
print('mean value of La_t', la_tMean)

##  values of la_t
la_tt = np.c_[convertedDays,la_t]; print('added day of year to la_t', la_tt)
print(la_tt.shape)
la_tMax = np.amax(la_tt); print('maximum value of la_t',la_tMax)
la_tMin = np.amin(la_tt); print('minimum value of la_t',la_tMin)

## ------------ period with low La_t values ------------ ##
# rows, cols = np.where(la_tt < 0.6); print('values below la_t=0.6',la_tt[rows])
# la_tL = la_tt[rows][:,0]; print('times of low la_t values: ', la_tL);
# rowvalL = np.in1d(time2, la_tL).nonzero()[0]; rowval_low = rowvalL.tolist(); print(rowval_low)
# [35, 44, 53, 62, 71, 80, 89, 98, 215, 791, 800, 1187, 1475, 1799, 1808, 1817, 1826, 1835, 2249, 2258, 
## 2339, 2348, 2393, 2870, 2879, 2888, 3050, 3059, 3563, 3572, 3635, 3644, 3662, 3671, 3680, 
### 3689, 3698, 3707, 4067, 4076, 4085, 4094, 4373, 4814, 4823, 5354, 5363, 5624, 5633]
rvL = [3059]
### plot
# for i in rvL:
# 	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
# 	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
# 	kpptime = time2[i]; 
# 	kepstime = ketime2[i]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_kpp day %s'%kpptime.values); plt.legend(loc=(0.45,0.35))
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc=(0.45,0.35))
# plt.xlabel('Temp -C'); plt.ylabel('Depth -m'); plt.title('GOTM_kpp_kep Low La_t values')
# plt.axis([11,18,-400,0])
# plt.savefig('GOTM_kpp_keps_Temp_lowLat_test.png')
# plt.show()

exit()
## ------------ period with high La_t values ------------ ##

la_tHigh = np.where(la_tt[:,1] >= [1.4]); # print('values above la_t=1.4',la_tHigh)
for x in la_tHigh:
	rowsH = la_tt[x]
print('extract high values in 2d :', rowsH)
la_tH = rowsH[:,0]; print('times of high la_t values: ', la_tH);
rowvalH = np.in1d(time2, la_tH).nonzero()[0]; rowval_high = rowvalH.tolist(); print(rowval_high)
## [1205, 1214, 1223, 1232, 1241, 1250, 1259, 1268, 1367, 1583, 1610, 1853, 1961, 2087, 2276, 2285, 
###  2294, 2303, 2411, 2447, 2906, 2915, 4121, 4220, 4229, 5669, 5678, 5687]
rvH = [2276]
## plot
# for j in rvH:
# 	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[j,:]; kppdepth = GTMdepth[:]
# 	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[j,:]; kepsdepth = GTMkedepth[:]
# 	kpptime = time2[j]; 
# 	kepstime = ketime2[j]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_kpp day %s'%kpptime.values); plt.legend(loc=(0.45,0.35))
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc=(0.45,0.35))
# plt.xlabel('Temp -C'); plt.ylabel('Depth -m'); plt.title('GOTM_kpp_keps High La_t values')
# plt.axis([11,18,-400,0])
# plt.savefig('GOTM_kpp_keps_Temp_highLat_test.png')
# plt.show()

print('low values of la_t (less than 1)', la_tt[0:30,:])

exit()