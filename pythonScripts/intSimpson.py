### ------------------- NUMERICAL INTEGRATION - SIMPSONS' RULE ---------------------- ###
import numpy as np
import scipy
from scipy import integrate
import math
import sys
import os.path
import xarray as xr
import pandas as pd
# from definitionsOSMOSIS import ctsolver, visualize, ctvariables
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
# np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.for1mat})

### ----------------------------- GLIDER ----------------------------- ##
# file1 = 'C:/home/users/mc837749/Documents/gotm-4.0.0/simulations/annualOSMOSISforcing_2502/glider_timeseries.nc'
gliderdata= Dataset('../glider_timeseries.nc'); glider = xr.open_dataset('../glider_timeseries.nc', decode_times=False)

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; Gldrdays = glider.dats; GldrSalt = glider.prac_sal;

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days

glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})

# print(gliderdata.variables.items())

### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("osmosis_annual_surface_forcingKPP/OSMOSIS_glider_comparison.nc", decode_times=False)
gotmkeps = xr.open_dataset("osmosis_annual_surface_forcingkeps/OSMOSIS_glider_comparison.nc", decode_times=False)

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; time = gotmkpp.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; ketime = gotmkeps.time
## Convert time
time2 = 261.0 + (time[:]/86400); ketime2 = 261.0 + (ketime[:]/86400)
# times = np.c_[time,time2]
# print('GOTM times :', time2.values);
gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})
#### --------------- testing CT difference ---------------- ####
# # print(gotmkpp.variables.items())
# print('shape of kpp temp', GTMtemp.shape); print('kpp depth',GTMdepth[0].values); 
# print('shape of kpp time', time2[:143].values)
# # print('shape of glider temp', GldrTemp[0:5,0:50].values, GldrTemp[0:5,0:50].shape);
# print('shape of kpp temp', GTMtemp[0:5,400:500].values, GTMtemp[0:5,400:500].shape);
# # print('shape of glider depth', GldrDepth[0:50].values, GldrDepth[0:50].shape); 
# print('shape of kpp depth', GTMdepth[400:500].values, GTMdepth[400:500].shape); 
# # print('shape of glider time', new_time[6:100], new_time[6:100].shape)
# print('shape of kpp time', time2[143:1007].values, time2[143:1011].shape)
# # Gldrtemp_new = GldrTemp[0:4,0:5]; Gldrz_new = GldrDepth[0:5]; Gldrtime_new = new_time[0:4]
# # print(Gldrtemp_new); print(Gldrz_new); print(Gldrtime_new)

## sellecting the temperature for the mixed layer depth
dz = np.diff(GTMdepth)[0];
kpptemp2dMLD = GTMtemp.isel(lat=0,lon=0,z=range(400,499));# print(kpptemp2dMLD.shape)
sumkppTempMLD = kpptemp2dMLD.sum(axis=1)*dz;# print(sumkppTempMLD.shape)

heatBudget = []
midtime = []
	# for i in range(1,Nz):
	# 	sumTempMLD[i] = np.sum(GTMtempA.isel(lon=0,lat=0)[i,400:500]) 
dt = np.diff(time2)[0]*86400;

for n in range(1, 52595):
	CT = (sumkppTempMLD[n+1].values-sumkppTempMLD[n].values)/dt
	# CT = (GTMtempA.isel(lon=0,lat=0)[n+1,:]-GTMtempA.isel(lon=0,lat=0)[n,:])/dt
	diffT = (time2[n+1].values+time2[n].values)/2
	heatBudget.append(CT)
	midtime.append(diffT)

# print(temp.shape); print(time.shape)

plt.plot(midtime, heatBudget, 'k-.'); plt.xlabel('Time -Julian days'); plt.ylabel('heat budget for mixed layer depth')
plt.axis([255,630,-0.0002,0.0003])
plt.savefig('heatFlux_osmosisAnnual/derived_heat_budget_mixed_layer_depth.png'); plt.show()


exit()
## call the definition
dt = 1
T = 626.25 # 52595
# dz = 1; L = 0
temp, time, t = ctsolver(dt,T)
print('rate of change of temperature in time:', temp.values)
print('change in time:', time)
exit()

#### --------------- temperature timeseries ---------------- ####

# Gldrdepth_ = GldrDepth[200].values; print(Gldrdepth_)

a = [0,25,50,75,100,125,150,175,200]
# a = np.arange(0,504)
for j in a:
	Gldrtemp1d = GldrTemp.isel()[:,j]; Gldrdepth = GldrDepth[j]
	Gldrtime = new_time[:] 
	plt.plot(Gldrtime,Gldrtemp1d,label='Glider depth %s'%Gldrdepth.values); plt.legend(loc='upper center', fontsize='xx-small')
	plt.xlabel('Julian -days');plt.ylabel('Glider temperature -C'); 
	# plt.title('Compared Temperature Profiles : Glider and GOTM w KPP,k-eps')
	plt.axis([260,600,10,25])
plt.savefig('Gldr__Temp_timeseries.png')
plt.show()

b = [299,324,349,374,399,424,449,474,499]
# b = np.arange(0,499)
for i in b:
	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
	kpptime = time2[:]
	plt.plot(kpptime,kpptemp1d,label = 'GOTM_KPP depth %s'%kppdepth.values); plt.legend(loc='upper center', fontsize='xx-small')
	plt.xlabel('Julian days -days');plt.ylabel('KPP GOTM temperature -C'); 
	plt.axis([250,650,8,25])
plt.savefig('GOTM_kpp__Temp_timeseries.png')
plt.show()

for i in b:
	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[:,i]; kepsdepth = GTMkedepth[i]
	kepstime = ketime2[:]
	plt.plot(kepstime,kepstemp1d,label = 'GOTM_k-eps depth %s'%kepsdepth.values); plt.legend(loc='upper center', fontsize='xx-small')
	plt.xlabel('Julian days -days');plt.ylabel('k-eps GOTM temperature -C'); 
	plt.axis([250,650,8,25])
plt.savefig('GOTM_keps__Temp_timeseries.png')
plt.show()


exit()
#### --------------- testing values ---------------- ####

# GTMdepth_ = GTMdepth.values; print(GTMdepth_)
GTMtemp400 = GTMtemp.isel(lon=0,lat=0)[299,:].values; print('temp values at -400m',GTMtemp400,len(GTMtemp400))
GTMtemp350 = GTMtemp.isel(lon=0,lat=0)[324,:].values; print('temp values at -350m',GTMtemp350,len(GTMtemp350))
GTMtemp300 = GTMtemp.isel(lon=0,lat=0)[349,:].values; print('temp values at -300m',GTMtemp300,len(GTMtemp300))
GTMtemp250 = GTMtemp.isel(lon=0,lat=0)[374,:].values; print('temp values at -250m',GTMtemp250,len(GTMtemp250))
GTMtemp200 = GTMtemp.isel(lon=0,lat=0)[399,:].values; print('temp values at -200m',GTMtemp200,len(GTMtemp200))
GTMtemp150 = GTMtemp.isel(lon=0,lat=0)[424,:].values; print('temp values at -150m',GTMtemp150,len(GTMtemp150))
GTMtemp100 = GTMtemp.isel(lon=0,lat=0)[449,:].values; print('temp values at -100m',GTMtemp100,len(GTMtemp100))
GTMtemp50  = GTMtemp.isel(lon=0,lat=0)[474,:].values; print('temp values at  -50m', GTMtemp50,len(GTMtemp50))
GTMtemp0   = GTMtemp.isel(lon=0,lat=0)[499,:].values; print('temp values at    0m',  GTMtemp0,len(GTMtemp0))

dx = (GTMdepth[499].values-GTMdepth[299].values)/50;# print(dx)
Simps = dx/3 * (GTMtemp400 + 4*GTMtemp350 + 2*GTMtemp300 + 4*GTMtemp250 + 2*GTMtemp200 + 4*GTMtemp150 + 2*GTMtemp100 + 4*GTMtemp50 + GTMtemp0)
print('integral of temp calculated by Simpsons rule: ',Simps)

exit()
#### -------------- PLOTS ------------------ ####

for j in range(len(Gldr_sel)):
	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
	Gldrtime = new_time[j] 
	for i in range(len(GTM_sel)):
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
		kpptime = time2[i]
		kepstime = ketime2[i]

exit()
def simps(f,a,b,N=50):
    '''Approximate the integral of f(x) from a to b by Simpson's rule.

    Simpson's rule approximates the integral int_a^b f(x) dx by the sum:
    (dx/3) sum_{k=1}^{N/2} (f(x_{2i-2} + 4f(x_{2i-1}) + f(x_{2i}))
    where x_i = a + i*dx and dx = (b - a)/N.

    Parameters
    ----------
    f : function
        Vectorized function of a single variable
    a , b : numbers
        Interval of integration [a,b]
    N : (even) integer
        Number of subintervals of [a,b]

    Returns
    -------
    float
        Approximation of the integral of f(x) from a to b using
        Simpson's rule with N subintervals of equal length.

    Examples
    --------
    >>> simps(lambda x : 3*x**2,0,1,10)
    1.0
    '''
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    dx = (b-a)/N
    x = np.linspace(a,b,N+1)
    y = f(x)
    S = dx/3 * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2])
    return S


ex1 = simps(lambda x : 3*x**2,0,1,10)
print(ex1)
ex2 = simps(lambda x : x**0.5,0,8,4)
print(ex2)


