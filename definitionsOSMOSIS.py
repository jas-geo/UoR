import numpy as np
import math
import sys
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

gliderdata= Dataset('../glider_timeseries.nc'); glider = xr.open_dataset("../glider_timeseries.nc", decode_times=False)

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; Gldrdays = glider.dats; GldrSalt = glider.prac_sal;

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days

glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})

### ----------------------------- GOTM ----------------------------- ###

## ANNUAL KPP ##
################

gotmkppA = xr.open_dataset("osmosis_annual_surface_forcingKPP/OSMOSIS_glider_comparison.nc", decode_times=False)

GTMtempA = gotmkppA.temp; GTMdepthA = gotmkppA.z; timeA = gotmkppA.time

kpptime2 = 261.0 + (timeA[:]/86400); ## Convert time
gotmkppA = gotmkppA.assign(kpptime2=("time", kpptime2)); gotmkppA = gotmkppA.swap_dims({"time" : "kpptime2"})

##MLD
## selecting the temperature for the mixed layer depth
kpptemp2dMLD = GTMtempA.isel(lat=0,lon=0,z=range(400,499))
sumkppTempMLD = kpptemp2dMLD.sum(axis=1)

## SPRING LOW LAT ##
####################

gotmkppSpr = xr.open_dataset("spring_lowLatKPP/OSMOSIS_glider_comparison.nc", decode_times=False)
gotmkepsSpr = xr.open_dataset("spring_lowLatkeps/OSMOSIS_glider_comparison.nc", decode_times=False)

GTMtempSpr = gotmkppSpr.temp; GTMdepthSpr = gotmkppSpr.z; timeSpr = gotmkppSpr.time
GTMketempSpr = gotmkepsSpr.temp; GTMkedepthSpr = gotmkepsSpr.z; ketimeSpr = gotmkepsSpr.time

time2 = 493.0 + (timeSpr[:]/86400); ketime2 = 493.0 + (ketimeSpr[:]/86400); ## Convert time
gotmkppSpr = gotmkppSpr.assign(time2=("time", time2)); gotmkppSpr = gotmkppSpr.swap_dims({"time" : "time2"})
gotmkepsSpr = gotmkepsSpr.assign(ketime2=("time", ketime2)); gotmkepsSpr = gotmkepsSpr.swap_dims({"time" : "ketime2"})


## SPRING high LAT ##
#####################

gotmkpp = xr.open_dataset("spring_highLatKPP/OSMOSIS_glider_comparison.nc", decode_times=False)
gotmkeps = xr.open_dataset("spring_highLatkeps/OSMOSIS_glider_comparison.nc", decode_times=False)

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; time = gotmkpp.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; ketime = gotmkeps.time

time2 = 493.0 + (time[:]/86400); ketime2 = 493.0 + (ketime[:]/86400); ## Convert time
gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})

## AUTUMN low LAT ## with KPP local ##
#####################

gotmkppAut = xr.open_dataset("autumn_lowLatKPP/OSMOSIS_cruise_winter.nc", decode_times=False)
gotmkpplocalAut = xr.open_dataset("autumn_lowLatKPP/OSMOSIS_cruise_winter_local.nc", decode_times=False)
gotmkepsAut = xr.open_dataset("autumn_lowLatkeps/OSMOSIS_cruise_winter.nc", decode_times=False)

GTMtempAut = gotmkppAut.temp; GTMdepthAut = gotmkppAut.z; timeAut = gotmkppAut.time
GTMtempLocAut = gotmkpplocalAut.temp; GTMdepthLocAut = gotmkpplocalAut.z; timeLocAut = gotmkpplocalAut.time
GTMketempAut = gotmkepsAut.temp; GTMkedepthAut = gotmkepsAut.z; ketimeAut = gotmkepsAut.time

time2Aut = 301.875 + (timeAut[:]/86400); timeL2Aut = 301.875 + (timeLocAut[:]/86400); ketime2Aut = 301.875 + (ketimeAut[:]/86400); ## Convert time
gotmkppAut = gotmkppAut.assign(time2Aut=("time", time2Aut)); gotmkppAut = gotmkppAut.swap_dims({"time" : "time2Aut"})
gotmkpplocalAut = gotmkpplocalAut.assign(timeL2Aut=("time", timeL2Aut)); gotmkpplocalAut = gotmkpplocalAut.swap_dims({"time" : "timeL2Aut"})
# gotmkepsAut = gotmkepsAut.assign(ketime2=("time", ketime2Aut)); gotmkepsAut = gotmkepsAut.swap_dims({"time" : "ketime2Aut"})

## AUTUMN high LAT ##
#####################



## ANNUAL KPP LOCAL FLUXES ##
#############################

gotmkpplocal = xr.open_dataset("KPPlocal_annual/OSMOSIS_glider_comparison.nc", decode_times=False)

GTMtempL = gotmkpplocal.temp; GTMdepthL = gotmkpplocal.z; timeL = gotmkpplocal.time

localtime2 = 261.0 + (timeL[:]/86400); ## Convert time
gotmkpplocal = gotmkpplocal.assign(localtime2=("time", localtime2)); gotmkpplocal = gotmkpplocal.swap_dims({"time" : "localtime2"})

######---------------------------------------------------------------######
### -------- initiate (define) new variables as data elements --------- ###
######---------------------------------------------------------------######

class GliderlowPoints:
	def __init__(self, Gldrtemp, Gldrdepth, Gldrtime):
	    self.Gldrtemp = Gldrtemp
	    self.Gldrdepth = Gldrdepth
	    self.Gldrtime = Gldrtime

class kepslowPoints:
	def __init__(self, kepstemp1d, kepsdepth, kepstime):
	    self.kepstemp1d = kepstemp1d
	    self.kepsdepth = kepsdepth
	    self.kepstime = kepstime

class KPPlowPoints:
	def __init__(self, kpptemp1d, kppdepth, kpptime):
	    self.kpptemp1d = kpptemp1d
	    self.kppdepth = kppdepth
	    self.kpptime = kpptime

class kepsAutlowPoints:
	def __init__(self,  kepstempAut, kepsdepthAut, kepstimeAut):
		self.kepstempAut = kepstempAut
		self.kepsdepthAut = kepsdepthAut
		self.kepstimeAut = kepstimeAut

class KPPAutlowPoints:
	def __init__(self, kpptempAut, kppdepthAut, kpptimeAut):
		self.kpptempAut = kpptempAut
		self.kppdepthAut = kppdepthAut
		self.kpptimeAut = kpptimeAut

class KPPLocAutlowPoints:
	def __init__(self, kpptempLocAut, kppdepthLocAut, kpptimeLocAut):
		self.kpptempLocAut = kpptempLocAut
		self.kppdepthLocAut = kppdepthLocAut
		self.kpptimeLocAut = kpptimeLocAut

class CTtempPoints:
	def __init__(self, CTtemp, CTtime):
	    self.CTtemp = CTtemp
	    self.CTtime = CTtime
	    # self.CTdepth = CTdepth

class localKPPPoints:
	def __init__(self,kpptempL1d, kppdepthL, kpptimeL, kpptemp1dA, kppdepthA, kpptimeA):
		self.kpptempL1d = kpptempL1d
		self.kppdepthL = kppdepthL
		self.kpptimeL = kpptimeL
		self.kpptemp1dA = kpptemp1dA
		self.kppdepthA = kppdepthA
		self.kpptimeA = kpptimeA


##########################################################################
### -------- OSMOSIS Spring observations with low Lat number --------- ###
##########################################################################

# define a function to select the corresponding values of temperature, depth and times from the array.
def gldrSprlowLat(a):
	data = []
	for i in a:
		Gldrtemp = GldrTemp.isel()[i,:]; 
		Gldrdepth = GldrDepth[:];
		Gldrtime = new_time[i]
		data.append(GliderlowPoints(Gldrtemp, Gldrdepth, Gldrtime))
	return data

# define a function to plot the temperature against the depth
def visualise_gldrSprlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.Gldrtemp, point.Gldrdepth, label = label%point.Gldrtime.values); 

### -------- spring k-epsilon model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kepsSprlowLat(b):
	data = []
	for i in b:
		kepstemp1d = GTMketempSpr.isel(lon=0,lat=0)[i,:]; 
		kepsdepth = GTMkedepthSpr[:]
		kepstime = ketimeSpr[i]
		data.append(kepslowPoints(kepstemp1d, kepsdepth, kepstime))
	return data

# define a function to plot the temperature against the depth
def visualise_kepsSprlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kepstemp1d, point.kepsdepth, label = label%point.kepstime.values)

### -------- spring KPP model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kppSprlowLat(b):
	data = []
	for i in b:
		kpptemp1d = GTMtempSpr.isel(lon=0,lat=0)[i,:]; 
		kppdepth = GTMdepthSpr[:]
		kpptime = timeSpr[i]
		data.append(KPPlowPoints(kpptemp1d, kppdepth, kpptime))
	return data

# define a function to plot the temperature against the depth
def visulaise_kppSprlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kpptemp1d, point.kppdepth, label = label%point.kpptime.values)

##########################################################################
### -------- OSMOSIS Spring observations with high Lat number --------- ###
##########################################################################

# define a function to select the corresponding values of temperature, depth and times from the array.
def gldrSprhighLat(a):
	data = []
	for i in a:
		Gldrtemp = GldrTemp.isel()[i,:]; 
		Gldrdepth = GldrDepth[:];
		Gldrtime = new_time[i]
		data.append(GliderlowPoints(Gldrtemp, Gldrdepth, Gldrtime))
	return data

# define a function to plot the temperature against the depth
def visualise_gldrSprhighLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.Gldrtemp, point.Gldrdepth, label = label%point.Gldrtime.values); 

### -------- spring k-epsilon model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kepsSprhighLat(b):
	data = []
	for i in b:
		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; 
		kepsdepth = GTMkedepth[:]
		kepstime = time2[i]
		data.append(kepslowPoints(kepstemp1d, kepsdepth, kepstime))
	return data

# define a function to plot the temperature against the depth
def visualise_kepsSprhighLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kepstemp1d, point.kepsdepth, label = label%point.kepstime.values)

### -------- spring KPP model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kppSprhighLat(b):
	data = []
	for i in b:
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; 
		kppdepth = GTMdepth[:]
		kpptime = ketime2[i]
		data.append(KPPlowPoints(kpptemp1d, kppdepth, kpptime))
	return data

# define a function to plot the temperature against the depth
def visualise_kppSprhighLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kpptemp1d, point.kppdepth, label = label%point.kpptime.values)

##########################################################################
### -------- OSMOSIS Autumn observations with low Lat number --------- ###
##########################################################################

# define a function to select the corresponding values of temperature, depth and times from the array.
def gldrAutlowLat(a):
	data = []
	for i in a:
		Gldrtemp = GldrTemp.isel()[i,:]; 
		Gldrdepth = GldrDepth[:];
		Gldrtime = new_time[i]
		data.append(GliderlowPoints(Gldrtemp, Gldrdepth, Gldrtime))
	return data

# define a function to plot the temperature against the depth
def visualise_gldAutlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.Gldrtemp, point.Gldrdepth, label = label%point.Gldrtime.values); 

### -------- autumn k-epsilon model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kepsAutlowLat(b):
	data = []
	for i in b:
		kepstempAut = GTMketempAut.isel(lon=0,lat=0)[i,:]; 
		kepsdepthAut = GTMkedepthAut[:]
		kepstimeAut = ketime2Aut[i]
		data.append(kepsAutlowPoints(kepstempAut, kepsdepthAut, kepstimeAut))
	return data

# define a function to plot the temperature against the depth
def visualise_kepAutlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kepstimeAut, point.kepsdepthAut, label = label%point.kepstimeAut.values)

### -------- autumn KPP model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kppAutlowLat(b):
	data = []
	for i in b:
		kpptempAut = GTMtempAut.isel(lon=0,lat=0)[i,:]; 
		kppdepthAut = GTMdepthAut[:]
		kpptimeAut = time2Aut[i]
		data.append(KPPAutlowPoints(kpptempAut, kppdepthAut, kpptimeAut))
	return data

# define a function to plot the temperature against the depth
def visualise_kppAutlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kpptempAut, point.kppdepthAut, label = label%point.kpptimeAut.values)

### -------- autumn KPP (local) model simulation with low Lat number --------- ###

# define a function to select the corresponding values of temperature, depth and times from the array.
def kppLocAutlowLat(b):
	data = []
	for i in b:
		kpptempLocAut = GTMtempLocAut.isel(lon=0,lat=0)[i,:]; 
		kppdepthLocAut = GTMdepthLocAut[:]
		kpptimeLocAut = timeL2Aut[i]
		data.append(KPPLocAutlowPoints(kpptempLocAut, kppdepthLocAut, kpptimeLocAut))
	return data

# define a function to plot the temperature against the depth
def visualise_kppLocAutlowLat(plt, data, label = ''):
	for point in data:
		plt.plot(point.kpptempLocAut, point.kppdepthLocAut, label = label%point.kpptimeLocAut.values)

#############################################################
### -------- OSMOSIS Annual surface heat budget --------- ###
#############################################################
def ctsolver(dt,T):
	"""
	Use centered difference method to obtain the rate of change in time from temperature,
	in range of the mixed layer depth
	"""
	Nt = int(round(T/dt))
	t = np.linspace(261, Nt*dt, Nt+1) # Mesh points in time	
	# Nz = int(round(L/dz))
	# z = np.linspace(-999.0, L, Nz+1) # Mesh points in space
	# temp = []
	# for n in range(1, Nt):
	# 	for i in range(0, Nz):
	# 		CT = (GTMtemp.isel(lon=0,lat=0)[:,i+1]-GTMtemp.isel(lon=0,lat=0)[:,i])
	# 	CT2 = CT[n]/dt
	# 	temp.append(CT2)
	temp = []
	time = []
	# for i in range(1,Nz):
	# 	sumTempMLD[i] = np.sum(GTMtempA.isel(lon=0,lat=0)[i,400:500]) 
	for n in range(1, Nt):
		CT = (sumkppTempMLD[n+1]-sumkppTempMLD[n])/dt
		# CT = (GTMtempA.isel(lon=0,lat=0)[n+1,:]-GTMtempA.isel(lon=0,lat=0)[n,:])/dt
		diffT = (kpptime2[n+1]+kpptime2[n])/2
		temp.append(CT)
		time.append(diffT)
	return temp, time, t

def ctvariables(count,temp,time):
	CTdata = []
	for i in count:
		CTtemp = temp[i]
		CTtime = time[i]
		CTdata.append(CTtempPoints(CTtemp,CTtime))
	return CTdata


def visualize(plt, CTdata):
	for point in CTdata:
		plt.plot(point.CTtime,point.CTtemp)
	# plt.xlabel('time')
	# plt.ylabel('temperature -C')
	plt.savefig('heatFlux_osmosisAnnual/GOTM_surface_heat_budget.png')
	plt.show()


##########################################################
### -------- OSMOSIS Annual local KPP fluxes --------- ###
##########################################################


# define a function to select the corresponding values of temperature, depth and times from the array.
def compareKPP(Gldr_sel,GTM_sel):
	kppdata = []
	kppdataL = []
	for j in Gldr_sel:
		Gldrtemp = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
		Gldrtime = new_time[j]
		kppdata.append(GliderlowPoints(Gldrtemp,Gldrdepth,Gldrtime))
		for i in GTM_sel:
			kpptemp1dA = GTMtempA.isel(lon=0,lat=0)[i,:]; 
			kppdepthA = GTMdepthA[:]
			kpptempL1d = GTMtempL.isel(lon=0,lat=0)[i,:]; 
			kppdepthL = GTMdepthL[:]
			kpptimeA = kpptime2[i]
			kpptimeL = localtime2[i]
			kppdataL.append(localKPPPoints(kpptempL1d, kppdepthL, kpptimeL, kpptemp1dA, kppdepthA, kpptimeA))
	return kppdata, kppdataL

# define a function to plot the temperature against the depth
def localKPPvisualize(plt, kppdata, kppdataL, label = ''):
	for point in kppdata:
		plt.plot(point.Gldrtemp, point.Gldrdepth, label = label%point.Gldrtime.values);
		for points in kppdataL:
			plt.plot(points.kpptemp1dA,points.kppdepthA,label = 'GOTM_KPP day %s'%points.kpptimeA.values); plt.legend(loc=(0.5,0.05), fontsize='x-small')
			plt.plot(points.kpptempL1d,points.kppdepthL,label = 'GOTM_KPP local day %s'%points.kpptimeL.values); plt.legend(loc=(0.5,0.05), fontsize='x-small')
	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
	# plt.title('Compare Temperatures : Glider & GOTM w KPP,KPP-local')
	plt.axis([10,21,-400,0])
	plt.savefig('KPPlocal_annual/Gldr_GOTM_kpp_Temp_profile_Nonlocal.png')
	plt.show()


##########################################################
### -------- OSMOSIS Annual local KPP fluxes --------- ###
##########################################################

