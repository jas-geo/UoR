## ------------- Written code to examine the periods with low h/L_L ------------- ##
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


##  ---------------------------  obtain initial conditions  --------------------------- ##
# print('Glider depth', GldrDepth.values)
print('Glider fixed time :', new_time[3243:4700])
# temperature
# print('Glider temp for time=464', GldrTemp[2340,:].values, len(GldrTemp[2340,:]))
Gldrtemp464 = GldrTemp[3243,:].values
t_profile464 = np.c_[GldrDepth,GldrTemp[3243,:].values]
print('Glider temp for time=464', t_profile464)
# salinity
# print('Glider salt for time=464', GldrSalt[2340,:].values, len(GldrSalt[2340,:]))
GldrSalt464 = GldrSalt[3243,:].values
s_profile464 = np.c_[GldrDepth,GldrSalt[3243,:].values]
print('Glider salt for time=464', s_profile464)
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

### ----------------------------- GOTM ----------------------------- ###


gotmkpp = xr.open_dataset("lowhL_KPP464/OSMOSIS_glider_KPP.nc", decode_times=False)
gotmkeps = xr.open_dataset("lowhL_KPP464/OSMOSIS_glider_keps.nc", decode_times=False)


GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMdepth1 = gotmkpp.z1; GTMu = gotmkpp.u; GTMv = gotmkpp.v; time = gotmkpp.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkedepth1 = gotmkeps.z1; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time

time2 = 464.0 + (time[:]/86400); ketime2 = 464.0 + (ketime[:]/86400)

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})

print('shape of Gldr temp', GldrTemp.shape)
print('shape of GTM temp', GTMtemp.shape)


## --------------------------- PLOTS- temperature --------------------------- ###

# ## -------------------------- Glider -------------------------- ##
# #choose the indices for each day
# # a = np.arange(482,554,10)
# a = [2341,2351,2363,2375,2387,2399,2410,2421]

# for i in a:
# 	Gldrtemp1d = GldrTemp.isel()[i,:]; Gldrdepth = GldrDepth[:]
# 	Gldrtime = new_time[i] 
# 	plt.plot(Gldrtemp1d,Gldrdepth,label='Glider Day %s'%Gldrtime.values);# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.legend(loc='best',fontsize ='small')
# 	plt.grid()
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	# plt.title('Glider Temperature Profiles for Spring high La_t values')
# 	# plt.axis([11.0,15.5,-150,0])
# 	plt.ylim([-400,0])
# plt.savefig('lowhL_KPP464/Glider_Temp_profile_lowhL.png')
# plt.show()

# ## -------------------------- KPP -------------------------- ##
# #choose the indices for each day
# # b = np.arange(0,1115,100)
# b = [0,144,288,432,576,720,864,1008,1151]


# for i in b:
# 	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
# 	kpptime = time2[i]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_KPP Day %s'%kpptime.values);# plt.legend(loc=(0.6,0.02),fontsize='xx-small')
# 	plt.legend(loc='best',fontsize='small')
# 	plt.grid()
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	# plt.title('GOTM KPP Temperature profiles for Spring high La_t values')
# 	# plt.axis([11.5,15.5,-150,0])
# 	plt.ylim([-400,0])
# plt.savefig('lowhL_KPP464/GOTM_kpp_Temp_profile_lowhL.png')
# plt.show()

# ## -------------------------- k-eps -------------------------- ##

# for i in b:
# 	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
# 	kepstime = ketime2[i]
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps Day %s'%kepstime.values);# plt.legend(loc=(0.6,0.02),fontsize='xx-small')
# 	plt.legend(loc='best',fontsize='small')
# 	plt.grid()
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	# plt.title('GOTM k-epsilon Temperature profiles for Spring high La_t values')
# 	# plt.axis([11.5,15.5,-150,0])
# 	plt.ylim([-400,0])
# plt.savefig('lowhL_KPP464/GOTM_keps_Temp_profile_lowhL.png')
# plt.show()

## --------------------------- PLOTS- velocities --------------------------- ###
# c = [499]

# for i in c:
# 	kppu1d = GTMu.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
# 	kppuLa = GTMuLa.isel(lon=0,lat=0)[:,i]; kppdepthLa = GTMdepthLa[i]
# 	kppu_heat = GTMu_heat.isel(lon=0,lat=0)[:,i]; kppdepth_heat = GTMdepth_heat[i]
# 	kppu_mh = GTMu_mh.isel(lon=0,lat=0)[:,i]; kppdepth_mh = GTMdepth_mh[i]
# 	kepsu1d = GTMkeu.isel(lon=0,lat=0)[:,i]; kepsdepth = GTMkedepth[i]
# 	kpptime = time2[:] 
# 	kpptimeLa = time2La[:]
# 	kpptime_heat = time2_heat[:]
# 	kpptime_mh = time2_mh[:]
# 	kepstime = ketime2[:]
# 	plt.plot(kpptime,kppu1d,'k-',label='KPP');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptimeLa,kppuLa,'b-',label='KPP EV');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_heat,kppu_heat,'g-',label='KPP HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_mh,kppu_mh,'r--',label='KPP EV-HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kepstime,kepsu1d,'k:',label='k-eps');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.legend(loc='best',fontsize ='small')
# 	plt.grid()
# 	plt.xlabel('Time -Julian days');plt.ylabel('Surface streamwise velocity -m/s'); 
# 	# plt.axis([11.0,15.5,-150,0])
# plt.savefig('osmosis_spring_variablePlots/GOTM_streamwise_velocity_timeseries_SpringhighLat.png')
# plt.show()

# for i in c:
# 	kppv1d = GTMv.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
# 	kppvLa = GTMvLa.isel(lon=0,lat=0)[:,i]; kppdepthLa = GTMdepthLa[i]
# 	kppv_heat = GTMv_heat.isel(lon=0,lat=0)[:,i]; kppdepth_heat = GTMdepth_heat[i]
# 	kppv_mh = GTMv_mh.isel(lon=0,lat=0)[:,i]; kppdepth_mh = GTMdepth_mh[i]
# 	kepsv1d = GTMkev.isel(lon=0,lat=0)[:,i]; kepsdepth = GTMkedepth[i]
# 	kpptime = time2[:] 
# 	kpptimeLa = time2La[:]
# 	kpptime_heat = time2_heat[:]
# 	kpptime_mh = time2_mh[:]
# 	kepstime = ketime2[:]
# 	plt.plot(kpptime,kppv1d,'k-',label='KPP');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptimeLa,kppvLa,'b-',label='KPP EV');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_heat,kppv_heat,'g-',label='KPP HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_mh,kppv_mh,'r--',label='KPP EV-HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kepstime,kepsv1d,'k:',label='k-eps');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.legend(loc='best',fontsize ='small')
# 	plt.grid()
# 	plt.xlabel('Time -Julian days');plt.ylabel('Surface spanwise velocity -m/s'); 
# 	# plt.axis([11.0,15.5,-150,0])
# plt.savefig('osmosis_spring_variablePlots/GOTM_spanwise_velocity_timeseries_SpringhighLat.png')
# plt.show()

## -------------------------- Comparisons ------------------------------ ##

GTM_sel = [0] # range [f:m:l] = [0:500:1151] # 590
Gldr_sel = [2340] # range [f:m:l] = [2341:2381:2421] # 3130

# Compared plot GLider and GOTM_KPP
for j in Gldr_sel:
	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
	Gldrtime = new_time[j] 
	for i in GTM_sel:
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
		# kpptempLoc = GTMtempLoc.isel(lon=0,lat=0)[i,:]; kppdepthLoc = GTMdepthLoc[:]
		# kpptempLa = GTMtempLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
		# kpptempLaLoc = GTMtempLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
		# kpptempLag05 = GTMtempLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
		# kpptempLag2 = GTMtempLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
		# kpptempLa05dlta = GTMtempLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
		# kpptempLa025dlta = GTMtempLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
		# kpptempLa125dlta = GTMtempLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
		# kpptempLa0625dlta = GTMtempLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
		# kpptemp_heat = GTMtemp_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
		# kpptemp_heatLoc = GTMtemp_heatLoc.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc = GTMdepth_heatLoc[:]
		# kpptemp_heatg05 = GTMtemp_heatg05.isel(lon=0,lat=0)[i,:]; kppdepth_heatg05 = GTMdepth_heatg05[:]
		# kpptemp_heatg2 = GTMtemp_heatg2.isel(lon=0,lat=0)[i,:]; kppdepth_heatg2 = GTMdepth_heatg2[:]
		# kpptemp_heat05dlta = GTMtemp_heat05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat05dlta = GTMdepth_heat05dlta[:]
		# kpptemp_heat025dlta = GTMtemp_heat025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat025dlta = GTMdepth_heat025dlta[:]
		# kpptemp_heat125dlta = GTMtemp_heat125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat125dlta = GTMdepth_heat125dlta[:]
		# kpptemp_heat0625dlta = GTMtemp_heat0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat0625dlta = GTMdepth_heat0625dlta[:]
		# kpptemp_mh = GTMtemp_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
		# kpptemp_mhLoc = GTMtemp_mhLoc.isel(lon=0,lat=0)[i,:]; kppdepth_mhLoc = GTMdepth_mhLoc[:]
		# kpptemp_mhg05 = GTMtemp_mhg05.isel(lon=0,lat=0)[i,:]; kppdepth_mhg05 = GTMdepth_mhg05[:]
		# kpptemp_mhg2 = GTMtemp_mhg2.isel(lon=0,lat=0)[i,:]; kppdepth_mhg2 = GTMdepth_mhg2[:]
		# kpptemp_mh05dlta = GTMtemp_mh05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh05dlta = GTMdepth_mh05dlta[:]
		# kpptemp_mh025dlta = GTMtemp_mh025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh025dlta = GTMdepth_mh025dlta[:]
		# kpptemp_mh125dlta = GTMtemp_mh125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh125dlta = GTMdepth_mh125dlta[:]
		# kpptemp_mh0625dlta = GTMtemp_mh0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh0625dlta = GTMdepth_mh0625dlta[:]
		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
		kpptime = time2[i]
		# kpptimeLoc = time2Loc[i]
		# kpptimeLa = time2La[i]
		# kpptimeLaLoc = time2LaLoc[i]
		# kpptimeLag05 = time2Lag05[i]
		# kpptimeLag2 = time2Lag2[i]
		# kpptimeLa05dlta = time2La05dlta[i]
		# kpptimeLa025dlta = time2La025dlta[i]
		# kpptimeLa125dlta = time2La125dlta[i]
		# kpptimeLa0625dlta = time2La0625dlta[i]
		# kpptime_heat = time2_heat[i]
		# kpptime_heatLoc = time2_heatLoc[i]
		# kpptime_heatg05 = time2_heatg05[i]
		# kpptime_heatg2 = time2_heatg2[i]
		# kpptime_heat05dlta = time2_heat05dlta[i]
		# kpptime_heat025dlta = time2_heat025dlta[i]
		# kpptime_heat125dlta = time2_heat125dlta[i]		
		# kpptime_heat0625dlta = time2_heat0625dlta[i]		
		# kpptime_mh = time2_mh[i]
		# kpptime_mhLoc = time2_mhLoc[i]
		# kpptime_mhg05 = time2_mhg05[i]
		# kpptime_mhg2 = time2_mhg2[i]
		# kpptime_mh05dlta = time2_mh05dlta[i]
		# kpptime_mh025dlta = time2_mh025dlta[i]
		# kpptime_mh125dlta = time2_mh125dlta[i]
		# kpptime_mh0625dlta = time2_mh0625dlta[i]
		kepstime = ketime2[i]
	plt.plot(kpptemp1d,kppdepth,color='gray',linestyle='solid',label = 'KPP d=%0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLoc,kppdepthLoc,color='gray',linestyle='dashed',label = 'KPP local day %0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLaLoc,kppdepthLaLoc,color='steelblue',linestyle='dashed',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heat,kppdepth_heat,color='firebrick',linestyle='solid',label = 'KPP HD d=%0.2f'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heatLoc,kppdepth_heatLoc,color='firebrick',linestyle='dashed',label = 'KPP HD loc d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heatg05,kppdepth_heatg05,color='firebrick',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_heatg05.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heatg2,kppdepth_heatg2,color='firebrick',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_heatg2.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heat05dlta,kppdepth_heat05dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_heat05dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heat025dlta,kppdepth_heat025dlta,color='firebrick',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_heat025dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heat125dlta,kppdepth_heat125dlta,color='firebrick',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_heat125dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heat0625dlta,kppdepth_heat0625dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_heat0625dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mh,kppdepth_mh,color='darkorchid',linestyle='solid',label = 'KPP EV-HD d=%0.2f'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
	# # plt.plot(kpptemp_mhLoc,kppdepth_mhLoc,color='darkorchid',linestyle='dashed',label = 'KPP EV-HD loc d=%0.2f'%kpptime_mhLoc.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mhg05,kppdepth_mhg05,color='darkorchid',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_mhg05.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mhg2,kppdepth_mhg2,color='darkorchid',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_mhg2.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mh05dlta,kppdepth_mh05dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_mh05dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mh025dlta,kppdepth_mh025dlta,color='darkorchid',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_mh025dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mh125dlta,kppdepth_mh125dlta,color='darkorchid',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_mh125dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_mh0625dlta,kppdepth_mh0625dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_mh0625dlta.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kepstemp1d,kepsdepth,color='goldenrod',linestyle='solid',label = 'k-eps d=%0.2f'%kepstime.values); plt.legend(loc='best', fontsize='small')
	plt.plot(Gldrtemp1d,Gldrdepth,color='olivedrab',linestyle='solid',label='OSMOSIS d=%0.2f'%Gldrtime.values); plt.legend(loc='best', fontsize='small')
	plt.grid()
	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
	# plt.title('Compared Temperature Profiles : Glider and GOTM w KPP,k-eps')
	# plt.axis([11,15.5,-150,0])
	plt.ylim(-400,0)
plt.savefig('lowhL_KPP464/Gldr_GOTM_Temperature_profile_lowhL.png')
plt.show()
exit()
# for i in GTM_sel:
# 	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
# 	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
# 	kpptime = time2[i]
# 	kepstime = ketime2[i]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_KPP Day %s'%kpptime.values); plt.legend(loc=(0.45,0.35))
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps Day %s'%kepstime.values); plt.legend(loc=(0.45,0.35))
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	plt.title('Compared GOTM Temperature Profiles : KPP and k-eps models')
# 	plt.axis([11,18,-150,0])
# plt.savefig('GOTM_kppkeps_Temp_profile_lowLat_test.png')
# plt.show()

## EDDY VISCOSITY
for i in GTM_sel:
	kppnum1d = GTMnum.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	# kppnumLa = GTMnumLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppnumLaLoc = GTMnumLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	# kppnumLag05 = GTMnumLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	# kppnumLag2 = GTMnumLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	# kppnumLa05dlta = GTMnumLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	# kppnumLa025dlta = GTMnumLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	# kppnumLa125dlta = GTMnumLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	# kppnumLa0625dlta = GTMnumLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
	# kppnum_heat = GTMnum_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
	# kppnum_heatLoc = GTMnum_heatLoc.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc = GTMdepth_heatLoc[:]
	# kppnum_heatg05 = GTMnum_heatg05.isel(lon=0,lat=0)[i,:]; kppdepth_heatg05 = GTMdepth_heatg05[:]
	# kppnum_heatg2 = GTMnum_heatg2.isel(lon=0,lat=0)[i,:]; kppdepth_heatg2 = GTMdepth_heatg2[:]
	# kppnum_heat05dlta = GTMnum_heat05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat05dlta = GTMdepth_heat05dlta[:]
	# kppnum_heat025dlta = GTMnum_heat025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat025dlta = GTMdepth_heat025dlta[:]
	# kppnum_heat125dlta = GTMnum_heat125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat125dlta = GTMdepth_heat125dlta[:]
	# kppnum_heat0625dlta = GTMnum_heat0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat0625dlta = GTMdepth_heat0625dlta[:]
	# kppnum_mh = GTMnum_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
	# kppnum_mhLoc = GTMnum_mhLoc.isel(lon=0,lat=0)[i,:]; kppdepth_mhLoc = GTMdepth_mhLoc[:]
	# kppnum_mhg05 = GTMnum_mhg05.isel(lon=0,lat=0)[i,:]; kppdepth_mhg05 = GTMdepth_mhg05[:]
	# kppnum_mhg2 = GTMnum_mhg2.isel(lon=0,lat=0)[i,:]; kppdepth_mhg2 = GTMdepth_mhg2[:]
	# kppnum_mh05dlta = GTMnum_mh05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh05dlta = GTMdepth_mh05dlta[:]
	# kppnum_mh025dlta = GTMnum_mh025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh025dlta = GTMdepth_mh025dlta[:]
	# kppnum_mh125dlta = GTMnum_mh125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh125dlta = GTMdepth_mh125dlta[:]
	# kppnum_mh0625dlta = GTMnum_mh0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh0625dlta = GTMdepth_mh0625dlta[:]
	kepsnum1d = GTMkenum.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kpptime = time2[i]
	# kpptimeLa = time2La[i]
	# kpptimeLaLoc = time2LaLoc[i]
	# kpptimeLag05 = time2Lag05[i]
	# kpptimeLag2 = time2Lag2[i]
	# kpptimeLa05dlta = time2La05dlta[i]
	# kpptimeLa025dlta = time2La025dlta[i]
	# kpptimeLa125dlta = time2La125dlta[i]
	# kpptimeLa0625dlta = time2La0625dlta[i]
	# kpptime1 = time21[i]
	# kpptime_heat = time2_heat[i]
	# kpptime_heatLoc = time2_heatLoc[i]
	# kpptime_heatg05 = time2_heatg05[i]
	# kpptime_heatg2 = time2_heatg2[i]
	# kpptime_heat05dlta = time2_heat05dlta[i]
	# kpptime_heat025dlta = time2_heat025dlta[i]
	# kpptime_heat125dlta = time2_heat125dlta[i]		
	# kpptime_heat0625dlta = time2_heat0625dlta[i]		
	# kpptime_mh = time2_mh[i]
	# kpptime_mhLoc = time2_mhLoc[i]
	# kpptime_mhg05 = time2_mhg05[i]
	# kpptime_mhg2 = time2_mhg2[i]
	# kpptime_mh05dlta = time2_mh05dlta[i]
	# kpptime_mh025dlta = time2_mh025dlta[i]
	# kpptime_mh125dlta = time2_mh125dlta[i]
	# kpptime_mh0625dlta = time2_mh0625dlta[i]
	kepstime = ketime2[i]
plt.plot(kppnum1d,kppdepth,color='gray',linestyle='solid',label = 'KPP d=%0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum1,kppdepth1,'c-.',label = 'KPP1 d=%0.2f'%kpptime1.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLa,kppdepthLa,color='steelblue',linestyle='solid' ,label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnumLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heat,kppdepth_heat,color='firebrick',linestyle='solid',label = 'KPP HD d=%0.2f'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnum_heatLoc,kppdepth_heatLoc,color='firebrick',linestyle='dashed',label = 'KPP HD loc d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnum_heatLoc1,kppdepth_heatLoc1,color='grey',label = 'KPP HD loc1 d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heatg2,kppdepth_heatg2,color='firebrick',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_heatg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heatg05,kppdepth_heatg05,color='firebrick',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_heatg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heat05dlta,kppdepth_heat05dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_heat05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heat025dlta,kppdepth_heat025dlta,color='firebrick',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_heat025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heat125dlta,kppdepth_heat125dlta,color='firebrick',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_heat125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heat0625dlta,kppdepth_heat0625dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_heat0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mh,kppdepth_mh,color='darkorchid',linestyle='solid',label = 'KPP EV-HD d=%0.2f'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnum_mhLoc,kppdepth_mhLoc,color='darkorchid',linestyle='dashed',label = 'KPP EV-HD loc d=%0.2f'%kpptime_mhLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mhg05,kppdepth_mhg05,color='darkorchid',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_mhg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mhg2,kppdepth_mhg2,color='darkorchid',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_mhg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mh05dlta,kppdepth_mh05dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_mh05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mh025dlta,kppdepth_mh025dlta,color='darkorchid',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_mh025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mh125dlta,kppdepth_mh125dlta,color='darkorchid',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_mh125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mh0625dlta,kppdepth_mh0625dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_mh0625dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnum1d,kepsdepth,color='goldenrod',linestyle='solid' ,label = 'k-eps d=%0.2f'%kepstime.values); plt.legend(loc='best', fontsize='small')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.grid()
plt.xlabel('Eddy Viscosity -m2/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-150,0])
plt.ylim(-400,0)
plt.savefig('lowhL_KPP464/GOTM_eddy_viscosity_profile_lowhL.png')
plt.show()

## HEAT DIFFUSIVITY
for i in GTM_sel:
	kppnuh1d = GTMnuh.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	# kppnuh1 = GTMnuh1.isel(lon=0,lat=0)[i,:]; kppdepth1 = GTMdepth1[:]
	# kppnuhLa = GTMnuhLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppnuhLaLoc = GTMnuhLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	# kppnuhLag05 = GTMnuhLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	# kppnuhLag2 = GTMnuhLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	# kppnuhLa05dlta = GTMnuhLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	# kppnuhLa025dlta = GTMnuhLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	# kppnuhLa125dlta = GTMnuhLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	# kppnuhLa0625dlta = GTMnuhLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
	# kppnuh_heat = GTMnuh_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
	# # kppnuh_heatLoc = GTMnuh_heatLoc.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc = GTMdepth_heatLoc[:]
	# # kppnuh_heatLoc1 = GTMnuh_heatLoc1.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc1 = GTMdepth_heatLoc1[:]
	# kppnuh_heatg05 = GTMnuh_heatg05.isel(lon=0,lat=0)[i,:]; kppdepth_heatg05 = GTMdepth_heatg05[:]
	# kppnuh_heatg2 = GTMnuh_heatg2.isel(lon=0,lat=0)[i,:]; kppdepth_heatg2 = GTMdepth_heatg2[:]
	# kppnuh_heat05dlta = GTMnuh_heat05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat05dlta = GTMdepth_heat05dlta[:]
	# kppnuh_heat025dlta = GTMnuh_heat025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat025dlta = GTMdepth_heat025dlta[:]
	# kppnuh_heat125dlta = GTMnuh_heat125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat125dlta = GTMdepth_heat125dlta[:]
	# kppnuh_heat0625dlta = GTMnuh_heat0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat0625dlta = GTMdepth_heat0625dlta[:]
	# kppnuh_mh = GTMnuh_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
	# kppnuh_mhLoc = GTMnuh_mhLoc.isel(lon=0,lat=0)[i,:]; kppdepth_mhLoc = GTMdepth_mhLoc[:]
	# kppnuh_mhg05 = GTMnuh_mhg05.isel(lon=0,lat=0)[i,:]; kppdepth_mhg05 = GTMdepth_mhg05[:]
	# kppnuh_mhg2 = GTMnuh_mhg2.isel(lon=0,lat=0)[i,:]; kppdepth_mhg2 = GTMdepth_mhg2[:]
	# kppnuh_mh05dlta = GTMnuh_mh05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh05dlta = GTMdepth_mh05dlta[:]
	# kppnuh_mh025dlta = GTMnuh_mh025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh025dlta = GTMdepth_mh025dlta[:]
	# kppnuh_mh125dlta = GTMnuh_mh125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh125dlta = GTMdepth_mh125dlta[:]
	# kppnuh_mh0625dlta = GTMnuh_mh0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh0625dlta = GTMdepth_mh0625dlta[:]
	kepsnuh1d = GTMkenuh.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kpptime = time2[i]
	# kpptimeLa = time2La[i]
	# kpptimeLag05 = time2Lag05[i]
	# kpptimeLag2 = time2Lag2[i]
	# kpptimeLa05dlta = time2La05dlta[i]
	# kpptimeLa025dlta = time2La025dlta[i]
	# kpptimeLa125dlta = time2La125dlta[i]
	# kpptimeLa0625dlta = time2La0625dlta[i]
	# kpptime1 = time21[i]
	# kpptime_heat = time2_heat[i]
	# kpptime_heatLoc = time2_heatLoc[i]
	# kpptime_heatg05 = time2_heatg05[i]
	# kpptime_heatg2 = time2_heatg2[i]
	# kpptime_heat05dlta = time2_heat05dlta[i]
	# kpptime_heat025dlta = time2_heat025dlta[i]
	# kpptime_heat125dlta = time2_heat125dlta[i]		
	# kpptime_heat0625dlta = time2_heat0625dlta[i]		
	# kpptime_mh = time2_mh[i]
	# kpptime_mhLoc = time2_mhLoc[i]
	# kpptime_mhg05 = time2_mhg05[i]
	# kpptime_mhg2 = time2_mhg2[i]
	# kpptime_mh05dlta = time2_mh05dlta[i]
	# kpptime_mh025dlta = time2_mh025dlta[i]
	# kpptime_mh125dlta = time2_mh125dlta[i]
	# kpptime_mh0625dlta = time2_mh0625dlta[i]
	kepstime = ketime2[i]
plt.plot(kppnuh1d,kppdepth,color='gray',linestyle='solid',label = 'KPP d=%0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh1,kppdepth1,'c-.',label = 'KPP1 d=%0.2f'%kpptime1.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuhLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heat,kppdepth_heat,color='firebrick',linestyle='solid',label = 'KPP HD d=%0.2f'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_heatLoc,kppdepth_heatLoc,color='olive',linestyle='dashed',label = 'KPP HD loc d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatLoc1,kppdepth_heatLoc1,color='grey',label = 'KPP HD loc1 d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg2,kppdepth_heatg2,color='firebrick',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_heatg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg05,kppdepth_heatg05,color='firebrick',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_heatg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heat05dlta,kppdepth_heat05dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_heat05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heat025dlta,kppdepth_heat025dlta,color='firebrick',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_heat025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heat125dlta,kppdepth_heat125dlta,color='firebrick',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_heat125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heat0625dlta,kppdepth_heat0625dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_heat0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh,kppdepth_mh,color='darkorchid',linestyle='solid',label = 'KPP EV-HD d=%0.2f'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_mhLoc,kppdepth_mhLoc,color='darkorchid',linestyle='dashed',label = 'KPP EV-HD loc d=%0.2f'%kpptime_mhLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mhg05,kppdepth_mhg05,color='darkorchid',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_mhg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mhg2,kppdepth_mhg2,color='darkorchid',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_mhg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh05dlta,kppdepth_mh05dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_mh05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh025dlta,kppdepth_mh025dlta,color='darkorchid',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_mh025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh125dlta,kppdepth_mh125dlta,color='darkorchid',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_mh125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh0625dlta,kppdepth_mh0625dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_mh0625dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnuh1d,kepsdepth,color='goldenrod',linestyle='solid' ,label = 'k-eps d=%0.2f'%kepstime.values); plt.legend(loc='best', fontsize='small')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.grid()
plt.xlabel('Eddy Diffusivity -m2/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-150,0])
plt.ylim(-400,0)
plt.savefig('lowhL_KPP464/GOTM_heat_diffusivity_profile_lowhL.png')
plt.show()


## Velocity profiles

# U
for i in GTM_sel:
	kppu1d = GTMu.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	# kppuLa = GTMuLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppuLaLoc = GTMuLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	# kppuLag05 = GTMuLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	# kppuLag2 = GTMuLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	# kppuLa05dlta = GTMuLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	# kppuLa025dlta = GTMuLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	# kppuLa125dlta = GTMuLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	# kppuLa0625dlta = GTMuLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
	# kppu_heat = GTMu_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
	# # kppu_heatLoc = GTMu_heatLoc.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc = GTMdepth_heatLoc[:]
	# # kppu_heatLoc1 = GTMu_heatLoc1.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc1 = GTMdepth_heatLoc1[:]
	# kppu_heatg05 = GTMu_heatg05.isel(lon=0,lat=0)[i,:]; kppdepth_heatg05 = GTMdepth_heatg05[:]
	# kppu_heatg2 = GTMu_heatg2.isel(lon=0,lat=0)[i,:]; kppdepth_heatg2 = GTMdepth_heatg2[:]
	# kppu_heat05dlta = GTMu_heat05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat05dlta = GTMdepth_heat05dlta[:]
	# kppu_heat025dlta = GTMu_heat025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat025dlta = GTMdepth_heat025dlta[:]
	# kppu_heat125dlta = GTMu_heat125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat125dlta = GTMdepth_heat125dlta[:]
	# kppu_heat0625dlta = GTMu_heat0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat0625dlta = GTMdepth_heat0625dlta[:]
	# kppu_mh = GTMu_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
	# kppu_mhLoc = GTMu_mhLoc.isel(lon=0,lat=0)[i,:]; kppdepth_mhLoc = GTMdepth_mhLoc[:]
	# kppu_mhg05 = GTMu_mhg05.isel(lon=0,lat=0)[i,:]; kppdepth_mhg05 = GTMdepth_mhg05[:]
	# kppu_mhg2 = GTMu_mhg2.isel(lon=0,lat=0)[i,:]; kppdepth_mhg2 = GTMdepth_mhg2[:]
	# kppu_mh05dlta = GTMu_mh05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh05dlta = GTMdepth_mh05dlta[:]
	# kppu_mh025dlta = GTMu_mh025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh025dlta = GTMdepth_mh025dlta[:]
	# kppu_mh125dlta = GTMu_mh125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh125dlta = GTMdepth_mh125dlta[:]
	# kppu_mh0625dlta = GTMu_mh0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh0625dlta = GTMdepth_mh0625dlta[:]
	kepsu1d = GTMkeu.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kpptime = time2[i]
	# kpptimeLa = time2La[i]
	# kpptimeLag05 = time2Lag05[i]
	# kpptimeLag2 = time2Lag2[i]
	# kpptimeLa05dlta = time2La05dlta[i]
	# kpptimeLa025dlta = time2La025dlta[i]
	# kpptimeLa125dlta = time2La125dlta[i]
	# kpptimeLa0625dlta = time2La0625dlta[i]
	# kpptime_heat = time2_heat[i]
	# kpptime_heatLoc = time2_heatLoc[i]
	# kpptime_heatg05 = time2_heatg05[i]
	# kpptime_heatg2 = time2_heatg2[i]
	# kpptime_heat05dlta = time2_heat05dlta[i]
	# kpptime_heat025dlta = time2_heat025dlta[i]
	# kpptime_heat125dlta = time2_heat125dlta[i]		
	# kpptime_heat0625dlta = time2_heat0625dlta[i]		
	# kpptime_mh = time2_mh[i]
	# kpptime_mhLoc = time2_mhLoc[i]
	# kpptime_mhg05 = time2_mhg05[i]
	# kpptime_mhg2 = time2_mhg2[i]
	# kpptime_mh05dlta = time2_mh05dlta[i]
	# kpptime_mh025dlta = time2_mh025dlta[i]
	# kpptime_mh125dlta = time2_mh125dlta[i]
	# kpptime_mh0625dlta = time2_mh0625dlta[i]
	kepstime = ketime2[i]
plt.plot(kppu1d,kppdepth,color='gray',linestyle='solid',label = 'KPP d=%0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppuLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heat,kppdepth_heat,color='firebrick',linestyle='solid',label = 'KPP HD d=%0.2f'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppu_heatLoc,kppdepth_heatLoc,color='olive',linestyle='dashed',label = 'KPP HD loc d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppu_heatLoc1,kppdepth_heatLoc1,color='grey',label = 'KPP HD loc1 d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heatg2,kppdepth_heatg2,color='firebrick',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_heatg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heatg05,kppdepth_heatg05,color='firebrick',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_heatg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heat05dlta,kppdepth_heat05dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_heat05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heat025dlta,kppdepth_heat025dlta,color='firebrick',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_heat025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heat125dlta,kppdepth_heat125dlta,color='firebrick',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_heat125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heat0625dlta,kppdepth_heat0625dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_heat0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mh,kppdepth_mh,color='darkorchid',linestyle='solid',label = 'KPP EV-HD d=%0.2f'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppu_mhLoc,kppdepth_mhLoc,color='darkorchid',linestyle='dashed',label = 'KPP EV-HD loc d=%0.2f'%kpptime_mhLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mhg05,kppdepth_mhg05,color='darkorchid',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_mhg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mhg2,kppdepth_mhg2,color='darkorchid',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_mhg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mh05dlta,kppdepth_mh05dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_mh05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mh025dlta,kppdepth_mh025dlta,color='darkorchid',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_mh025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mh125dlta,kppdepth_mh125dlta,color='darkorchid',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_mh125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mh0625dlta,kppdepth_mh0625dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_mh0625dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsu1d,kepsdepth,color='goldenrod',linestyle='solid' ,label = 'k-eps d=%0.2f'%kepstime.values); plt.legend(loc='best', fontsize='small')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.grid()
plt.xlabel('Streamwise velocity -m/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-150,0])
plt.ylim(-400,0)
plt.savefig('lowhL_KPP464/GOTM_streamwise_velocity_profile_lowhL.png')
plt.show()

# V
for i in GTM_sel:
	kppv1d = GTMv.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	# kppvLa = GTMvLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppuLaLoc = GTMvLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	# kppvLag05 = GTMvLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	# kppvLag2 = GTMvLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	# kppvLa05dlta = GTMvLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	# kppvLa025dlta = GTMvLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	# kppvLa125dlta = GTMvLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	# kppvLa0625dlta = GTMvLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
	# kppv_heat = GTMv_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
	# # kppv_heatLoc = GTMv_heatLoc.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc = GTMdepth_heatLoc[:]
	# # kppv_heatLoc1 = GTMv_heatLoc1.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc1 = GTMdepth_heatLoc1[:]
	# kppv_heatg05 = GTMv_heatg05.isel(lon=0,lat=0)[i,:]; kppdepth_heatg05 = GTMdepth_heatg05[:]
	# kppv_heatg2 = GTMv_heatg2.isel(lon=0,lat=0)[i,:]; kppdepth_heatg2 = GTMdepth_heatg2[:]
	# kppv_heat05dlta = GTMv_heat05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat05dlta = GTMdepth_heat05dlta[:]
	# kppv_heat025dlta = GTMv_heat025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat025dlta = GTMdepth_heat025dlta[:]
	# kppv_heat125dlta = GTMv_heat125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat125dlta = GTMdepth_heat125dlta[:]
	# kppv_heat0625dlta = GTMv_heat0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat0625dlta = GTMdepth_heat0625dlta[:]
	# kppv_mh = GTMv_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
	# kppv_mhLoc = GTMv_mhLoc.isel(lon=0,lat=0)[i,:]; kppdepth_mhLoc = GTMdepth_mhLoc[:]
	# kppv_mhg05 = GTMv_mhg05.isel(lon=0,lat=0)[i,:]; kppdepth_mhg05 = GTMdepth_mhg05[:]
	# kppv_mhg2 = GTMv_mhg2.isel(lon=0,lat=0)[i,:]; kppdepth_mhg2 = GTMdepth_mhg2[:]
	# kppv_mh05dlta = GTMv_mh05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh05dlta = GTMdepth_mh05dlta[:]
	# kppv_mh025dlta = GTMv_mh025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh025dlta = GTMdepth_mh025dlta[:]
	# kppv_mh125dlta = GTMv_mh125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh125dlta = GTMdepth_mh125dlta[:]
	# kppv_mh0625dlta = GTMv_mh0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh0625dlta = GTMdepth_mh0625dlta[:]
	kepsv1d = GTMkeu.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kpptime = time2[i]
	# kpptimeLa = time2La[i]
	# kpptimeLag05 = time2Lag05[i]
	# kpptimeLag2 = time2Lag2[i]
	# kpptimeLa05dlta = time2La05dlta[i]
	# kpptimeLa025dlta = time2La025dlta[i]
	# kpptimeLa125dlta = time2La125dlta[i]
	# kpptimeLa0625dlta = time2La0625dlta[i]
	# kpptime_heat = time2_heat[i]
	# kpptime_heatLoc = time2_heatLoc[i]
	# kpptime_heatg05 = time2_heatg05[i]
	# kpptime_heatg2 = time2_heatg2[i]
	# kpptime_heat05dlta = time2_heat05dlta[i]
	# kpptime_heat025dlta = time2_heat025dlta[i]
	# kpptime_heat125dlta = time2_heat125dlta[i]		
	# kpptime_heat0625dlta = time2_heat0625dlta[i]		
	# kpptime_mh = time2_mh[i]
	# kpptime_mhLoc = time2_mhLoc[i]
	# kpptime_mhg05 = time2_mhg05[i]
	# kpptime_mhg2 = time2_mhg2[i]
	# kpptime_mh05dlta = time2_mh05dlta[i]
	# kpptime_mh025dlta = time2_mh025dlta[i]
	# kpptime_mh125dlta = time2_mh125dlta[i]
	# kpptime_mh0625dlta = time2_mh0625dlta[i]
	kepstime = ketime2[i]
plt.plot(kppv1d,kppdepth,color='gray',linestyle='solid',label = 'KPP d=%0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppvLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppvLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heat,kppdepth_heat,color='firebrick',linestyle='solid',label = 'KPP HD d=%0.2f'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppv_heatLoc,kppdepth_heatLoc,color='olive',linestyle='dashed',label = 'KPP HD loc d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppv_heatLoc1,kppdepth_heatLoc1,color='grey',label = 'KPP HD loc1 d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heatg2,kppdepth_heatg2,color='firebrick',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_heatg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heatg05,kppdepth_heatg05,color='firebrick',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_heatg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heat05dlta,kppdepth_heat05dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_heat05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heat025dlta,kppdepth_heat025dlta,color='firebrick',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_heat025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heat125dlta,kppdepth_heat125dlta,color='firebrick',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_heat125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_heat0625dlta,kppdepth_heat0625dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_heat0625dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mh,kppdepth_mh,color='darkorchid',linestyle='solid',label = 'KPP EV-HD d=%0.2f'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppv_mhLoc,kppdepth_mhLoc,color='darkorchid',linestyle='dashed',label = 'KPP EV-HD loc d=%0.2f'%kpptime_mhLoc.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mhg05,kppdepth_mhg05,color='darkorchid',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_mhg05.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mhg2,kppdepth_mhg2,color='darkorchid',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_mhg2.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mh05dlta,kppdepth_mh05dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_mh05dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mh025dlta,kppdepth_mh025dlta,color='darkorchid',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_mh025dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mh125dlta,kppdepth_mh125dlta,color='darkorchid',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_mh125dlta.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppv_mh0625dlta,kppdepth_mh0625dlta,color='darkorchid',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_mh0625dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsv1d,kepsdepth,color='goldenrod',linestyle='solid' ,label = 'k-eps d=%0.2f'%kepstime.values); plt.legend(loc='best', fontsize='small')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.grid()
plt.xlabel('Spanwise velocity -m/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-150,0])
plt.ylim(-400,0)
plt.savefig('spring_highLatKPP/GOTM_kpp_kppLa_kppheat_kppmh_keps_spanwise_velocity_profile.png')
plt.show()


exit()
# Compared plot GLider and GOTM_KPP
for j in Gldr_sel:
	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
	Gldrtime = new_time[j] 
	for i in GTM_sel:
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
		kpptempLoc = GTMtempLoc.isel(lon=0,lat=0)[i,:]; kppdepthLoc = GTMdepthLoc[:]
		kpptempLa = GTMtempLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
		kpptemp_heat = GTMtemp_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
		kpptemp_mh = GTMtemp_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
		kpptime = time2[i]
		kpptimeLoc = time2Loc[i]
		kpptimeLa = time2La[i]
		kpptime_heat = time2_heat[i]
		kpptime_mh = time2_mh[i]
		kepstime = ketime2[i]
	plt.plot(kpptemp1d,kppdepth,'r-',label = 'GOTM_KPP day %s'%kpptime.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLoc,kppdepthLoc,'k--',label = 'GOTM_KPP local day %s'%kpptime.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptempLa,kppdepthLa,'b--' ,label = 'GOTM_KPP edit Day %s'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat,kppdepth_heat,'b-',label = 'GOTM_KPP diff day %s'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_mh,kppdepth_mh,'k:',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kepstemp1d,kepsdepth,'k:',label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc='best', fontsize='small')
	plt.plot(Gldrtemp1d,Gldrdepth,'g-',label='Glider day %s'%Gldrtime.values); plt.legend(loc='best', fontsize='small')
	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
	# plt.title('Compared Temperature Profiles : Glider and GOTM w KPP,k-eps')
	plt.axis([11,15.5,-150,0])
plt.savefig('spring_highLatKPP/Gldr_GOTM_kpp_kpploc_kppLa_keps_Temp_profile_SpringhighLat.png')
plt.show()

## EDDY VISCOSITY
for i in GTM_sel:
	kppnum1d = GTMnum.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kppnumLa = GTMnumLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	kppnum_heat = GTMnum_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
	kppnum_mh = GTMnum_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
	# kppnum_heatg5 = GTMnum_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
	# kppnum_heatg15 = GTMnum_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
	kepsnum1d = GTMkenum.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kpptime = time2[i]
	kpptimeLa = time2La[i]
	kpptime_heat = time2_heat[i]
	kpptime_mh = time2_mh[i]
	# kpptime_heatg5 = time2_heatg5[i]
	# kpptime_heatg15 = time2_heatg15[i]
	kepstime = ketime2[i]
plt.plot(kppnum1d,kppdepth,'k-',label = 'GOTM_KPP day %s'%kpptime.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLa,kppdepthLa,'b-' ,label = 'GOTM_KPP visc day %s'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum_heat,kppdepth_heat,'g--',label = 'GOTM_KPP diff day %s'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum_mh,kppdepth_mh,'r--',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heatg5,kppdepth_heatg5,'r-',label = 'GOTM_KPP diff gam=0.5 day %s'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heatg15,kppdepth_heatg15,'r--',label = 'GOTM_KPP diff gam=1.5 day %s'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnum1d,kepsdepth,'k:' ,label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc='best', fontsize='small')
plt.xlabel('Eddy Viscosity -m2/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-150,0])
plt.ylim(-150,0)
plt.savefig('spring_highLatKPP/GOTM_kpp_kpploc_kppvisc_kppheat_kppmh_keps_eddy_viscosity_profile.png')
plt.show()

## HEAT DIFFUSIVITY
for i in GTM_sel:
	kppnuh1d = GTMnuh.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kppnuhLa = GTMnuhLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	kppnuh_heat = GTMnuh_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
	kppnuh_mh = GTMnuh_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
	# kppnuh_heatg5 = GTMnuh_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
	# kppnuh_heatg15 = GTMnuh_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
	kepsnuh1d = GTMkenuh.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kpptime = time2[i]
	kpptimeLa = time2La[i]
	kpptime_heat = time2_heat[i]
	kpptime_mh = time2_mh[i]
	# kpptime_heatg5 = time2_heatg5[i]
	# kpptime_heatg15 = time2_heatg15[i]
	kepstime = ketime2[i]
plt.plot(kppnuh1d,kppdepth,'k-',label = 'GOTM_KPP day %s'%kpptime.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa,kppdepthLa,'b-' ,label = 'GOTM_KPP visc day %s'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh_heat,kppdepth_heat,'g--',label = 'GOTM_KPP diff day %s'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh_mh,kppdepth_mh,'r--',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg5,kppdepth_heatg5,'r-',label = 'GOTM_KPP diff gam=0.5 day %s'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg15,kppdepth_heatg15,'r--',label = 'GOTM_KPP diff gam=1.5 day %s'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnuh1d,kepsdepth,'k:' ,label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc='best', fontsize='small')
plt.xlabel('Heat Diffusivity -m2/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-150,0])
plt.ylim(-150,0)
plt.savefig('spring_highLatKPP/GOTM_kpp_kpploc_kppLa_kppheat_kppmh_keps_heat_diffusivity_profile.png')
plt.show()
exit()


exit()