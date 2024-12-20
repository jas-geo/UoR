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

plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

### ----------------------------- GLIDER ----------------------------- ##
# file1 = 'C:/home/users/mc837749/Documents/gotm-4.0.0/simulations/annualOSMOSISforcing_2502/glider_timeseries.nc'
# gliderdata= Dataset('../glider_timeseries.nc'); 
# gliderdata= xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc'); 
glider = xr.open_dataset('~/Documents/gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; GldrDens = glider.pot_den; 
Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen
# print('original Gldr temp', GldrTemp.shape)
# GldrTempNaN = GldrTemp.isnull(); print('number of missing values in GldrTemp', GldrTempNaN)

# GldrTemp = GldrTemp.to_dataframe().dropna(how='all'); print('Gldr temp with dropped NaN values?', GldrTemp, GldrTemp.shape)
# GldrTemp = GldrTemp.to_dataframe();

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year 

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



# # # check variables in file
# print(gliderdata.variables.items())
# print(GldrDens.values)
# exit()
### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
## LaT addition ##
gotmkppLaRib = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat.nc",decode_times=False)
# gotmkppnum2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_num15.nc",decode_times=False)
# gotmkppnuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nuh15.nc",decode_times=False)
# # gotmkppRibnum2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_LaT_Rib_num2.nc",decode_times=False)
# gotmkppRib05num2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_Rib_num15.nc",decode_times=False)
# # gotmkppRibnuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_LaT_Rib_nuh2.nc",decode_times=False)
# gotmkppRib2nuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_Rib_nuh15.nc",decode_times=False)
## h/L correction ##
gotmkppLaRibhLPveP100 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP100.nc",decode_times=False)
# gotmkppLaRibhLPveP200 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP200.nc",decode_times=False)
# gotmkppLaRibhLPveP1000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1000.nc",decode_times=False)
# gotmkppLaRibhLPveP1300 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1300.nc",decode_times=False)
# gotmkppLaRibhLPveP1350 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1350.nc",decode_times=False)
# gotmkppLaRibhLPveP1400 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1400.nc",decode_times=False)
gotmkppLaRibhLPveP1500 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP1500.nc",decode_times=False)
# gotmkppLaRibhLPveP1600 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1600.nc",decode_times=False)
# gotmkppLaRibhLPveP1750 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1750.nc",decode_times=False)
# gotmkppLaRibhLPveP2000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP2000.nc",decode_times=False)
# gotmkppLaRibhLPveP4000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP4000.nc",decode_times=False)
# gotmkppLaRibhLPveP5000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP5000.nc",decode_times=False)
gotmkppLaRibhLPveP10000 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP10000.nc",decode_times=False)
# gotmkppLaRibhLPveP10000shear = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_hLPveP10000_salinity.nc",decode_times=False)

# # check variables in files
# print(gotmkpp.variables.items())

#######################
## collect variables ##
#######################

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMsal = gotmkpp.salt; GTMo2 = gotmkpp.o2_obs; GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMxflx = gotmkpp.gamu; GTMyflx = gotmkpp.gamv; time = gotmkpp.time
# GTMtempL = gotmkppL.temp; GTMdepthL = gotmkppL.z; GTMnumL = gotmkppL.num; GTMnuhL = gotmkppL.nuh; GTMsalL = gotmkppL.salt; GTMo2L = gotmkppL.o2_obs; GTMuL = gotmkppL.u; GTMvL = gotmkppL.v; GTMxflxL = gotmkppL.gamu; GTMyflxL = gotmkppL.gamv; timeL = gotmkppL.time
# GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time
# GTMtempLa = gotmkppLa.temp; GTMdepthLa = gotmkppLa.z; GTMnumLa = gotmkppLa.num; GTMnuhLa = gotmkppLa.nuh; GTMuLa = gotmkppLa.u; GTMvLa = gotmkppLa.v; timeLa = gotmkppLa.time
GTMtempLaRib = gotmkppLaRib.temp; GTMdepthLaRib = gotmkppLaRib.z; GTMnumLaRib = gotmkppLaRib.num; GTMnuhLaRib = gotmkppLaRib.nuh; GTMuLaRib = gotmkppLaRib.u; GTMvLaRib = gotmkppLaRib.v; timeLaRib = gotmkppLaRib.time
# GTMtempLaRibL = gotmkppLaRibL.temp; GTMdepthLaRibL = gotmkppLaRibL.z; GTMnumLaRibL = gotmkppLaRibL.num; GTMnuhLaRibL = gotmkppLaRibL.nuh; GTMuLaRibL = gotmkppLaRibL.u; GTMvLaRibL = gotmkppLaRibL.v; timeLaRibL = gotmkppLaRibL.time
# GTMtempnum2 = gotmkppnum2.temp; GTMdepthnum2 = gotmkppnum2.z; GTMnumnum2 = gotmkppnum2.num; GTMnuhnum2 = gotmkppnum2.nuh; GTMunum2 = gotmkppnum2.u; GTMvnum2 = gotmkppnum2.v; timenum2 = gotmkppnum2.time
# GTMtempnuh2 = gotmkppnuh2.temp; GTMdepthnuh2 = gotmkppnuh2.z; GTMnumnuh2 = gotmkppnuh2.num; GTMnuhnuh2 = gotmkppnuh2.nuh; GTMunuh2 = gotmkppnuh2.u; GTMvnuh2 = gotmkppnuh2.v; timenuh2 = gotmkppnuh2.time
# # GTMtempRibnum2 = gotmkppRibnum2.temp; GTMdepthRibnum2 = gotmkppRibnum2.z; GTMnumRibnum2 = gotmkppRibnum2.num; GTMnuhRibnum2 = gotmkppRibnum2.nuh; GTMuRibnum2 = gotmkppRibnum2.u; GTMvRibnum2 = gotmkppRibnum2.v; timeRibnum2 = gotmkppRibnum2.time
# GTMtempRib05num2 = gotmkppRib05num2.temp; GTMdepthRib05num2 = gotmkppRib05num2.z; GTMnumRib05num2 = gotmkppRib05num2.num; GTMnuhRib05num2 = gotmkppRib05num2.nuh; GTMuRib05num2 = gotmkppRib05num2.u; GTMvRib05num2 = gotmkppRib05num2.v; timeRib05num2 = gotmkppRib05num2.time
# # GTMtempRibnuh2 = gotmkppRibnuh2.temp; GTMdepthRibnuh2 = gotmkppRibnuh2.z; GTMnumRibnuh2 = gotmkppRibnuh2.num; GTMnuhRibnuh2 = gotmkppRibnuh2.nuh; GTMuRibnuh2 = gotmkppRibnuh2.u; GTMvRibnuh2 = gotmkppRibnuh2.v; timeRibnuh2 = gotmkppRibnuh2.time
# GTMtempRib2nuh2 = gotmkppRib2nuh2.temp; GTMdepthRib2nuh2 = gotmkppRib2nuh2.z; GTMnumRib2nuh2 = gotmkppRib2nuh2.num; GTMnuhRib2nuh2 = gotmkppRib2nuh2.nuh; GTMuRib2nuh2 = gotmkppRib2nuh2.u; GTMvRib2nuh2 = gotmkppRib2nuh2.v; timeRib2nuh2 = gotmkppRib2nuh2.time

# GTMtempLaRibhLPveN100 = gotmkppLaRibhLPveN100.temp; GTMdepthLaRibhLPveN100 = gotmkppLaRibhLPveN100.z; GTMnumLaRibhLPveN100 = gotmkppLaRibhLPveN100.num; GTMnuhLaRibhLPveN100 = gotmkppLaRibhLPveN100.nuh; GTMuLaRibhLPveN100 = gotmkppLaRibhLPveN100.u; GTMvLaRibhLPveN100 = gotmkppLaRibhLPveN100.v; timeLaRibhLPveN100 = gotmkppLaRibhLPveN100.time
GTMtempLaRibhLPveP100 = gotmkppLaRibhLPveP100.temp; GTMdepthLaRibhLPveP100 = gotmkppLaRibhLPveP100.z; GTMnumLaRibhLPveP100 = gotmkppLaRibhLPveP100.num; GTMnuhLaRibhLPveP100 = gotmkppLaRibhLPveP100.nuh; GTMuLaRibhLPveP100 = gotmkppLaRibhLPveP100.u; GTMvLaRibhLPveP100 = gotmkppLaRibhLPveP100.v; timeLaRibhLPveP100 = gotmkppLaRibhLPveP100.time
# GTMtempLaRibhLPveP200 = gotmkppLaRibhLPveP200.temp; GTMdepthLaRibhLPveP200 = gotmkppLaRibhLPveP200.z; GTMnumLaRibhLPveP200 = gotmkppLaRibhLPveP200.num; GTMnuhLaRibhLPveP200 = gotmkppLaRibhLPveP200.nuh; GTMuLaRibhLPveP200 = gotmkppLaRibhLPveP200.u; GTMvLaRibhLPveP200 = gotmkppLaRibhLPveP200.v; timeLaRibhLPveP200 = gotmkppLaRibhLPveP200.time
# GTMtempLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.temp; GTMdepthLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.z; GTMnumLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.num; GTMnuhLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.nuh; GTMuLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.u; GTMvLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.v; timeLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.time
# GTMtempLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.temp; GTMdepthLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.z; GTMnumLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.num; GTMnuhLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.nuh; GTMuLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.u; GTMvLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.v; timeLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.time
# GTMtempLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.temp; GTMdepthLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.z; GTMnumLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.num; GTMnuhLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.nuh; GTMuLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.u; GTMvLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.v; timeLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.time
# GTMtempLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.temp; GTMdepthLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.z; GTMnumLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.num; GTMnuhLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.nuh; GTMuLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.u; GTMvLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.v; timeLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.time
GTMtempLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.temp; GTMdepthLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.z; GTMnumLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.num; GTMnuhLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.nuh; GTMuLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u; GTMvLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.v; timeLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.time
# GTMtempLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.temp; GTMdepthLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.z; GTMnumLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.num; GTMnuhLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.nuh; GTMuLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.u; GTMvLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.v; timeLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.time
# GTMtempLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.temp; GTMdepthLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.z; GTMnumLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.num; GTMnuhLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.nuh; GTMuLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.u; GTMvLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.v; timeLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.time
# GTMtempLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.temp; GTMdepthLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.z; GTMnumLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.num; GTMnuhLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.nuh; GTMuLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.u; GTMvLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.v; timeLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.time
# GTMtempLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.temp; GTMdepthLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.z; GTMnumLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.num; GTMnuhLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.nuh; GTMuLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.u; GTMvLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.v; timeLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.time
# GTMtempLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.temp; GTMdepthLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.z; GTMnumLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.num; GTMnuhLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.nuh; GTMuLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.u; GTMvLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.v; timeLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.time
GTMtempLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.temp; GTMdepthLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.z; GTMnumLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.num; GTMnuhLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.nuh; GTMuLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u; GTMvLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.v; timeLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.time
# GTMtempLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.temp; GTMdepthLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.z; GTMnumLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.num; GTMnuhLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.nuh; GTMuLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.u; GTMvLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.v; timeLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.time


# print('GTM KPP time :', time.values)
# exit()
##################
## convert time ##
##################

time2 = 267.75 + (time[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
# time2L = 267.75 + (timeL[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
# time2La = 445.0 + (timeLa[:]/86400); 
time2LaRib = 267.75 + (timeLaRib[:]/86400); 
# time2LaRibL = 267.75 + (timeLaRibL[:]/86400); 
# time2num2 = 267.75 + (timenum2[:]/86400); 
# time2nuh2 = 267.75 + (timenuh2[:]/86400); 
# # time2Ribnum2 = 267.75 + (timeRibnum2[:]/86400); 
# time2Rib05num2 = 267.75 + (timeRib05num2[:]/86400); 
# # time2Ribnuh2 = 267.75 + (timeRibnuh2[:]/86400); 
# time2Rib2nuh2 = 267.75 + (timeRib2nuh2[:]/86400); 

time2LaRibhLPveP100 = 267.75 + (timeLaRibhLPveP100[:]/86400); 
# time2LaRibhLPveP200 = 267.75 + (timeLaRibhLPveP200[:]/86400); 
# time2LaRibhLPveP1000 = 267.75 + (timeLaRibhLPveP1000[:]/86400); 
# time2LaRibhLPveP1300 = 267.75 + (timeLaRibhLPveP1300[:]/86400); 
# time2LaRibhLPveP1350 = 267.75 + (timeLaRibhLPveP1350[:]/86400); 
# time2LaRibhLPveP1400 = 267.75 + (timeLaRibhLPveP1400[:]/86400); 
time2LaRibhLPveP1500 = 267.75 + (timeLaRibhLPveP1500[:]/86400); 
# time2LaRibhLPveP1600 = 267.75 + (timeLaRibhLPveP1600[:]/86400); 
# time2LaRibhLPveP1750 = 267.75 + (timeLaRibhLPveP1750[:]/86400); 
# time2LaRibhLPveP2000 = 267.75 + (timeLaRibhLPveP2000[:]/86400); 
# time2LaRibhLPveP4000 = 267.75 + (timeLaRibhLPveP4000[:]/86400); 
# time2LaRibhLPveP5000 = 267.75 + (timeLaRibhLPveP5000[:]/86400); 
time2LaRibhLPveP10000 = 267.75 + (timeLaRibhLPveP10000[:]/86400); 
# time2LaRibhLPveP10000shear = 267.75 + (timeLaRibhLPveP10000shear[:]/86400); 



datetime2 = []

def julian_to_GTMdate():
	initial_date = "2012-01-01 00:00:00"
	global datetime2
	for i in time2LaRib.values:
		date = pd.to_datetime(initial_date) + pd.DateOffset(days=i)
		datetime2.append(date)
		# print(date)
	return print('Time format is in Date-Time == SUCCESS')

julian_to_GTMdate()






###################
# time allocation #
###################

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
# gotmkppL = gotmkppL.assign(time2L=("time", time2L)); gotmkppL = gotmkppL.swap_dims({"time" : "time2L"})
# gotmkppLa = gotmkppLa.assign(time2La=("time", time2La)); gotmkppLa = gotmkppLa.swap_dims({"time" : "time2La"})
gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})
# gotmkppLaRibL = gotmkppLaRibL.assign(time2LaRibL=("time", time2LaRibL)); gotmkppLaRibL = gotmkppLaRibL.swap_dims({"time" : "time2LaRibL"})
# gotmkppnum2 = gotmkppnum2.assign(time2num2=("time", time2num2)); gotmkppnum2 = gotmkppnum2.swap_dims({"time" : "time2num2"})
# gotmkppnuh2 = gotmkppnuh2.assign(time2nuh2=("time", time2nuh2)); gotmkppnuh2 = gotmkppnuh2.swap_dims({"time" : "time2nuh2"})
# # gotmkppRibnum2 = gotmkppRibnum2.assign(time2Ribnum2=("time", time2Ribnum2)); gotmkppRibnum2 = gotmkppRibnum2.swap_dims({"time" : "time2Ribnum2"})
# gotmkppRib05num2 = gotmkppRib05num2.assign(time2Rib05num2=("time", time2Rib05num2)); gotmkppRib05num2 = gotmkppRib05num2.swap_dims({"time" : "time2Rib05num2"})
# # gotmkppRibnuh2 = gotmkppRibnuh2.assign(time2Ribnuh2=("time", time2Ribnuh2)); gotmkppRibnuh2 = gotmkppRibnuh2.swap_dims({"time" : "time2Ribnuh2"})
# gotmkppRib2nuh2 = gotmkppRib2nuh2.assign(time2Rib2nuh2=("time", time2Rib2nuh2)); gotmkppRib2nuh2 = gotmkppRib2nuh2.swap_dims({"time" : "time2Rib2nuh2"})

gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.assign(time2LaRibhLPveP100=("time", time2LaRibhLPveP100)); gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.swap_dims({"time" : "time2LaRibhLPveP100"})
# gotmkppLaRibhLPveP200 = gotmkppLaRibhLPveP200.assign(time2LaRibhLPveP200=("time", time2LaRibhLPveP200)); gotmkppLaRibhLPveP200 = gotmkppLaRibhLPveP200.swap_dims({"time" : "time2LaRibhLPveP200"})
# gotmkppLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.assign(time2LaRibhLPveP1000=("time", time2LaRibhLPveP1000)); gotmkppLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.swap_dims({"time" : "time2LaRibhLPveP1000"})
# gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.assign(time2LaRibhLPveP2000=("time", time2LaRibhLPveP2000)); gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.swap_dims({"time" : "time2LaRibhLPveP2000"})
# gotmkppLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.assign(time2LaRibhLPveP1400=("time", time2LaRibhLPveP1400)); gotmkppLaRibhLPveP1400 = gotmkppLaRibhLPveP1400.swap_dims({"time" : "time2LaRibhLPveP1400"})
gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.assign(time2LaRibhLPveP1500=("time", time2LaRibhLPveP1500)); gotmkppLaRibhLPveP1500 = gotmkppLaRibhLPveP1500.swap_dims({"time" : "time2LaRibhLPveP1500"})
# gotmkppLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.assign(time2LaRibhLPveP1600=("time", time2LaRibhLPveP1600)); gotmkppLaRibhLPveP1600 = gotmkppLaRibhLPveP1600.swap_dims({"time" : "time2LaRibhLPveP1600"})
# gotmkppLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.assign(time2LaRibhLPveP1300=("time", time2LaRibhLPveP1300)); gotmkppLaRibhLPveP1300 = gotmkppLaRibhLPveP1300.swap_dims({"time" : "time2LaRibhLPveP1300"})
# gotmkppLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.assign(time2LaRibhLPveP1750=("time", time2LaRibhLPveP1750)); gotmkppLaRibhLPveP1750 = gotmkppLaRibhLPveP1750.swap_dims({"time" : "time2LaRibhLPveP1750"})
# gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.assign(time2LaRibhLPveP2000=("time", time2LaRibhLPveP2000)); gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.swap_dims({"time" : "time2LaRibhLPveP2000"})
# gotmkppLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.assign(time2LaRibhLPveP4000=("time", time2LaRibhLPveP4000)); gotmkppLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.swap_dims({"time" : "time2LaRibhLPveP4000"})
# gotmkppLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.assign(time2LaRibhLPveP1350=("time", time2LaRibhLPveP1350)); gotmkppLaRibhLPveP1350 = gotmkppLaRibhLPveP1350.swap_dims({"time" : "time2LaRibhLPveP1350"})
# gotmkppLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.assign(time2LaRibhLPveP5000=("time", time2LaRibhLPveP5000)); gotmkppLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.swap_dims({"time" : "time2LaRibhLPveP5000"})
gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.assign(time2LaRibhLPveP10000=("time", time2LaRibhLPveP10000)); gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.swap_dims({"time" : "time2LaRibhLPveP10000"})
# gotmkppLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.assign(time2LaRibhLPveP10000shear=("time", time2LaRibhLPveP10000shear)); gotmkppLaRibhLPveP10000shear = gotmkppLaRibhLPveP10000shear.swap_dims({"time" : "time2LaRibhLPveP10000shear"})

# print(gotmkppLaRibhLPveP2000.variables.items())

#### --------------- testing CT difference ---------------- ####
# print(glider.variables.items())
# print('shape of kpp temp', GTMtemp.shape); print('kpp depth',GTMdepth[0].values); 
# print('shape of kpp time', time2[:143].values)
# print('shape of glider temp', GldrTemp[0:5,0:50].values, GldrTemp.shape);
# print('shape of kpp temp', GTMtemp[0:5,400:500].values, GTMtemp[0:5,400:500].shape);
# print('shape of glider depth', GldrDepth[0:200], GldrDepth[0:200].shape); 
# print('shape of kpp depth', GTMdepth[400:500].values, GTMdepth[400:500].shape); 
# # print('shape of glider time', new_time[6:100], new_time[6:100].shape)
# print('shape of kpp time', time2[143:1007].values, time2[143:1011].shape)
# # Gldrtemp_new = GldrTemp[0:4,0:5]; Gldrz_new = GldrDepth[0:5]; Gldrtime_new = new_time[0:4]
# # print(Gldrtemp_new); print(Gldrz_new); print(Gldrtime_new)
# print('GTM time', time2[100:300])
# print('Gldr temperature', GldrTemp.values)
# print('GOTM temp shape', GTMtemp)
# print('GOTM depth shape', GTMdepth.values)

# ###############################################################################################
# ### ----------------------------- calculate isothermal layer depth ----------------------------- ###
# ###############################################################################################

# ## using KPP model
# ild_densGldr = GldrDens.sel(pressure=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GldrDens-ild_densGldr).argmin('z'); # print(z_indexes)
# # print('GOTM temp MLD indexes', z_indexes)
# ild_Gldr = GldrDens.z[z_indexes]; 
# # ild_Gldr = ild_Gldr.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


###############################################################################################
### ----------------------------- calculate mixed layer depth ----------------------------- ###
###############################################################################################

## Psedocode 
# 1. extract temperature values for depth z_ref
# 2. new_temp = temp(z_ref) - delta_t
# 3. depth at new_temp

# MLD_DT02 = depth where (T = T(z_ref) ± delta_t)
# z_ref = -10 -- reference depth;   T(z_ref) -- temperature at reference depth          
# delta_t = 0.2 -- ;    

# df = GldrTemp.to_series()
# print(df)
# exit()
# interpGldrTemp = df.interpolate(method='linear', limit_direction='forward')
# print('interpolate NaN GldrTemp', interpGldrTemp.shape)
# exit()

# ## using OSMOSIS
# mld_tempOS = GldrTemp.sel(pressure=11) - 0.2;
# # print('OSMOSIS MLD',mld_tempOS)
# z_indexes = abs(interpGldrTemp-mld_tempOS).argmin('pressure');
# exit()
# mldOS = sliceGldrTemp.pressure[z_indexes];
# mldOS = mldOS.isel();
# # print('OSMOSIS MLD', mldOS)

## using KPP model
mld_tempKPP = GTMtemp.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMtemp-mld_tempKPP).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
mldKPP = GTMtemp.z[z_indexes]; 
mldKPP = mldKPP.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model LOCAL
# mld_tempKPP_Local = GTMtempL.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempL-mld_tempKPP_Local).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP_Local = GTMtempL.z[z_indexes]; 
# mldKPP_Local = mldKPP_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# plt.plot(time2,mldKPP); plt.xlabel('time -julian day'); plt.ylabel('mixed layer depth -m')
# plt.axis([260,341,-140,0]); plt.title('GOTM Mixed Layer Depth')
# plt.savefig('GOTM_Mixed_Layer_Depth_test6.png'); plt.show()

# ## using KPP EV model
# mld_tempKPPLa = GTMtempLa.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLa-mld_tempKPPLa).argmin('z'); # print(z_indexes)
# mldKPPLa = GTMtempLa.z[z_indexes]; 
# mldKPPLa = mldKPPLa.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original
mld_tempKPPLaRib = GTMtempLaRib.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMtempLaRib-mld_tempKPPLaRib).argmin('z'); # print(z_indexes)
mldKPPLaRib = GTMtempLaRib.z[z_indexes]; 
mldKPPLaRib = mldKPPLaRib.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original
# mld_tempKPPLaRib_Local = GTMtempLaRibL.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibL-mld_tempKPPLaRib_Local).argmin('z'); # print(z_indexes)
# mldKPPLaRib_Local = GTMtempLaRib.z[z_indexes]; 
# mldKPPLaRib_Local = mldKPPLaRib_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use num*2
# mld_tempKPPnum2 = GTMtempnum2.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempnum2-mld_tempKPPnum2).argmin('z'); # print(z_indexes)
# mldKPPnum2 = GTMtempnum2.z[z_indexes]; 
# mldKPPnum2 = mldKPPnum2.isel(lat=0,lon=0); # print('mixed yer depth :',mld, mld.shape)

# ## using KPP model - use nuh*2
# mld_tempKPPnuh2 = GTMtempnuh2.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempnuh2-mld_tempKPPnuh2).argmin('z'); # print(z_indexes)
# mldKPPnuh2 = GTMtempnuh2.z[z_indexes]; 
# mldKPPnuh2 = mldKPPnuh2.isel(lat=0,lon=0); # print('mixed yer depth :',mld, mld.shape)

# # ## using KPP model - use t and Rib num*2
# # mld_tempKPPRibnum2 = GTMtempRibnum2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtempRibnum2-mld_tempKPPRibnum2).argmin('z'); # print(z_indexes)
# # mldKPPRibnum2 = GTMtempRibnum2.z[z_indexes]; 
# # mldKPPRibnum2 = mldKPPRibnum2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use t and Rib num*2
# mld_tempKPPRib05num2 = GTMtempRib05num2.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempRib05num2-mld_tempKPPRib05num2).argmin('z'); # print(z_indexes)
# mldKPPRib05num2 = GTMtempRib05num2.z[z_indexes]; 
# mldKPPRib05num2 = mldKPPRib05num2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# # ## using KPP model - use t and Rib nuh*2
# # mld_tempKPPRibnuh2 = GTMtempRibnuh2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtempRibnuh2-mld_tempKPPRibnuh2).argmin('z'); # print(z_indexes)
# # mldKPPRibnuh2 = GTMtempRibnuh2.z[z_indexes]; 
# # mldKPPRibnuh2 = mldKPPRibnuh2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use t and Rib nuh*2
# mld_tempKPPRib2nuh2 = GTMtempRib2nuh2.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempRib2nuh2-mld_tempKPPRib2nuh2).argmin('z'); # print(z_indexes)
# mldKPPRib2nuh2 = GTMtempRib2nuh2.z[z_indexes]; 
# mldKPPRib2nuh2 = mldKPPRib2nuh2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100)
mld_tempKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMtempLaRibhLPveP100-mld_tempKPPLaRibhLPveP100).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.z[z_indexes]; 
mldKPPLaRibhLPveP100 = mldKPPLaRibhLPveP100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>200) for gamma
# mld_tempKPPLaRibhLPveP200 = GTMtempLaRibhLPveP200.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP200-mld_tempKPPLaRibhLPveP200).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP200 = GTMtempLaRibhLPveP200.z[z_indexes]; 
# mldKPPLaRibhLPveP200 = mldKPPLaRibhLPveP200.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100)
# mld_tempKPPhLPveP100 = GTMtemphLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemphLPveP100-mld_tempKPPhLPveP100).argmin('z'); # print(z_indexes)
# mldKPPhLPveP100 = GTMtemphLPveP100.z[z_indexes]; 
# mldKPPhLPveP100 = mldKPPhLPveP100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>1000)
# mld_tempKPPLaRibhLPveP1000 = GTMtempLaRibhLPveP1000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP1000-mld_tempKPPLaRibhLPveP1000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1000 = GTMtempLaRibhLPveP1000.z[z_indexes]; 
# mldKPPLaRibhLPveP1000 = mldKPPLaRibhLPveP1000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>2000) for gamma
# mld_tempKPPLaRibhLPveP2000 = GTMtempLaRibhLPveP2000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP2000-mld_tempKPPLaRibhLPveP2000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP2000 = GTMtempLaRibhLPveP2000.z[z_indexes]; 
# mldKPPLaRibhLPveP2000 = mldKPPLaRibhLPveP2000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP1400 = GTMtempLaRibhLPveP1400.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP1400-mld_tempKPPLaRibhLPveP1400).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1400 = GTMtempLaRibhLPveP1400.z[z_indexes]; 
# mldKPPLaRibhLPveP1400 = mldKPPLaRibhLPveP1400.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>1500) for gamma
mld_tempKPPLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMtempLaRibhLPveP1500-mld_tempKPPLaRibhLPveP1500).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP1500 = GTMtempLaRibhLPveP1500.z[z_indexes]; 
mldKPPLaRibhLPveP1500 = mldKPPLaRibhLPveP1500.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>1300) for gamma
# mld_tempKPPLaRibhLPveP1600 = GTMtempLaRibhLPveP1600.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP1600-mld_tempKPPLaRibhLPveP1600).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1600 = GTMtempLaRibhLPveP1600.z[z_indexes]; 
# mldKPPLaRibhLPveP1600 = mldKPPLaRibhLPveP1600.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP1300 = GTMtempLaRibhLPveP1300.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP1300-mld_tempKPPLaRibhLPveP1300).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1300 = GTMtempLaRibhLPveP1300.z[z_indexes]; 
# mldKPPLaRibhLPveP1300 = mldKPPLaRibhLPveP1300.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>500)
# mld_tempKPPLaRibhLPveP1750 = GTMtempLaRibhLPveP1750.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP1750-mld_tempKPPLaRibhLPveP1750).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1750 = GTMtempLaRibhLPveP1750.z[z_indexes]; 
# mldKPPLaRibhLPveP1750 = mldKPPLaRibhLPveP1750.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>2000) for gamma
# mld_tempKPPLaRibhLPveP2000 = GTMtempLaRibhLPveP2000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP2000-mld_tempKPPLaRibhLPveP2000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP2000 = GTMtempLaRibhLPveP2000.z[z_indexes]; 
# mldKPPLaRibhLPveP2000 = mldKPPLaRibhLPveP2000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP4000 = GTMtempLaRibhLPveP4000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP4000-mld_tempKPPLaRibhLPveP4000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP4000 = GTMtempLaRibhLPveP4000.z[z_indexes]; 
# mldKPPLaRibhLPveP4000 = mldKPPLaRibhLPveP4000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP1350 = GTMtempLaRibhLPveP1350.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP1350-mld_tempKPPLaRibhLPveP1350).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1350 = GTMtempLaRibhLPveP1350.z[z_indexes]; 
# mldKPPLaRibhLPveP1350 = mldKPPLaRibhLPveP1350.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>500)
# mld_tempKPPLaRibhLPveP5000 = GTMtempLaRibhLPveP5000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP5000-mld_tempKPPLaRibhLPveP5000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP5000 = GTMtempLaRibhLPveP5000.z[z_indexes]; 
# mldKPPLaRibhLPveP5000 = mldKPPLaRibhLPveP5000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.sel(z=-11.0, method='nearest') - 0.2; 
z_indexes = abs(GTMtempLaRibhLPveP10000-mld_tempKPPLaRibhLPveP10000).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.z[z_indexes]; 
mldKPPLaRibhLPveP10000 = mldKPPLaRibhLPveP10000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP10000shear = GTMtempLaRibhLPveP10000shear.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtempLaRibhLPveP10000shear-mld_tempKPPLaRibhLPveP10000shear).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP10000shear = GTMtempLaRibhLPveP10000shear.z[z_indexes]; 
# mldKPPLaRibhLPveP10000shear = mldKPPLaRibhLPveP10000shear.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## testing for collapse of MLD

MLDLow = np.where((mldKPPLaRibhLPveP10000 >= -1))
# print(r'collapse of MLD :', MLDLow)
# print(gotmkpp.variables.items())


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
	# datess= np.append(datess, datetime.strptime(DT_m[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_date)
# yearmonth = str(yearM) + "-" + str(monthM); 
#print('year and month', yearmonth)
hoursM = np.empty((1465,0), int)
for ti in range(len(mtime)):
	timeDaysM = datetime.strptime(mtime[ti], '%H:%M:%S')#.split(':')
	hoursM = np.append(hoursM, timeDaysM.hour)

hourstodaysM = hoursM/24
# print(hourstodaysM)
convertedDaysM = yearDaysM + hourstodaysM
# print('Converted Days (momentum): ', convertedDaysM[0])

for data in range(len(convertedDaysM)):
	if convertedDaysM[data] < 268.75:
		convertedDaysM[data] =convertedDaysM[data] + 366.25
	else: 
		convertedDaysM[data] == convertedDaysM[data]
convertedDaysM[1464]=633.0

momFlux = np.c_[convertedDaysM,mflux]
# print('converted days (momentum)', convertedDaysM[712:1087])
convertedDaysMlist=convertedDaysM.astype(np.float)

# ### ----------------------------- calculate u_s and u_* ----------------------------- ###
u_s = np.empty((1465,0), float)
fric_u = np.zeros(1465)
## Density of air == 1.225 kg/m**3 
## Density of water == 1000 kg/m**3
# for time in range(len(momFlux)):
# 	u_s = np.append(u_s, math.sqrt((momFlux[time][5])**2 + (momFlux[time][6])**2))
# 	fric_u = np.append(fric_u, math.sqrt(math.sqrt((momFlux[time][3])**2 + (momFlux[time][4])**2)/1000))
# # print('Stokes drift veloctiy :',u_s); print(u_s.shape)
# # print('friction velocity :',fric_u); print(fric_u.shape)

# ### ----------------------------- calculate la_t ----------------------------- ###
la_t = np.zeros(1465)
# for i in range(len(u_s)):
	# la_t = np.append(la_t, math.sqrt(fric_u[i]/u_s[i]))
# print('turbulent Langmuir number :',la_t, la_t.shape)
# print('converted Days :',convertedDays.shape)


### ----------------------------- moving averages h/L ----------------------------- ###

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

# GTM_hL = gotmkpp.hL; GTM_hLR = gotmkpp.hLR; GTM_Lat = gotmkpp.La_t; GTM_tdotus = gotmkpp.tdotus; 
GTM_hL_LaRib = gotmkppLaRib.hL; GTM_hLR_LaRib = gotmkppLaRib.hLR; GTM_Lat_LaRib = gotmkppLaRib.La_t; GTM_tdotus_LaRib = gotmkppLaRib.tdotus; 
GTM_hL_LaRibhLPveP100 = gotmkppLaRibhLPveP100.hL; GTM_hLR_LaRibhLPveP100 = gotmkppLaRibhLPveP100.hLR; GTM_Lat_LaRibhLPveP100 = gotmkppLaRibhLPveP100.La_t; GTM_tdotus_LaRibhLPveP100 = gotmkppLaRibhLPveP100.tdotus; 
# GTM_hL_LaRibhLPveN100 = gotmkppLaRibhLPveN100.hL; GTM_hLR_LaRibhLPveN100 = gotmkppLaRibhLPveN100.hLR; GTM_Lat_LaRibhLPveN100 = gotmkppLaRibhLPveN100.La_t; GTM_tdotus_LaRibhLPveN100 = gotmkppLaRibhLPveN100.tdotus; 
# GTM_hL_LaRibhLPveN200 = gotmkppLaRibhLPveN200.hL; GTM_hLR_LaRibhLPveN200 = gotmkppLaRibhLPveN200.hLR; GTM_Lat_LaRibhLPveN200 = gotmkppLaRibhLPveN200.La_t; GTM_tdotus_LaRibhLPveN200 = gotmkppLaRibhLPveN200.tdotus; 
# GTM_hL_LaRibhLPveN500 = gotmkppLaRibhLPveN500.hL; GTM_hLR_LaRibhLPveN500 = gotmkppLaRibhLPveN500.hLR; GTM_Lat_LaRibhLPveN500 = gotmkppLaRibhLPveN500.La_t; GTM_tdotus_LaRibhLPveN500 = gotmkppLaRibhLPveN500.tdotus; 
# GTM_hL_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.hL; GTM_hLR_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.hLR; GTM_Lat_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.La_t; GTM_tdotus_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.tdotus; 
# GTM_hL_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.hL; GTM_hLR_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.hLR; GTM_Lat_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.La_t; GTM_tdotus_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.tdotus; GTM_Theatflux_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.total; GTM_fric_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.u_taus; GTM_us_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.u_s; 
# GTM_hL_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.hL; GTM_hLR_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.hLR; GTM_Lat_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.La_t; GTM_tdotus_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.tdotus; GTM_Theatflux_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.total; GTM_fric_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.u_taus; GTM_us_LaRibhLPveP1300 = gotmkppLaRibhLPveP1300.u_s; 
# GTM_hL_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.hL; GTM_hLR_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.hLR; GTM_Lat_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.La_t; GTM_tdotus_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.tdotus; GTM_Theatflux_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.total; GTM_fric_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.u_taus; GTM_us_LaRibhLPveP1350 = gotmkppLaRibhLPveP1350.u_s; 
# GTM_hL_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.hL; GTM_hLR_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.hLR; GTM_Lat_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.La_t; GTM_tdotus_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.tdotus; GTM_Theatflux_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.total; GTM_fric_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.u_taus; GTM_us_LaRibhLPveP1400 = gotmkppLaRibhLPveP1400.u_s;
GTM_hL_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.hL; GTM_hLR_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.hLR; GTM_Lat_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.La_t; GTM_tdotus_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.tdotus; GTM_Theatflux_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.total; GTM_fric_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u_taus; GTM_us_LaRibhLPveP1500 = gotmkppLaRibhLPveP1500.u_s; 
# GTM_hL_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.hL; GTM_hLR_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.hLR; GTM_Lat_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.La_t; GTM_tdotus_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.tdotus; GTM_Theatflux_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.total; GTM_fric_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.u_taus; GTM_us_LaRibhLPveP1600 = gotmkppLaRibhLPveP1600.u_s; 
# GTM_hL_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.hL; GTM_hLR_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.hLR; GTM_Lat_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.La_t; GTM_tdotus_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.tdotus; GTM_Theatflux_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.total; GTM_fric_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.u_taus; GTM_us_LaRibhLPveP1750 = gotmkppLaRibhLPveP1750.u_s; 
# GTM_hL_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.hL; GTM_hLR_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.hLR; GTM_Lat_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.La_t; GTM_tdotus_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.tdotus; GTM_Theatflux_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.total; GTM_fric_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.u_taus; GTM_us_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.u_s; 
# GTM_hL_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.hL; GTM_hLR_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.hLR; GTM_Lat_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.La_t; GTM_tdotus_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.tdotus; GTM_Theatflux_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.total; GTM_fric_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.u_taus; GTM_us_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.u_s; 
# GTM_hL_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.hL; GTM_hLR_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.hLR; GTM_Lat_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.La_t; GTM_tdotus_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.tdotus; GTM_Theatflux_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.total; GTM_fric_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.u_taus; 
GTM_hL_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.hL; GTM_hLR_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.hLR; GTM_Lat_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.La_t; GTM_tdotus_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.tdotus; GTM_Theatflux_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.total; GTM_fric_LaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u_taus;




############################################################################
# REDIFINE h/L_L WITH MULTIPLE FACTOR OF 0.0606165
############################################################################
# print(GTM_hLR_LaRib[0:10,:,:])

GTM_hLR_LaRib = GTM_hLR_LaRib[:]*0.0606165

# print(GTM_hLR_LaRib[0:10,:,:])


################################
# extract data for hL analysis #
################################
# GTM_hLR                  = GTM_hLR.isel(lon=0,lat=0)
GTM_hLR_LaRib            = GTM_hLR_LaRib.isel(lon=0,lat=0)
# GTM_hLR_LaRibhL2 = GTM_hLR_LaRibhL2.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP100   = GTM_hLR_LaRibhLPveP100.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1000 = GTM_hLR_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1300 = GTM_hLR_LaRibhLPveP1300.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1350 = GTM_hLR_LaRibhLPveP1350.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1400 = GTM_hLR_LaRibhLPveP1400.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP1500  = GTM_hLR_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1600 = GTM_hLR_LaRibhLPveP1600.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1750 = GTM_hLR_LaRibhLPveP1750.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP2000 = GTM_hLR_LaRibhLPveP2000.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP4000 = GTM_hLR_LaRibhLPveP4000.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP5000 = GTM_hLR_LaRibhLPveP5000.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP10000 = GTM_hLR_LaRibhLPveP10000.isel(lon=0,lat=0)

#######################################
# La_t - extract data for hL analysis #
#######################################

GTM_Lat_LaRib = GTM_Lat_LaRib.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP1000 = GTM_Lat_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP1300 = GTM_Lat_LaRibhLPveP1300.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP1350 = GTM_Lat_LaRibhLPveP1350.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP1400 = GTM_Lat_LaRibhLPveP1400.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP1500 = GTM_Lat_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP1600 = GTM_Lat_LaRibhLPveP1600.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP1750 = GTM_Lat_LaRibhLPveP1750.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP2000 = GTM_Lat_LaRibhLPveP2000.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP4000 = GTM_Lat_LaRibhLPveP4000.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP5000 = GTM_Lat_LaRibhLPveP5000.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP10000 = GTM_Lat_LaRibhLPveP10000.isel(lon=0,lat=0)

#################################################################
# Total heat flux with radiation - extract data for hL analysis #
#################################################################

# GTM_Theatflux_LaRibhLPveP1000 = GTM_Theatflux_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP1300 = GTM_Theatflux_LaRibhLPveP1300.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP1350 = GTM_Theatflux_LaRibhLPveP1350.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP1400 = GTM_Theatflux_LaRibhLPveP1400.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP1500 = GTM_Theatflux_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP1600 = GTM_Theatflux_LaRibhLPveP1600.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP1750 = GTM_Theatflux_LaRibhLPveP1750.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP2000 = GTM_Theatflux_LaRibhLPveP2000.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP4000 = GTM_Theatflux_LaRibhLPveP4000.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP5000 = GTM_Theatflux_LaRibhLPveP5000.isel(lon=0,lat=0)

##########################################################################
# windStress to Stokes drfit flux density - extract data for hL analysis #
##########################################################################

# GTM_tdotus_LaRibhLPveP1000 = GTM_tdotus_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP1300 = GTM_tdotus_LaRibhLPveP1300.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP1350 = GTM_tdotus_LaRibhLPveP1350.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP1400 = GTM_tdotus_LaRibhLPveP1400.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP1500 = GTM_tdotus_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP1600 = GTM_tdotus_LaRibhLPveP1600.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP1750 = GTM_tdotus_LaRibhLPveP1750.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP2000 = GTM_tdotus_LaRibhLPveP2000.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP4000 = GTM_tdotus_LaRibhLPveP4000.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP5000 = GTM_tdotus_LaRibhLPveP5000.isel(lon=0,lat=0)

####################################################
# friction velocity - extract data for hL analysis #
####################################################

# GTM_fric_LaRibhLPveP1000 = GTM_fric_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP1300 = GTM_fric_LaRibhLPveP1300.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP1350 = GTM_fric_LaRibhLPveP1350.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP1400 = GTM_fric_LaRibhLPveP1400.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP1500 = GTM_fric_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP1600 = GTM_fric_LaRibhLPveP1600.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP1750 = GTM_fric_LaRibhLPveP1750.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP2000 = GTM_fric_LaRibhLPveP2000.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP4000 = GTM_fric_LaRibhLPveP4000.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP5000 = GTM_fric_LaRibhLPveP5000.isel(lon=0,lat=0)

########################################################
# Stokes drift velocity - extract data for hL analysis #
########################################################

# GTM_us_LaRibhLPveP1000 = GTM_us_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP1300 = GTM_us_LaRibhLPveP1300.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP1350 = GTM_us_LaRibhLPveP1350.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP1400 = GTM_us_LaRibhLPveP1400.isel(lon=0,lat=0)
GTM_us_LaRibhLPveP1500 = GTM_us_LaRibhLPveP1500.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP1600 = GTM_us_LaRibhLPveP1600.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP1750 = GTM_us_LaRibhLPveP1750.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP2000 = GTM_us_LaRibhLPveP2000.isel(lon=0,lat=0)
# GTM_us_LaRibhLPveP4000 = GTM_us_LaRibhLPveP4000.isel(lon=0,lat=0)

#########################################
# Viscosity - extract data for profiles #
#########################################
# print(GTMnum)

# GTMdepth1 = gotmkpp.z1
# GTMnum = GTMnum.isel(lon=0,lat=0)
# GTMnum = GTMnum.sel(time=1200)
# GTMnum = GTMnum.sel(z1=-2,method='nearest')

# GTMdepth1num2 = gotmkppnum2.z1
# GTMnumnum2 = GTMnumnum2.isel(lon=0,lat=0)
# GTMnumnum2 = GTMnumnum2.sel(time=1200)
# # GTMnumnum2 = GTMnumLanum2.sel(z1=-2,method='nearest')

# GTMdepth1nuh2 = gotmkppnuh2.z1
# GTMnumnuh2 = GTMnumnuh2.isel(lon=0,lat=0)
# GTMnumnuh2 = GTMnumnuh2.sel(time=1200)
# # GTMnumnuh2 = GTMnumnuh2.sel(z1=-2,method='nearest')

# GTMdepth1Ribnum2 = gotmkppRibnum2.z1
# GTMnumRibnum2 = GTMnumRibnum2.isel(lon=0,lat=0)
# GTMnumRibnum2 = GTMnumRibnum2.sel(time=1200)

# GTMdepth1Ribnuh2 = gotmkppRibnuh2.z1
# GTMnumRibnuh2 = GTMnumRibnuh2.isel(lon=0,lat=0)
# GTMnumRibnuh2 = GTMnumRibnuh2.sel(time=1200)
# # GTMnumRibnuh2 = GTMnumRibnuh2.sel(z1=-2,method='nearest')

# print('depth : ',GTMdepth1num2)
# print('viscosity : ',GTMnumnum2)


###########################################
# Diffusivity - extract data for profiles #
###########################################

# GTMnuh = GTMnuh.isel(lon=0,lat=0)
# GTMnuh = GTMnuh.sel(time=1200)
# GTMnuh = GTMnuh.sel(z1=-2,method='nearest')

# GTMnuhnum2 = GTMnuhnum2.isel(lon=0,lat=0)
# GTMnuhnum2 = GTMnuhnum2.sel(time=1200)
# # GTMnuhnum2 = GTMnuhnum2.sel(z1=-2,method='nearest')

# GTMnuhnuh2 = GTMnuhnuh2.isel(lon=0,lat=0)
# GTMnuhnuh2 = GTMnuhnuh2.sel(time=1200)
# # GTMnuhnuh2 = GTMnuhnuh2.sel(z1=-2,method='nearest')

# GTMnuhRibnum2 = GTMnuhRibnum2.isel(lon=0,lat=0)
# GTMnuhRibnum2 = GTMnuhRibnum2.sel(time=1200)
# # GTMnuhRibnum2 = GTMnuhRibnum2.sel(z1=-2,method='nearest')

# GTMnuhRibnuh2 = GTMnuhRibnuh2.isel(lon=0,lat=0)
# GTMnuhRibnuh2 = GTMnuhRibnuh2.sel(time=1200)
# # GTMnuhRibnuh2 = GTMnuhRibnuh2.sel(z1=-2,method='nearest')


# # apply moving average for 1 week
# time2LaRibhL2_1015 = moving_average(time2LaRibhL2,1015)
# GTM_hL_LaRibhL2_1015 = moving_average(GTM_hL_LaRibhL2.isel(lon=0,lat=0),1015)
# GTM_hLR_LaRibhL2_1015 = moving_average(GTM_hLR_LaRibhL2,1015)

#############################################################################################################################################
#---------------------------------------------#
#          hL parameter -- SMA=2wks           #
#---------------------------------------------#
# # hL2_test
# time2LaRibhL2_2030 = moving_average(time2LaRibhL2,2030)
# GTM_hL_LaRibhL2_2030 = moving_average(GTM_hL_LaRibhL2.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhL2_2030 = moving_average(GTM_hLR_LaRibhL2,2030)
# # Original KPP model
# time2_2030 = moving_average(time2,2030)
# GTM_hL_2030 = moving_average(GTM_hL.isel(lon=0,lat=0),2030)
# GTM_hLR_2030 = moving_average(GTM_hLR,2030)
# nu_LaT_Rib Ch=Cm=1
time2LaRib_2030 = moving_average(time2LaRib,2030)
GTM_hL_LaRib_2030 = moving_average(GTM_hL_LaRib.isel(lon=0,lat=0),2030)
GTM_hLR_LaRib_2030 = moving_average(GTM_hLR_LaRib,2030)
# # nu_LaT_Rib Ch=Cm=1 laT>0.3
# time2LaRibLaTPve03_2030 = moving_average(time2LaRibLaTPve03,2030)
# GTM_hL_LaRibLaTPve03_2030 = moving_average(GTM_hL_LaRibLaTPve03.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibLaTPve03_2030 = moving_average(GTM_hLR_LaRibLaTPve03,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>0
# time2LaRibhLPve_2030 = moving_average(time2LaRibhLPve,2030)
# GTM_hL_LaRibhLPve_2030 = moving_average(GTM_hL_LaRibhLPve.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPve_2030 = moving_average(GTM_hLR_LaRibhLPve,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>0
# time2LaRibhLgamPve_2030 = moving_average(time2LaRibhLgamPve,2030)
# GTM_hL_LaRibhLgamPve_2030 = moving_average(GTM_hL_LaRibhLgamPve.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLgamPve_2030 = moving_average(GTM_hLR_LaRibhLgamPve,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>100
# time2LaRibhLPveP100_2030 = moving_average(time2LaRibhLPveP100,2030)
# GTM_hL_LaRibhLPveP100_2030 = moving_average(GTM_hL_LaRibhLPveP100.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP100_2030 = moving_average(GTM_hLR_LaRibhLPveP100,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>100
# time2LaRibhLgamPveP100_2030 = moving_average(time2LaRibhLgamPveP100,2030)
# GTM_hL_LaRibhLgamPveP100_2030 = moving_average(GTM_hL_LaRibhLgamPveP100.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLgamPveP100_2030 = moving_average(GTM_hLR_LaRibhLgamPveP100,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL<0
# time2LaRibhLNve_2030 = moving_average(time2LaRibhLNve,2030)
# GTM_hL_LaRibhLNve_2030 = moving_average(GTM_hL_LaRibhLNve.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLNve_2030 = moving_average(GTM_hLR_LaRibhLNve,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>-100
# time2LaRibhLPveN100_2030 = moving_average(time2LaRibhLPveN100,2030)
# GTM_hL_LaRibhLPveN100_2030 = moving_average(GTM_hL_LaRibhLPveN100.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveN100_2030 = moving_average(GTM_hLR_LaRibhLPveN100,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>-100
# time2LaRibhLgamPveN100_2030 = moving_average(time2LaRibhLgamPveN100,2030)
# GTM_hL_LaRibhLgamPveN100_2030 = moving_average(GTM_hL_LaRibhLgamPveN100.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLgamPveN100_2030 = moving_average(GTM_hLR_LaRibhLgamPveN100,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>-200
# time2LaRibhLPveN200_2030 = moving_average(time2LaRibhLPveN200,2030)
# GTM_hL_LaRibhLPveN200_2030 = moving_average(GTM_hL_LaRibhLPveN200.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveN200_2030 = moving_average(GTM_hLR_LaRibhLPveN200,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>-500
# time2LaRibhLPveN500_2030 = moving_average(time2LaRibhLPveN500,2030)
# GTM_hL_LaRibhLPveN500_2030 = moving_average(GTM_hL_LaRibhLPveN500.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveN500_2030 = moving_average(GTM_hLR_LaRibhLPveN500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>-1000
# time2LaRibhLPveN1000_2030 = moving_average(time2LaRibhLPveN1000,2030)
# GTM_hL_LaRibhLPveN1000_2030 = moving_average(GTM_hL_LaRibhLPveN1000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveN1000_2030 = moving_average(GTM_hLR_LaRibhLPveN1000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1000
# time2LaRibhLPveP1000_2030 = moving_average(time2LaRibhLPveP1000,2030)
# GTM_hL_LaRibhLPveP1000_2030 = moving_average(GTM_hL_LaRibhLPveP1000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1000_2030 = moving_average(GTM_hLR_LaRibhLPveP1000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
# time2LaRibhLPveP1300_2030 = moving_average(time2LaRibhLPveP1300,2030)
# GTM_hL_LaRibhLPveP1300_2030 = moving_average(GTM_hL_LaRibhLPveP1300.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1300_2030 = moving_average(GTM_hLR_LaRibhLPveP1300,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1350
# time2LaRibhLPveP1350_2030 = moving_average(time2LaRibhLPveP1350,2030)
# GTM_hL_LaRibhLPveP1350_2030 = moving_average(GTM_hL_LaRibhLPveP1350.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1350_2030 = moving_average(GTM_hLR_LaRibhLPveP1350,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3250
# time2LaRibhLPveP1400_2030 = moving_average(time2LaRibhLPveP1400,2030)
# GTM_hL_LaRibhLPveP1400_2030 = moving_average(GTM_hL_LaRibhLPveP1400.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1400_2030 = moving_average(GTM_hLR_LaRibhLPveP1400,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1500
time2LaRibhLPveP1500_2030 = moving_average(time2LaRibhLPveP1500,2030)
GTM_hL_LaRibhLPveP1500_2030 = moving_average(GTM_hL_LaRibhLPveP1500.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP1500_2030 = moving_average(GTM_hLR_LaRibhLPveP1500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
# time2LaRibhLPveP1600_2030 = moving_average(time2LaRibhLPveP1600,2030)
# GTM_hL_LaRibhLPveP1600_2030 = moving_average(GTM_hL_LaRibhLPveP1600.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1600_2030 = moving_average(GTM_hLR_LaRibhLPveP1600,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>500
# time2LaRibhLPveP1750_2030 = moving_average(time2LaRibhLPveP1750,2030)
# GTM_hL_LaRibhLPveP1750_2030 = moving_average(GTM_hL_LaRibhLPveP1750.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1750_2030 = moving_average(GTM_hLR_LaRibhLPveP1750,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>2000
# time2LaRibhLPveP2000_2030 = moving_average(time2LaRibhLPveP2000,2030)
# GTM_hL_LaRibhLPveP2000_2030 = moving_average(GTM_hL_LaRibhLPveP2000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP2000_2030 = moving_average(GTM_hLR_LaRibhLPveP2000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
# time2LaRibhLPveP4000_2030 = moving_average(time2LaRibhLPveP4000,2030)
# GTM_hL_LaRibhLPveP4000_2030 = moving_average(GTM_hL_LaRibhLPveP4000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP4000_2030 = moving_average(GTM_hLR_LaRibhLPveP4000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>5000
# time2LaRibhLPveP5000_2030 = moving_average(time2LaRibhLPveP5000,2030)
# GTM_hL_LaRibhLPveP5000_2030 = moving_average(GTM_hL_LaRibhLPveP5000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP5000_2030 = moving_average(GTM_hLR_LaRibhLPveP5000,2030)



print('ERA orig start time: ',mdate[0],mtime[0])
print('GOTM orig start timeL ', timeLaRib[0])
print('GOTM start time',time2LaRib[0])
print('GOTM 2wks start time',time2LaRib_2030[0])
print('GOTM datetime time',datetime2[0])


#---------------------------------------------#
#          LaT parameter -- SMA=2wks          #
#---------------------------------------------#
# # Original KPP
# GTM_Lat_2030 = moving_average(GTM_Lat,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
GTM_Lat_LaRib_2030 = moving_average(GTM_Lat_LaRib,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
# GTM_Lat_LaRibhLPveP1000_2030 = moving_average(GTM_Lat_LaRibhLPveP1000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
# GTM_Lat_LaRibhLPveP1300_2030 = moving_average(GTM_Lat_LaRibhLPveP1300,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1350
# GTM_Lat_LaRibhLPveP1350_2030 = moving_average(GTM_Lat_LaRibhLPveP1350,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3250
# GTM_Lat_LaRibhLPveP1400_2030 = moving_average(GTM_Lat_LaRibhLPveP1400,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1500
GTM_Lat_LaRibhLPveP1500_2030 = moving_average(GTM_Lat_LaRibhLPveP1500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
# GTM_Lat_LaRibhLPveP1600_2030 = moving_average(GTM_Lat_LaRibhLPveP1600,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3750
# GTM_Lat_LaRibhLPveP1750_2030 = moving_average(GTM_Lat_LaRibhLPveP1750,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>2000
# GTM_Lat_LaRibhLPveP2000_2030 = moving_average(GTM_Lat_LaRibhLPveP2000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
# GTM_Lat_LaRibhLPveP4000_2030 = moving_average(GTM_Lat_LaRibhLPveP4000,2030)

#--------------------------------------------------------#
#          Heat flux with radiation -- SMA=2wks          #
#--------------------------------------------------------#

# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1000
# GTM_Theatflux_LaRibhLPveP1000_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
# GTM_Theatflux_LaRibhLPveP1300_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1300,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1350
# GTM_Theatflux_LaRibhLPveP1350_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1350,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3250
# GTM_Theatflux_LaRibhLPveP1400_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1400,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1500
GTM_Theatflux_LaRibhLPveP1500_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1600
# GTM_Theatflux_LaRibhLPveP1600_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1600,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3750
# GTM_Theatflux_LaRibhLPveP1750_2030 = moving_average(GTM_Theatflux_LaRibhLPveP1750,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>2000
# GTM_Theatflux_LaRibhLPveP2000_2030 = moving_average(GTM_Theatflux_LaRibhLPveP2000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
# GTM_Theatflux_LaRibhLPveP4000_2030 = moving_average(GTM_Theatflux_LaRibhLPveP4000,2030)

#-------------------------------------------------#
#          friction velocity -- SMA=2wks          #
#-------------------------------------------------#
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
# GTM_fric_LaRibhLPveP2000_2030 = moving_average(GTM_fric_LaRibhLPveP2000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1500
GTM_fric_LaRibhLPveP1500_2030 = moving_average(GTM_fric_LaRibhLPveP1500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1600
# GTM_fric_LaRibhLPveP1600_2030 = moving_average(GTM_fric_LaRibhLPveP1600,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
# GTM_fric_LaRibhLPveP1300_2030 = moving_average(GTM_fric_LaRibhLPveP1300,2030)
# # # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
# GTM_fric_LaRibhLPveP4000_2030 = moving_average(GTM_fric_LaRibhLPveP4000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1350
# GTM_fric_LaRibhLPveP1350_2030 = moving_average(GTM_fric_LaRibhLPveP1350,2030)

#-----------------------------------------------------#
#          Stokes drift velocity -- SMA=2wks          #
#-----------------------------------------------------#
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
# GTM_us_LaRibhLPveP2000_2030 = moving_average(GTM_us_LaRibhLPveP2000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1500
GTM_us_LaRibhLPveP1500_2030 = moving_average(GTM_us_LaRibhLPveP1500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1600
# GTM_us_LaRibhLPveP1600_2030 = moving_average(GTM_us_LaRibhLPveP1600,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1300
# GTM_us_LaRibhLPveP1300_2030 = moving_average(GTM_us_LaRibhLPveP1300,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
# GTM_us_LaRibhLPveP4000_2030 = moving_average(GTM_us_LaRibhLPveP4000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1350
# GTM_us_LaRibhLPveP1350_2030 = moving_average(GTM_us_LaRibhLPveP1350,2030)

#-----------------------------------------------------#
#          Surface sea temperature - degreeC          #
#-----------------------------------------------------#

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
### Plot of MLD from OSMOSIS Observations
# plt.plot(j,np.negative(i),color='olivedrab',linestyle='dotted'); 
# plt.xlabel('Time -Julian day'); plt.ylabel('Mixed Layer Depth -m')
# plt.grid()
# # plt.axis([260,341,-400,0]); 
# # plt.title('OSMOSIS Mixed Layer Depth')
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


# print('manualy moving average :', df_mld)
# df_mld = pd.read_csv('mld.csv')
# print(df_mld.info())
# print(df_mld)
# df_mld = pd.DataFrame(i)
# df_mld = i.rolling(3,min_periods=1).mean()
# print('moving average MLD', df_mld)

# #################
# ## 2D meshplot ##
# #################
# print('GTM variables : ', gotmkpp.variables.items() )
# # plt.pcolormesh(GTMnuhLaRibhLPveP10000,time2LaRibhLPveP10000,GTMdepthLaRibhLPveP10000)

# depthT, timeT = np.meshgrid(GTMdepthLaRibhLPveP10000,time2LaRibhLPveP10000)
# temp = GTMtempLaRibhLPveP10000
# plt.contour(depthT,timeT,temp,20,cmap='inferno_r')
# plt.colorbar()
# plt.xlabel(r'Depth -$m$'); plt.ylabel('Time - Julian days')
# plt.title(r'Temperature contour plot $h/L >10 000$')
# plt.figure()

# depthT, timeT = np.meshgrid(GTMdepth,time2)
# temp = GTMtemp
# plt.contour(depthT,timeT,temp,20,cmap='inferno_r')
# plt.colorbar()
# plt.xlabel(r'Depth -$m$'); plt.ylabel('Time - Julian days')
# plt.title(r'Temperature contour plot')
# plt.figure()


# GTMnumLaRibhLPveP10000 = GTMnumLaRibhLPveP10000.isel(lon=0,lat=0)
# depthM, timeM = np.meshgrid(GTMdepthLaRibhLPveP10000,time2LaRibhLPveP10000)
# nuM = GTMnumLaRibhLPveP10000
# plt.contour(depthM,timeM,nuM,20,cmap='inferno_r')
# plt.colorbar()
# plt.xlabel(r'Depth -$m$'); plt.ylabel('Time - Julian days')
# plt.title(r'Viscosity contour plot $h/L >10 000$')
# plt.figure()

# depthM, timeM = np.meshgrid(GTMdepth,time2)
# nuM = GTMnum.isel(lat=0,lon=0)
# plt.contour(depthM,timeM,nuM,20,cmap='inferno_r')
# plt.colorbar()
# plt.xlabel(r'Depth -$m$'); plt.ylabel('Time - Julian days')
# plt.title(r'Viscosity contour plot')
# plt.figure()

# GTMnuhLaRibhLPveP10000 = GTMnuhLaRibhLPveP10000.isel(lon=0,lat=0)
# depthH, timeH = np.meshgrid(GTMdepthLaRibhLPveP10000,time2LaRibhLPveP10000)
# nuH = GTMnuhLaRibhLPveP10000
# plt.contour(depthH,timeH,nuH,20,cmap='inferno_r')
# plt.colorbar()
# plt.xlabel(r'Depth -$m$'); plt.ylabel('Time - Julian days')
# plt.title(r'Diffusivity contour plot $h/L >10 000$')
# plt.figure()

# depthH, timeH = np.meshgrid(GTMdepth,time2)
# nuH = GTMnuh.isel(lat=0,lon=0)
# plt.contour(depthH,timeH,nuH,20,cmap='inferno_r')
# plt.colorbar()
# plt.xlabel(r'Depth -$m$'); plt.ylabel('Time - Julian days')
# plt.title(r'Diffusivity contour plot')
# plt.figure()

# plt.show()

# exit()
# ###################################################
# ## Comparison of MLDs with various hL thresholds ##
# ###################################################

# fig, ax = plt.subplots(2, 2)

# ax[0, 0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[0, 0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[0, 0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# # ax[0, 0].plot(time2LaRibhLPve,mldKPPLaRibhLPve,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > 0') 
# ax[0, 0].grid(); ax[0, 0].legend(loc='best',fontsize='small')
# # ax[0, 0].set_xlim(490,635); ax[0, 0].set_ylim(-150,0)

# ax[0, 1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[0, 1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[0, 1].plot(time2LaRibnum2,mldKPPLaRibnum2,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 num*2')
# # ax[0, 1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# # ax[0, 1].plot(time2LaRibhLPveP100,mldKPPLaRibhLPveP100,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > 100')
# ax[0, 1].grid(); ax[0, 1].legend(loc='best',fontsize='small')
# # ax[0, 1].set_xlim(490,635); ax[0, 1].set_ylim(-150,0)

# ax[1, 0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5') 
# ax[1, 0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[1, 0].plot(time2LaRibnuh2,mldKPPLaRibnuh2,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 nuh*2')
# # ax[1, 0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# # ax[1, 0].plot(time2LaRibhLPveP500,mldKPPLaRibhLPveP500,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > 500')
# ax[1, 0].grid(); ax[1, 0].legend(loc='best',fontsize='small')
# # ax[1, 0].set_xlim(490,635); ax[1, 0].set_ylim(-150,0)

# ax[1, 1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[1, 1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[1, 1].plot(time2LaRibnum2nuh2,mld_tempKPPLaRibnum2nuh2,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 num*2 nuh*2')
# # ax[1, 1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# # ax[1, 1].plot(time2LaRibhLPveP1000,mldKPPLaRibhLPveP1000,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > 1000')
# ax[1, 1].grid(); ax[1, 1].legend(loc='best',fontsize='small')
# # ax[1, 1].set_xlim(490,635); ax[1, 1].set_ylim(-150,0)


# # fig.legend((l1, l2), ('Line 1', 'Line 2'), 'upper left')
# # fig.legend((l3, l4), ('Line 3', 'Line 4'), 'upper right')

# fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0, 0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0,0].transAxes)
# ax[0, 1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[0,1].transAxes)
# ax[1, 0].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[1,0].transAxes)
# ax[1, 1].text(0.06, 0.08, '(d)', horizontalalignment='center', verticalalignment='center', transform=ax[1,1].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs with various hL thresholds')
# plt.show()
# exit()

# ###################################################################################################################
# ## Compare evolution of MLDs with NON-LOCAL MIXING vs LOCAL MIXING using hL thresholds ranging [3000,3500,4000]  ## 
# ###################################################################################################################

# fig, ax = plt.subplots(2,1)

# ax[0].plot(df_dats5,np.negative(df_mld5), color='maroon',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
# ax[0].plot(time2LaRibhLPveP4000,mldKPPLaRibhLPveP4000,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
# ax[0].set_ylabel('MLD (NON-LOCAL enabled) -m');
# ax[0].set_xlim(265,650)
# # ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# ax[1].plot(df_dats5,np.negative(df_mld5), color='maroon',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
# ax[1].plot(time2LaRibhLPveP1350,mldKPPLaRibhLPveP4000_Local,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
# ax[1].set_ylabel('MLD (LOCAL) -m');
# ax[1].set_xlim(265,650)
# # ax[1].set_xlim(480,650)
# # ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# # fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
# ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# # plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
# plt.show()

# exit()

# ###############################################################################
# # Comparison of MLDs using hL thresholds w NON-LOCAL mixing  ##
# ###############################################################################

# compare_MLD = np.c_[time2,mldKPP,mldKPPLaRibhLPveP10000]

# # with open('compare_MLD.dat','w') as f: 
# #     for row in compare_MLD:
# #         f.write(row)

# fig, ax = plt.subplots(3,1,sharex=True,figsize=(40,18))

# # # testing factors nuh num
# # ax.plot(df_dats5,np.negative(df_mld5), color='crimson',linestyle='solid',label = 'OSMOSIS')
# # ax.plot(time2,mldKPP,color='black',linestyle='solid',label = r'KPP $La_{T}$ OFF $h/L_l$ OFF')
# # ax.plot(time2num2,mldKPPnum2,color='dodgerblue',linestyle='solid',label = r'KPP $\nu_M*1.5$ $La_{T}$ OFF')
# # # ax.plot(time2Rib05num2,mldKPPRib05num2,color='darkmagenta',linestyle='solid',label = r'KPP $\nu_M*1.5 \ Ri_{B}/1.5 \ La_{T}$ OFF')
# # ax.plot(time2nuh2,mldKPPnuh2,color='limegreen',linestyle='solid',label = r'KPP $\nu_H*1.5$ $La_{T}$ OFF')
# # # ax.plot(time2Rib2nuh2,mldKPPRib2nuh2,color='chocolate',linestyle='solid',label = r'KPP $\nu_H*1.5 \ Ri_{B}*1.5 \ La_{T}$ OFF')
# # ax.grid(); ax.legend(loc='lower right',fontsize='small')
# # # ax.set_xlim(250,650);# ax.set_ylim(-200,0)
# # # ax.set_xlim(480,540); ax.set_ylim(-200,0)

# ax0=ax[0].twiny(); 
# ax0.tick_params(axis='both', which='major', labelsize=18)
# ax0.set_axisbelow(True)
# # ax1.set_xlim(mdate[0],mdate[-1])
# xa0 = ax0.get_xaxis(); 
# xa0.axis_date(); 
# xa0.set_minor_locator(mdates.DayLocator(interval=1)); 
# xa0.set_major_locator(mdates.MonthLocator()); 
# xa0.set_major_formatter(mdates.DateFormatter('%b-%y'))
# ax0.plot(new_time2,GldrTemp,color='white',linestyle='dotted')

# # # SST
# # ax.plot(new_time,GldrTemp,color='indigo',linestyle='dashed',label ='OSMOSIS')
# # ax.plot(time2,GTMtemp,color='black',linestyle='solid',label =r'KPP $La_{T}$ OFF $h/L_L$ OFF')
# # ax.plot(time2,GTMtempLaRib,color='goldenrod',linestyle='solid',label =r'KPP $La_{T}$ $\gamma = \delta = 1$; $h/L_L$ OFF')
# # ax.plot(time2LaRibhLPveP100,GTMtempLaRibhLPveP100,color='limegreen',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L_L>6.10$')
# # ax.plot(time2LaRibhLPveP1500,GTMtempLaRibhLPveP1500,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L_L>90.90$')
# # ax.plot(time2LaRibhLPveP10000,GTMtempLaRibhLPveP10000,color='dodgerblue',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L_L>606.20$')
# # ax.grid(); ax.legend(loc='lower right',fontsize=25)
# # ax.tick_params(axis='both', which='major', labelsize=20)
# # # ax.set_xlim(400,635);# ax.set_ylim(-200,0)
# # # ax.set_xlim(480,540);# ax.set_ylim(-15,0)
# # # ax[1].set_ylim(11,23)
# # ax.set_ylabel(r'SST -$\degree C$',fontsize=30)
# # ax.set_xlabel(r'Time -Julian days',fontsize=30)

# # # Isothermal Layer Depth
# # ax.plot(new_time,ild_Gldr, color='black',linestyle='solid',label = 'OSMOSIS')
# # ax.grid(); ax.legend(loc='lower right',fontsize='small')
# # # # ax.set_xlim(480,540);# ax.set_ylim(-150,0)
# # ax.set_ylabel(r'ILD -$m$')


# ax[0].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = 'OSMOSIS')
# ax[0].plot(time2,mldKPP,color='black',linestyle='solid',label = r'KPP $La_{T}$ OFF $h/L_L$ OFF')
# ax[0].plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'KPP $La_{T}$ $\gamma = \delta = 1$; $h/L_L$ OFF')
# ax[0].plot(time2LaRibhLPveP100,mldKPPLaRibhLPveP100,color='limegreen',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1; \ \delta=0 \ for \ h/L_L>6.10$')
# ax[0].grid(); ax[0].legend(loc='lower right',fontsize=15)
# ax[0].tick_params(axis='both', which='major', labelsize=20)
# # ax[0].set_xlim(250,650);# ax[0].set_ylim(-200,0)
# ax[0].set_xlim(480,550); ax[0].set_ylim(-150,0)

# ax[1].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = '_nolegend_')
# ax[1].plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
# ax[1].plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'_nolegend_')
# ax[1].plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L_L>90.90$')
# ax[1].grid(); ax[1].legend(loc='lower right',fontsize=15)
# ax[1].set_ylabel('Mixed Layer Depth -m',fontsize=25)
# ax[1].tick_params(axis='both', which='major', labelsize=20)
# # ax[1].set_xlim(250,650);# ax[1].set_ylim(-200,0)
# ax[1].set_xlim(480,550); ax[1].set_ylim(-150,0)

# ax[2].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = '_nolegend_')
# ax[2].plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
# ax[2].plot(time2LaRib,mldKPPLaRib,color='goldenrod',linestyle='solid',label = r'_nolegend_')
# ax[2].plot(time2LaRibhLPveP10000,mldKPPLaRibhLPveP10000,color='dodgerblue',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L_L>606.20$')
# ax[2].grid(); ax[2].legend(loc='lower right',fontsize=15)
# ax[2].set_xlabel('Time -Julian days',fontsize=25)
# ax[2].tick_params(axis='both', which='major', labelsize=20)
# # ax[2].set_xlim(250,650);# ax[2].set_ylim(-200,0)
# ax[2].set_xlim(480,550); ax[2].set_ylim(-150,0)

# # ## test conditions

# # ax[0].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='dashed',label = 'OSMOSIS')
# # ax[0].plot(time2,mldKPP,color='black',linestyle='solid',label = r'KPP $La_{T}$ OFF $h/L_L$ OFF')
# # # ax[0].plot(time2LaRib,mldKPPLaRib,color='dodgerblue',linestyle='solid',label = r'_nolegend_')
# # # ax[0].plot(time2LaRibhLPveP100,mldKPPLaRibhLPveP100,color='limegreen',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1; \ \delta=0 \ for \ h/L>100$')
# # # ax[0].plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='dodgerblue',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L>1500$')
# # ax[0].plot(time2LaRibhLPveP10000,mldKPPLaRibhLPveP10000,color='orangered',linestyle='solid',label = r'KPP $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L>10 000$')
# # ax[0].grid(); ax[0].legend(loc='lower right',fontsize='small')
# # # ax[0].set_xlim(250,650);# ax[0].set_ylim(-200,0)
# # # ax[0].set_xlim(480,600); ax[0].set_ylim(-150,0)
# # ax[0].set_ylabel(r'Mixed Layer Depth -$m$')

# # ax[1].plot(df_dats5,np.negative(df_mld5), color='indigo',linestyle='solid',label = '_nolegend_')
# # ax[1].plot(time2,mldKPP,color='black',linestyle='solid',label = '_nolegend_')
# # # ax[1].plot(time2LaRib,mldKPPLaRib,color='dodgerblue',linestyle='solid',label = r'_nolegend_')
# # # ax[1].plot(time2LaRibhLPveP10000shear,mldKPPLaRibhLPveP10000shear,color='orangered',linestyle='solid',label = r'KPP $\nu_S$ $La_{T}$ $\gamma=1 \ \delta=0 \ for \ h/L>10 000$')
# # ax[1].grid(); ax[1].legend(loc='lower right',fontsize='small')
# # # ax[1].set_xlim(250,650);# ax[1].set_ylim(-200,0)
# # # ax[1].set_xlim(480,600); ax[1].set_ylim(-150,0)

# # plt.legend(loc='best',fontsize='small')

# # fig.text(0.5, 0.01, r'Time -Julian days', ha='center', va='center')
# # # fig.text(0.008, 0.5, r'Mixed Layer Depth -$m$', ha='center', va='center', rotation='vertical')
# # # # # ax.text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
# ax[0].text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes,fontsize=20)
# ax[1].text(0.02, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes,fontsize=20)
# ax[2].text(0.02, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes,fontsize=20)

# # plt.xlabel('Time -Julian days',fontsize=12); plt.ylabel('Mixed Layer Depth -m',fontsize=12)

# # # # plt.legend(loc='best',fontsize='small')
# # # # plt.suptitle(r'Comparison of mixed layer depths with $La_T$ and $h/L$ thresholds')
# plt.show()
# plt.savefig('compare_MLD_thresholds_480_600.png')



# #######################################################
# ## Comparison of viscosity and diffusivity profiles  ##
# #######################################################

# fig, ax = plt.subplots(1,2)

# ax[0].plot(GTMnum,GTMdepth1,color='black',linestyle='solid',label = r'_nolegend_')
# ax[0].plot(GTMnumLanum4,GTMdepth1Lanum4,color='limegreen',linestyle='solid',label = r'_nolegend_')
# ax[0].plot(GTMnumLanuh4,GTMdepth1Lanuh4,color='saddlebrown',linestyle='solid',label = r'_nolegend_')
# ax[0].plot(GTMnumLanum2,GTMdepth1Lanum2,color='dodgerblue',linestyle='solid',label = r'_nolegend_')
# ax[0].plot(GTMnumLanuh2,GTMdepth1Lanuh2,color='orangered',linestyle='solid',label = r'_nolegend_')
# ax[0].grid(); ax[0].legend(loc='lower right',fontsize='small')
# ax[0].set_xlabel(r'$\nu_M$')
# ax[0].set_ylabel('Depth -m')
# # ax[0].set_xlim(475,635); 
# ax[0].set_ylim(-50,0)

# ax[1].plot(GTMnuh,GTMdepth1,color='black',linestyle='solid',label = r'KPP $La_{T}$ OFF $h/L$ OFF')
# ax[1].plot(GTMnuhLanum4,GTMdepth1Lanum4,color='limegreen',linestyle='solid',label = r'KPP $\nu_M*0.5$ $La_{T}$ OFF')
# ax[1].plot(GTMnuhLanuh4,GTMdepth1Lanuh4,color='saddlebrown',linestyle='solid',label = r'KPP $\nu_H*0.5$ $La_{T}$ OFF')
# ax[1].plot(GTMnuhLanum2,GTMdepth1Lanum2,color='dodgerblue',linestyle='solid',label = r'KPP $\nu_M*2$; $Ri_{B}$ $La_{T}$ OFF')
# ax[1].plot(GTMnuhLanuh2,GTMdepth1Lanuh2,color='orangered',linestyle='solid',label = r'KPP $\nu_H*2$; $Ri_{B}$ $La_{T}$ OFF')
# ax[1].grid(); ax[1].legend(loc='lower right',fontsize='small')
# ax[1].set_xlabel(r'$\nu_H$')
# # ax[1].set_ylabel('Depth -m')
# # ax[1].set_xlim(475,635); 
# ax[1].set_ylim(-50,0)


# # fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# # fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0].text(0.06, 0.08, '(a)', horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes)
# ax[1].text(0.06, 0.08, '(b)', horizontalalignment='left', verticalalignment='center', transform=ax[1].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# # plt.suptitle(r'The profiles of eddy viscosity $(\nu_M)$ and heat diffusivity $(\nu_H)$')
# # plt.suptitle(r'The profiles of eddy viscosity $(\nu_M)$ and heat diffusivity $(\nu_H)$')


# plt.show()



# exit()
# #############################################################################
# ## Comparison of MLDs using hL thresholds w LOCAL mixing [3000,3500,4000]  ##
# #############################################################################

# fig, ax = plt.subplots(3,1)

# ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid')#,label = 'OSMOSIS moving average n=5')
# ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid')#,label = 'GOTM KPP')
# ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid')#,label = 'KPP Lat Rib Cm=Ch=1')
# ax[0].plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='maroon',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>1500$')
# ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
# ax[0].set_xlim(475,635); ax[0].set_ylim(-150,0)

# ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid')#,label = 'OSMOSIS moving average n=5')
# ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid')#,label = 'GOTM KPP')
# ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid')#,label = 'KPP Lat Rib Cm=Ch=1')
# ax[1].plot(time2LaRibhLPveP1300,mldKPPLaRibhLPveP1600_Local,color='maroon',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
# ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
# ax[1].set_xlim(475,635); ax[1].set_ylim(-150,0)

# ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
# ax[2].plot(time2LaRibhLPveP1350,mldKPPLaRibhLPveP4000_Local,color='maroon',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
# ax[2].set_xlim(475,635); ax[2].set_ylim(-150,0)

# fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0].text(0.06, 0.08, '(a)', horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes)
# ax[1].text(0.06, 0.08, '(b)', horizontalalignment='left', verticalalignment='center', transform=ax[1].transAxes)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='left', verticalalignment='center', transform=ax[2].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds w LOCAL mixing [3000,3500,4000]')
# plt.show()

# exit()

###############################
## h/L vs LaT regime diagram ##
###############################
## regions in the regime diagram
Langmuir_x = []; Langmuir_y = [];
with open('../regimeLangmuir.dat', 'r') as f:
	for row in f:
		b, c = row.split()
		Langmuir_x.append(b)
		Langmuir_y.append(c)
Langmuir_x = np.asarray(Langmuir_x); Langmuir_x = Langmuir_x.astype('float64')
Langmuir_y = np.asarray(Langmuir_y); Langmuir_y = Langmuir_y.astype('float64')

Wind_x = []; Wind_y = [];
with open('../regimeWind.dat', 'r') as f:
	for row in f:
		d, e = row.split()
		Wind_x.append(d)
		Wind_y.append(e)
Wind_x = np.asarray(Wind_x); Wind_x = Wind_x.astype('float64')
Wind_y = np.asarray(Wind_y); Wind_y = Wind_y.astype('float64')

Convection_x = []; Convection_y = [];
with open('../regimeConvection.dat', 'r') as f:
	for row in f:
		g, h = row.split()
		Convection_x.append(g)
		Convection_y.append(h)
Convection_x = np.asarray(Convection_x); Convection_x = Convection_x.astype('float64')
Convection_y = np.asarray(Convection_y); Convection_y = Convection_y.astype('float64')

# ## scatter plots
# # with radiation
# fig, ax = plt.subplots(1,2, figsize=(30,15))
# # plt.rcParams['figure.figsize'] = [20,10]
# ax[0].scatter(GTM_Lat_LaRib,GTM_hLR_LaRib,s=1,color='orangered')
# ax[0].plot(Langmuir_x, Langmuir_y, 'k'); ax[0].plot(Wind_x, Wind_y, 'k'); ax[0].plot(Convection_x,Convection_y, 'k')
# ax[0].text(0.15, 0.005, 'LANGMUIR', ha='center', va='center',fontsize=25)
# ax[0].text(0.20, 10000, 'CONVECTION', ha='center', va='center',fontsize=25)
# ax[0].text(6, 0.005, 'WIND', ha='center', va='center',fontsize=25)
# ax[0].set_ylabel(r'$h/L_L$',fontsize=30); ax[0].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
# ax[0].set_xlabel(r'$La_T$',fontsize=30); ax[0].set_xscale('log')
# ax[0].set_ylim([0.001,1000000])
# ax[0].set_xlim([0.1,10])
# ax[0].tick_params(axis='both', which='major', labelsize=20)
# ax[0].grid()
# # ax[0].axis([0.1,10,0.001,1000000])
# # ax[0].set_yticks([0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000])


# ax[1].scatter(GTM_Lat_LaRib,np.negative(GTM_hLR_LaRib),s=1,color='orangered')
# ax[1].set_ylabel(r'-$h/L_L$',fontsize=30); ax[1].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
# ax[1].set_xlabel(r'$La_T$',fontsize=30); ax[1].set_xscale('log')
# ax[1].set_ylim([0.001,1000000])
# ax[1].set_xlim([0.1,10])
# ax[1].tick_params(axis='both', which='major', labelsize=20)
# ax[1].grid()
# # ax[0].axis([0.1,10,0.001,1000000])
# # ax[1].get_yticks([0.01,0,100,10000,1000000])
# # ax[1].get_yaxis().set_major_formatter(LogFormatter())

# ax[0].text(0.03, 0.95, '(a)', horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes, fontsize=20)
# ax[1].text(0.03, 0.95, '(b)', horizontalalignment='left', verticalalignment='center', transform=ax[1].transAxes, fontsize=20)

# plt.savefig("regime_diagram_w_solar_rad.png")
# # # without radiation
# # fig, ax = plt.subplots(1,2, figsize=(30,15))
# # # plt.rcParams['figure.figsize'] = [20,10]
# # ax[0].scatter(GTM_Lat_LaRibhLPveP1500,GTM_hL_LaRibhLPveP1500,s=1,color='orangered')
# # ax[0].plot(Langmuir_x, Langmuir_y, 'k'); ax[0].plot(Wind_x, Wind_y, 'k'); ax[0].plot(Convection_x,Convection_y, 'k')
# # ax[0].text(0.15, 0.005, 'LANGMUIR', ha='center', va='center')
# # ax[0].text(0.15, 10000, 'CONVECTION', ha='center', va='center')
# # ax[0].text(6, 0.005, 'WIND', ha='center', va='center')
# # ax[0].set_ylabel(r'$h/L_L$ without solar radiation',fontsize=12); ax[0].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
# # ax[0].set_xlabel(r'$La_T$',fontsize=12); ax[0].set_xscale('log')
# # ax[0].set_ylim([0.001,1000000])
# # ax[0].set_xlim([0.1,10])
# # ax[0].tick_params(axis='both', which='major', labelsize=20)
# # ax[0].grid()
# # # ax[0].axis([0.1,10,0.001,1000000])
# # # ax[0].set_yticks([0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000])


# # ax[1].scatter(GTM_Lat_LaRibhLPveP1500,np.negative(GTM_hL_LaRibhLPveP1500),s=1,color='orangered')
# # ax[1].set_ylabel(r'-$h/L_L$ without solar radiation',fontsize=12); ax[1].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
# # ax[1].set_xlabel(r'$La_T$',fontsize=12); ax[1].set_xscale('log')
# # ax[1].set_ylim([0.001,1000000])
# # ax[1].set_xlim([0.1,10])
# # ax[1].tick_params(axis='both', which='major', labelsize=20)
# # ax[1].grid()
# # # ax[0].axis([0.1,10,0.001,1000000])
# # # ax[1].get_yticks([0.01,0,100,10000,1000000])
# # # ax[1].get_yaxis().set_major_formatter(LogFormatter())

# # ax[0].text(0.03, 0.95, '(a)', horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes)
# # ax[1].text(0.03, 0.95, '(b)', horizontalalignment='left', verticalalignment='center', transform=ax[1].transAxes)


# # plt.savefig("regime_diagram.png")


# plt.show()

#################################################################################
## Evolution of u*, u_s and w'theta' using hL thresholds ranging [3000,3500,4000]  ##
#################################################################################

fig, ax = plt.subplots(3,1,sharex=True,figsize=(30,15))

ax[0].plot(time2LaRibhLPveP1500,GTM_fric_LaRibhLPveP1500,color='black',linestyle='solid',label = r'ERA-I')
ax[0].plot(time2LaRibhLPveP1500_2030,GTM_fric_LaRibhLPveP1500_2030,color='orangered',linestyle='solid',label = r'ERA-I averaged by 2 wks')
ax[0].grid(); ax[0].legend(loc='best',fontsize=23)
ax[0].set_ylabel(r'$u_*$ -$ms^{-1}$',fontsize=25);
ax[0].set_xlim(250,650)
ax[0].tick_params(axis='both', which='major', labelsize=20)
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(time2LaRibhLPveP1500,GTM_us_LaRibhLPveP1500,color='black',linestyle='solid')#,label = r'GOTM')
ax[1].plot(time2LaRibhLPveP1500_2030,GTM_us_LaRibhLPveP1500_2030,color='orangered',linestyle='solid')#,label = r'GOTM SMA = 2 wks')
ax[1].grid();# ax[1].legend(loc='best',fontsize='small')
ax[1].set_ylabel(r'$u_s$ -$ms^{-1}$',fontsize=25); ax[1].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
ax[1].tick_params(axis='both', which='major', labelsize=20)
ax[1].set_xlim(265,635)
# ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

ax[2].plot(time2LaRibhLPveP1500,GTM_Theatflux_LaRibhLPveP1500,color='black',linestyle='solid')#,label = r'GOTM')
ax[2].plot(time2LaRibhLPveP1500_2030,GTM_Theatflux_LaRibhLPveP1500_2030,color='orangered',linestyle='solid')#,label = r'GOTM SMA = 2 wks')
ax[2].grid();# ax[2].legend(loc='best',fontsize='small')
ax[2].set_ylabel(r'$R_n + H + \lambda E$ -$Wm^{-2}$',fontsize=25);
ax[2].set_xlabel('Time -Julian days',fontsize=25);
ax[2].tick_params(axis='both', which='major', labelsize=20)

# ax[2].set_ylim(0,1) $w \prime \theta \prime$
# ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

ax0=ax[0].twiny(); 
ax0.set_axisbelow(True)
# ax1.set_xlim(mdate[0],mdate[-1])
xa0 = ax0.get_xaxis(); 
xa0.axis_date(); 
xa0.set_minor_locator(mdates.DayLocator(interval=1)); 
xa0.set_major_locator(mdates.MonthLocator()); 
xa0.set_major_formatter(mdates.DateFormatter('%b-%y'))
ax0.tick_params(axis='both', which='major', labelsize=18)
ax0.plot(yearmonth,fric_u,color='white',linestyle='dotted')


# fig.text(0.5, 0.01, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.02, 0.94, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes,fontsize=20)
ax[1].text(0.02, 0.94, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes,fontsize=20)
ax[2].text(0.02, 0.94, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes,fontsize=20)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
plt.savefig("timeseries_u*_us_heat.png")

plt.show()

#################################################################################
## Evolution of MLDs, hL and Lat using hL thresholds ranging [3000,3500,4000]  ##
#################################################################################

fig, ax = plt.subplots(2,1,sharex=True,figsize=(30,15))

# ax[0].plot(df_dats5,np.negative(df_mld5), color='maroon',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[0].plot(time2LaRibL,mldKPPLaRib_Local,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
# ax[0].plot(time2LaRibhLPveP1350,mldKPPLaRibhLPveP4000_Local,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
# ax[0].set_ylabel('Mixed layer depth -m');
# ax[0].set_xlim(250,650)
# # ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# ax.scatter(time2LaRibhLPveP1500,GTM_hLR_LaRibhLPveP1500,s=2,color='black')
ax[1].plot(time2LaRib,GTM_hLR_LaRib,color='black',linestyle='solid')#,label = r'GOTM$')
ax[1].plot(time2LaRib_2030,GTM_hLR_LaRib_2030,color='orangered',linestyle='solid')#,label = r'GOTM SMA=2 wks')
ax[1].grid();# ax[1].legend(loc='best',fontsize='small')
ax[1].set_ylabel(r'$h/L_L$',fontsize=25); ax[1].set_yscale('symlog',nonposy='clip')#, linthreshy=0.01)
ax[1].set_xlabel(r'Time -Julian days',fontsize=25);
ax[1].set_yticks((-100000000,-1000000,-10000,-100,0,100,10000,1000000,100000000))
ax[1].set_xlim(265,635)
ax[1].tick_params(axis='both', which='major', labelsize=20)

# ax[1].set_yticks((-10000000000,-100000000,-1000000,-10000,-100,0,0.01,1,100,10000,1000000,100000000,10000000000))
# ax[1].set_yticks(-10000000000,10000000000,100)
# ax[1].get_yaxis().set_major_locator(LogLocator(base=2))
# ax[1].set_ylim(-10^10,10^10)

# ax[1].set_xlim(490,635); 

ax0=ax[0].twiny(); 
ax0.set_axisbelow(True)
# ax0.set_xlim(yearmonth[0],yearmonth[-0])
xa0 = ax0.get_xaxis(); 
xa0.axis_date(); 
xa0.set_minor_locator(mdates.DayLocator(interval=1)); 
xa0.set_major_locator(mdates.MonthLocator()); 
xa0.set_major_formatter(mdates.DateFormatter('%b-%y'))
# ax0.plot(datetime2,GTM_Lat_LaRib,color='')
ax0.plot(yearmonth,la_t,color='white')
ax0.tick_params(axis='both', which='major', labelsize=20)


ax[0].plot(time2LaRib,GTM_Lat_LaRib,color='black',linestyle='solid',label = r'ERA-I')
ax[0].plot(time2LaRib_2030,GTM_Lat_LaRib_2030,color='orangered',linestyle='solid',label = r'ERA-I averaged by 2 wks')
ax[0].grid(); ax[0].legend(loc='best',fontsize=23)
ax[0].set_ylabel(r'$La_T$',fontsize=30);
ax[0].set_ylim(0,1)
ax[0].tick_params(axis='both', which='major', labelsize=20)


# fig.text(0.5, 0.01, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes,fontsize=20)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes,fontsize=20)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
plt.savefig("timeseries_hL_LaT.png")
plt.show()


exit()



exit()

###################################################################################################################
## Compare evolution of MLDs with NON-LOCAL MIXING vs LOCAL MIXING using hL thresholds ranging [3000,3500,4000]  ## 
###################################################################################################################

fig, ax = plt.subplots(2,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='maroon',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibhLPveP3000,mldKPPLaRibhLPveP3000,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_ylabel('MLD (NON-LOCAL enabled) -m');
ax[0].set_xlim(265,650)
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='maroon',linestyle='solid',label = 'OSMOSIS SMA n=5')
ax[1].plot(time2L,mldKPP_Local,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRibL,mldKPPLaRib_Local,color='steelblue',linestyle='solid',label = r'KPP $La_T$ $Ri_B$ $\gamma=\alpha=1$')
ax[1].plot(time2LaRibhLPveP1500,mldKPPLaRibhLPveP1500,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_ylabel('MLD (LOCAL) -m');
ax[1].set_xlim(265,650)
# ax[1].set_xlim(480,650)
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
plt.show()

exit()
# ########################################################################################
# ## Evolution of tdotus, u* and heatflux using hL thresholds ranging [3000,3500,4000]  ##
# ########################################################################################

# fig, ax = plt.subplots(3,1)

# ax[0].plot(time2LaRibhLPveP1350,GTM_fric_LaRibhLPveP1350,color='black',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[0].plot(time2LaRibhLPveP1350_2030,GTM_fric_LaRibhLPveP1350_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[0].grid();# ax[0].legend(loc='best',fontsize='small')
# ax[0].set_ylabel(r'$u_{*}^2$ - $ms^{-1}$');
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# ax[1].plot(time2LaRibhLPveP1350,GTM_Theatflux_LaRibhLPveP1350,color='black',linestyle='solid',label = r'KPP original  $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[1].plot(time2LaRibhLPveP1350_2030,GTM_Theatflux_LaRibhLPveP1350_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[1].grid();# ax[1].legend(loc='best',fontsize='small')
# ax[1].set_ylabel(r'$[H + \lambda E + R]$ - $Wm^{-2}$');
# ax[1].set_xlim(490, 635); ax[1].set_ylim(-150,0)

# ax[2].plot(time2LaRibhLPveP1350,GTM_us_LaRibhLPveP1350,color='black',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[2].plot(time2LaRibhLPveP1350_2030,GTM_us_LaRibhLPveP1350_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
# ax[2].set_ylabel(r'$u_s$ - $ms^{-1}$');
# ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

# fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# # fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
# ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# # plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
# plt.show()

# exit()
# ##############################################################################
# ## Comparison of hL time series w & wout SMA ranging [1500,3000,3500,4000]  ##
# ##############################################################################


# fig, ax = plt.subplots(4,1)

# ax[0].plot(time2LaRibhLPveP1500,GTM_hLR_LaRibhLPveP2000,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>2000$')
# ax[0].plot(time2LaRibhLPveP2000_2030,GTM_hLR_LaRibhLPveP2000_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>2000$')
# ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# ax[1].plot(time2LaRibhLPveP3000,GTM_hLR_LaRibhLPveP3000,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
# ax[1].plot(time2LaRibhLPveP3000_2030,GTM_hLR_LaRibhLPveP3000_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
# ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
# ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

# ax[2].plot(time2LaRibhLPveP1600,GTM_hLR_LaRibhLPveP1600,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
# ax[2].plot(time2LaRibhLPveP1600_2030,GTM_hLR_LaRibhLPveP1600_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
# ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
# ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

# ax[3].plot(time2LaRibhLPveP4000,GTM_hLR_LaRibhLPveP4000,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[3].plot(time2LaRibhLPveP4000_2030,GTM_hLR_LaRibhLPveP4000_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[3].grid(); ax[3].legend(loc='best',fontsize='small')
# ax[3].set_xlim(490,635); ax[3].set_ylim(-150,0)

# fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
# ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
# ax[3].text(0.06, 0.08, '(d)', horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of hL time series for simualtions with & without SMA ranging [2000,3000,3500,4000]')
# plt.show()

# exit()
#######################################################################################
## Comparison of MLDs with multiplied factor on original diffusivities gamma=alpha=0 ##
#######################################################################################


fig, ax = plt.subplots(3,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibnum15,mldKPPLaRibnum15,color='darkorange',linestyle='dashed',label = r'KPP $\nu_{m}*1.5 $')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(480,635); 
# ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibnuh15,mldKPPLaRibnuh15,color='darkorange',linestyle='dashed',label = r'KPP $\nu_{h}*1.5 $')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(480,635); 
# ax[1].set_ylim(-150,0)

ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[2].plot(time2LaRibnum15nuh15,mldKPPLaRibnum15nuh15,color='darkorange',linestyle='dashed',label = r'KPP $\nu_{h}*1.5 \nu_{m}*1.5$')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_xlim(480,635); 
# ax[2].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs with multiplied factor on original diffusivities')
plt.show()


################################################################
## Comparison of MLDs with different values for gamma & alpha ##
################################################################


fig, ax = plt.subplots(2,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibg0a1,mldKPPLaRibg0a1,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w g=0 a=1')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(480,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibg1a0,mldKPPLaRibg1a0,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w g=1 a=0')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(480,635); ax[1].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs with different values for gamma & alpha')
plt.show()
exit()
########################################################################
## Comparison of MLDs with threshold condition hL>A on viscosity ##
########################################################################


fig, ax = plt.subplots(3,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibhLgamPveN100,mldKPPLaRibhLgamPveN100,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w gam hL > -100')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibhLgamPve,mldKPPLaRibhLgamPve,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w gam hL > 0')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[2].plot(time2LaRibhLgamPveP100,mldKPPLaRibhLgamPveP100,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w gam hL > 100')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs with threshold conditions of hL on viscosity')
plt.show()


################################################################################
## Comparison of MLDs with threshold condition hL>A on original difusivities ##
################################################################################


fig, ax = plt.subplots(3,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2hLPveN100,mldKPPhLPveN100,color='darkorange',linestyle='dashed',label = 'KPP w hL > -100')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2hLPve,mldKPPhLPve,color='darkorange',linestyle='dashed',label = 'KPP w hL > 0')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[2].plot(time2hLPveP100,mldKPPhLPveP100,color='darkorange',linestyle='dashed',label = 'KPP w hL > 100')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs with threshold conditions of hL on viscosity')
plt.show()

exit()
########################################################################
## Time series of Mixed Layer Depths and hL with condition hL>-A ##
########################################################################


fig, ax = plt.subplots(3,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibLaTPve03,mldKPPLaRibLaTPve03,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w LaT>0.3')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_ylabel('Mixed Layer Depth -m')
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)
ax[0].set_xlim(250,635);

ax[1].plot(time2LaRibLaTPve03_2030,GTM_hLR_LaRibLaTPve03_2030,color='darkorange',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1 w LaT>0.3')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_ylabel('h/L SMA_n=2wks')
ax[1].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
# ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)
ax[1].set_xlim(250,635);

ax[2].plot(time2LaRibLaTPve03,GTM_Lat_LaRibLaTPve03,color='darkorange',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1 w LaT>0.3')
ax[2].grid(); ax[1].legend(loc='best',fontsize='small')
ax[2].set_ylabel('La_T')
# ax[2].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
# ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)
ax[2].set_xlim(250,635); ax[2].set_ylim(0,1)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Time series of MLDs and hL with condition hL>-A')
plt.show()

exit()

########################################################################
## Comparison of MLDs with and without shear instability with hL>-100 ##
########################################################################


fig, ax = plt.subplots(2,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibhLPveN100,mldKPPLaRibhLPveN100,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > -100')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibhLPveSSN100,mldKPPLaRibhLPveSSN100,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > -100 w SS')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs with and without shear instability with hL>-100')
plt.show()

exit()
###########################################################################
## Comparison of MLDs varying the constant (Cm, Ch) values with hL>-1000 ##
###########################################################################

fig, ax = plt.subplots(4,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibhLPveN1000,mldKPPLaRibhLPveN1000,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -1000')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibhLa1g2PveN1000,mldKPPLaRibhLa1g2PveN1000,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=2, Ch=1 w hL cr < -1000')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[2].plot(time2LaRibhLa2g1PveN1000,mldKPPLaRibhLa2g1PveN1000,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=1, Ch=2 w hL cr < -1000')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

ax[3].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[3].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[3].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[3].plot(time2LaRibhLPveN1500,mldKPPLaRibhLPveN1500,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -1500')
ax[3].grid(); ax[3].legend(loc='best',fontsize='small')
ax[3].set_xlim(490,635); ax[3].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs varying the constant (Cm, Ch) values with hL>-1000')
plt.show()

###########################################################################
## Comparison of MLDs varying the constant (Cm, Ch) values with hL>-500 ##
###########################################################################

fig, ax = plt.subplots(3,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibhLPveN500,mldKPPLaRibhLPveN500,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -500')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibhLa1g2PveN500,mldKPPLaRibhLa1g2PveN500,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=2, Ch=1 w hL cr < -500')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
ax[2].plot(time2LaRibhLa2g1PveN500,mldKPPLaRibhLa2g1PveN500,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=1, Ch=2 w hL cr < -500')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of MLDs varying the constant (Cm, Ch) values with hL>-500')
