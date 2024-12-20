## analysis of salinity ##
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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, ScalarFormatter)

# plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

### ----------------------------- GLIDER ----------------------------- ##
gliderdata= xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc'); 
glider = xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days
glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})

### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
## LaT addition ##
gotmkppLaRib = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib.nc",decode_times=False)
gotmkppnum2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_num2.nc",decode_times=False)
gotmkppnuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nuh2.nc",decode_times=False)
# gotmkppRibnum2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_LaT_Rib_num2.nc",decode_times=False)
gotmkppRib05num2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_Rib05_num2.nc",decode_times=False)
# gotmkppRibnuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_LaT_Rib_nuh2.nc",decode_times=False)
gotmkppRib2nuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_Rib2_nuh2.nc",decode_times=False)
## h/L correction ##
gotmkppLaRibhLPveP100 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP100.nc",decode_times=False)
# gotmkppLaRibhLPveP200 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP200.nc",decode_times=False)
# gotmkppLaRibhLPveP1000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1000.nc",decode_times=False)
# gotmkppLaRibhLPveP1300 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1300.nc",decode_times=False)
# gotmkppLaRibhLPveP1350 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1350.nc",decode_times=False)
# gotmkppLaRibhLPveP1400 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1400.nc",decode_times=False)
gotmkppLaRibhLPveP1500 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1500.nc",decode_times=False)
# gotmkppLaRibhLPveP1600 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1600.nc",decode_times=False)
# gotmkppLaRibhLPveP1750 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1750.nc",decode_times=False)
# gotmkppLaRibhLPveP2000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP2000.nc",decode_times=False)
# gotmkppLaRibhLPveP4000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP4000.nc",decode_times=False)
# gotmkppLaRibhLPveP5000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP5000.nc",decode_times=False)
gotmkppLaRibhLPveP10000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP10000.nc",decode_times=False)

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
GTMtempnum2 = gotmkppnum2.temp; GTMdepthnum2 = gotmkppnum2.z; GTMnumnum2 = gotmkppnum2.num; GTMnuhnum2 = gotmkppnum2.nuh; GTMunum2 = gotmkppnum2.u; GTMvnum2 = gotmkppnum2.v; timenum2 = gotmkppnum2.time
GTMtempnuh2 = gotmkppnuh2.temp; GTMdepthnuh2 = gotmkppnuh2.z; GTMnumnuh2 = gotmkppnuh2.num; GTMnuhnuh2 = gotmkppnuh2.nuh; GTMunuh2 = gotmkppnuh2.u; GTMvnuh2 = gotmkppnuh2.v; timenuh2 = gotmkppnuh2.time
# GTMtempRibnum2 = gotmkppRibnum2.temp; GTMdepthRibnum2 = gotmkppRibnum2.z; GTMnumRibnum2 = gotmkppRibnum2.num; GTMnuhRibnum2 = gotmkppRibnum2.nuh; GTMuRibnum2 = gotmkppRibnum2.u; GTMvRibnum2 = gotmkppRibnum2.v; timeRibnum2 = gotmkppRibnum2.time
GTMtempRib05num2 = gotmkppRib05num2.temp; GTMdepthRib05num2 = gotmkppRib05num2.z; GTMnumRib05num2 = gotmkppRib05num2.num; GTMnuhRib05num2 = gotmkppRib05num2.nuh; GTMuRib05num2 = gotmkppRib05num2.u; GTMvRib05num2 = gotmkppRib05num2.v; timeRib05num2 = gotmkppRib05num2.time
# GTMtempRibnuh2 = gotmkppRibnuh2.temp; GTMdepthRibnuh2 = gotmkppRibnuh2.z; GTMnumRibnuh2 = gotmkppRibnuh2.num; GTMnuhRibnuh2 = gotmkppRibnuh2.nuh; GTMuRibnuh2 = gotmkppRibnuh2.u; GTMvRibnuh2 = gotmkppRibnuh2.v; timeRibnuh2 = gotmkppRibnuh2.time
GTMtempRib2nuh2 = gotmkppRib2nuh2.temp; GTMdepthRib2nuh2 = gotmkppRib2nuh2.z; GTMnumRib2nuh2 = gotmkppRib2nuh2.num; GTMnuhRib2nuh2 = gotmkppRib2nuh2.nuh; GTMuRib2nuh2 = gotmkppRib2nuh2.u; GTMvRib2nuh2 = gotmkppRib2nuh2.v; timeRib2nuh2 = gotmkppRib2nuh2.time

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
time2num2 = 267.75 + (timenum2[:]/86400); 
time2nuh2 = 267.75 + (timenuh2[:]/86400); 
# time2Ribnum2 = 267.75 + (timeRibnum2[:]/86400); 
time2Rib05num2 = 267.75 + (timeRib05num2[:]/86400); 
# time2Ribnuh2 = 267.75 + (timeRibnuh2[:]/86400); 
time2Rib2nuh2 = 267.75 + (timeRib2nuh2[:]/86400); 

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

###################
# time allocation #
###################

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
# gotmkppL = gotmkppL.assign(time2L=("time", time2L)); gotmkppL = gotmkppL.swap_dims({"time" : "time2L"})
# gotmkppLa = gotmkppLa.assign(time2La=("time", time2La)); gotmkppLa = gotmkppLa.swap_dims({"time" : "time2La"})
gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})
# gotmkppLaRibL = gotmkppLaRibL.assign(time2LaRibL=("time", time2LaRibL)); gotmkppLaRibL = gotmkppLaRibL.swap_dims({"time" : "time2LaRibL"})
gotmkppnum2 = gotmkppnum2.assign(time2num2=("time", time2num2)); gotmkppnum2 = gotmkppnum2.swap_dims({"time" : "time2num2"})
gotmkppnuh2 = gotmkppnuh2.assign(time2nuh2=("time", time2nuh2)); gotmkppnuh2 = gotmkppnuh2.swap_dims({"time" : "time2nuh2"})
# gotmkppRibnum2 = gotmkppRibnum2.assign(time2Ribnum2=("time", time2Ribnum2)); gotmkppRibnum2 = gotmkppRibnum2.swap_dims({"time" : "time2Ribnum2"})
gotmkppRib05num2 = gotmkppRib05num2.assign(time2Rib05num2=("time", time2Rib05num2)); gotmkppRib05num2 = gotmkppRib05num2.swap_dims({"time" : "time2Rib05num2"})
# gotmkppRibnuh2 = gotmkppRibnuh2.assign(time2Ribnuh2=("time", time2Ribnuh2)); gotmkppRibnuh2 = gotmkppRibnuh2.swap_dims({"time" : "time2Ribnuh2"})
gotmkppRib2nuh2 = gotmkppRib2nuh2.assign(time2Rib2nuh2=("time", time2Rib2nuh2)); gotmkppRib2nuh2 = gotmkppRib2nuh2.swap_dims({"time" : "time2Rib2nuh2"})

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


######################
## -- timeseries -- ## 
######################

# select surface 
GldrSalt = GldrSalt.isel(pressure=0)
print('Sea Surface Salinity', GldrSalt)

plt.plot(new_time,GldrSalt,color='crimson',linestyle='solid',label='OSMOSIS')
plt.xlabel('Time -Julian days'); plt.ylabel('SSS')
plt.xlim([240,635]); plt.legend(loc='best',fontsize='small')
plt.grid();
plt.show();