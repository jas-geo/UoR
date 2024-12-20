## ------------- Written code to examine the periods with low La_t ------------- ##
import numpy as np
import math
import os
import csv
import sys
import xarray as xr
import pandas as pd
from scipy import interpolate
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
from operator import itemgetter
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  


### ----------------------------- GLIDER ----------------------------- ##

gliderdata= Dataset('../glider_timeseries.nc'); glider = xr.open_dataset("../glider_timeseries.nc", decode_times=False)

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
# print('Glider time until 2012-09-22',new_time[144:202].values)
# print('Glider time after 2012-10-03',new_time[316:382].values)
# print(gliderdata.variables.keys())
GldrSalt = glider.prac_sal;# print(GldrSalt[309,:])

# print('glider time = 277', new_time[309])
# gliderTemp277 = np.c_[GldrDepth[0:200],GldrTemp[309,0:200]], 
# print('glider temp values at time = 277',gliderTemp277);
# gliderSalt277 = np.c_[GldrDepth[0:200],GldrSalt[309,0:200]], 
# print('glider temp values at time = 277',gliderSalt277)

gliderTemp_Depth0 = np.c_[np.arange(308,383),GldrTemp[308:383,0]]
print('Glider depth', GldrDepth[0].values)
print('Glider temp for depth=-1', GldrTemp[308:383,0].values)
print('Glider temp for depth=-1', gliderTemp_Depth0)
print(np.nanmin(GldrTemp[308:383,0].values),np.nanmax(GldrTemp[308:383,0].values),np.nanmedian(GldrTemp[308:383,0].values),len(GldrTemp[308:383,0].values))

gliderTemp_Depth0 = gliderTemp_Depth0[np.argsort(gliderTemp_Depth0[:, 1])]
print('glider Temp at depth = -1 m', gliderTemp_Depth0)

# arr = np.array([356,349,358,365,364,314,340,311,378,369,319,329,327,316,382])
# arr = np.array([356,358,364,340,378,319,327,382])
# arr = arr[np.argsort(arr)]
# print(arr)
### ----------------------------- GOTM ----------------------------- ###

# with (PAP-SO)
gotmkpp = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_comparison_PAP.nc", decode_times=False)
gotmkppLoc = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_comparisonLoc_PAP.nc", decode_times=False)
gotmkeps = xr.open_dataset("autumn_highLatkeps/OSMOSIS_cruise_winter.nc", decode_times=False)
# changed the definition of eddy viscosity (num(k)) in kpp.F90
gotmkppLa = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_PAP.nc", decode_times=False)
# changed the definition of eddy viscosity (num(k)) in kpp.F90 wihtout non-local fluxes
gotmkppLaLoc = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyViscLoc_PAP.nc", decode_times=False)
gotmkppLag05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_0.5_PAP.nc", decode_times=False)
gotmkppLag2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_2_PAP.nc", decode_times=False)
gotmkppLa05dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_0.5dlta_PAP.nc", decode_times=False)
gotmkppLa025dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_0.25dlta_PAP.nc", decode_times=False)
gotmkppLa125dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_0.125dlta_PAP.nc", decode_times=False)
gotmkppLa0625dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newEddyVisc_0.0625dlta_PAP.nc", decode_times=False)
# changed the definition of heat diffusivity (nuh(k)) in kpp.F90
gotmkpp_heat = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_PAP.nc", decode_times=False)
# changed the definition of heat diffusivity (nuh(k)) in kpp.F90 without non-local fluxes
gotmkpp_heatLoc = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiffLoc_PAP.nc", decode_times=False)
gotmkpp_heatg05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_0.5_PAP.nc", decode_times=False)
gotmkpp_heatg2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_2_PAP.nc", decode_times=False)
gotmkpp_heat05dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_0.5dlta_PAP.nc", decode_times=False)
gotmkpp_heat025dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_0.25dlta_PAP.nc", decode_times=False)
gotmkpp_heat125dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_0.125dlta_PAP.nc", decode_times=False)
gotmkpp_heat0625dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newHeatDiff_0.0625dlta_PAP.nc", decode_times=False)
# changed the definition of heat diffusitivity (nuh(k)) in kpp.F90 ##
gotmkpp_m_h = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_new_ViscDiff_PAP.nc", decode_times=False)
gotmkpp_m_hLoc = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_new_ViscDiffLoc_PAP.nc", decode_times=False)
gotmkpp_m_hg05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newViscDiff_0.5_PAP.nc", decode_times=False)
gotmkpp_m_hg2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newViscDiff_2_PAP.nc", decode_times=False)
gotmkpp_m_h05dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newViscDiff_0.5dlta_PAP.nc", decode_times=False)
gotmkpp_m_h025dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newViscDiff_0.25dlta_PAP.nc", decode_times=False)
gotmkpp_m_h125dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newViscDiff_0.125dlta_PAP.nc", decode_times=False)
gotmkpp_m_h0625dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/autumn_highLatKPP/OSMOSIS_newViscDiff_0.0625dlta_PAP.nc", decode_times=False)

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMdepth1 = gotmkpp.z1; GTMu = gotmkpp.u; GTMv = gotmkpp.v; time = gotmkpp.time
GTMtempLoc = gotmkppLoc.temp; GTMdepthLoc = gotmkppLoc.z; timeLoc = gotmkppLoc.time
GTMtempLa = gotmkppLa.temp; GTMdepthLa = gotmkppLa.z; GTMnumLa = gotmkppLa.num; GTMnuhLa = gotmkppLa.nuh; GTMdepth1La = gotmkppLa.z1; GTMuLa = gotmkppLa.u; GTMvLa = gotmkppLa.v;  timeLa = gotmkppLa.time
GTMtempLaLoc = gotmkppLaLoc.temp; GTMdepthLaLoc = gotmkppLaLoc.z; GTMnumLaLoc = gotmkppLaLoc.num; GTMnuhLaLoc = gotmkppLaLoc.nuh; GTMdepth1LaLoc = gotmkppLaLoc.z1; GTMuLaLoc = gotmkppLaLoc.u; GTMvLaLoc = gotmkppLaLoc.v;  timeLaLoc = gotmkppLaLoc.time
GTMtempLag05 = gotmkppLag05.temp; GTMdepthLag05 = gotmkppLag05.z; GTMnumLag05 = gotmkppLag05.num; GTMnuhLag05 = gotmkppLag05.nuh; GTMdepth1Lag05 = gotmkppLag05.z1; GTMuLag05 = gotmkppLag05.u; GTMvLag05 = gotmkppLag05.v;  timeLag05 = gotmkppLag05.time
GTMtempLag2 = gotmkppLag2.temp; GTMdepthLag2 = gotmkppLag2.z; GTMnumLag2 = gotmkppLag2.num; GTMnuhLag2 = gotmkppLag2.nuh; GTMdepth1Lag2 = gotmkppLag2.z1; GTMuLag2 = gotmkppLag2.u; GTMvLag2 = gotmkppLag2.v;  timeLag2 = gotmkppLag2.time
GTMtempLa05dlta = gotmkppLa05dlta.temp; GTMdepthLa05dlta = gotmkppLa05dlta.z; GTMnumLa05dlta = gotmkppLa05dlta.num; GTMnuhLa05dlta = gotmkppLa05dlta.nuh; GTMdepth1La05dlta = gotmkppLa05dlta.z1; GTMuLa05dlta = gotmkppLa05dlta.u; GTMvLa05dlta = gotmkppLa05dlta.v;  timeLa05dlta = gotmkppLa05dlta.time
GTMtempLa025dlta = gotmkppLa025dlta.temp; GTMdepthLa025dlta = gotmkppLa025dlta.z; GTMnumLa025dlta = gotmkppLa025dlta.num; GTMnuhLa025dlta = gotmkppLa025dlta.nuh; GTMdepth1La025dlta = gotmkppLa025dlta.z1; GTMuLa025dlta = gotmkppLa025dlta.u; GTMvLa025dlta = gotmkppLa025dlta.v;  timeLa025dlta = gotmkppLa025dlta.time
GTMtempLa125dlta = gotmkppLa125dlta.temp; GTMdepthLa125dlta = gotmkppLa125dlta.z; GTMnumLa125dlta = gotmkppLa125dlta.num; GTMnuhLa125dlta = gotmkppLa125dlta.nuh; GTMdepth1La125dlta = gotmkppLa125dlta.z1; GTMuLa125dlta = gotmkppLa125dlta.u; GTMvLa125dlta = gotmkppLa125dlta.v;  timeLa125dlta = gotmkppLa125dlta.time
GTMtempLa0625dlta = gotmkppLa0625dlta.temp; GTMdepthLa0625dlta = gotmkppLa0625dlta.z; GTMnumLa0625dlta = gotmkppLa0625dlta.num; GTMnuhLa0625dlta = gotmkppLa0625dlta.nuh; GTMdepth1La0625dlta = gotmkppLa0625dlta.z1; GTMuLa0625dlta = gotmkppLa0625dlta.u; GTMvLa0625dlta = gotmkppLa0625dlta.v;  timeLa0625dlta = gotmkppLa0625dlta.time
GTMtemp_heat = gotmkpp_heat.temp; GTMdepth_heat = gotmkpp_heat.z; GTMnum_heat = gotmkpp_heat.num; GTMdepth1_heat = gotmkpp_heat.z1; GTMnuh_heat = gotmkpp_heat.nuh; GTMu_heat = gotmkpp_heat.u; GTMv_heat = gotmkpp_heat.v; time_heat = gotmkpp_heat.time
GTMtemp_heatLoc = gotmkpp_heatLoc.temp; GTMdepth_heatLoc = gotmkpp_heatLoc.z; GTMnum_heatLoc = gotmkpp_heatLoc.num; GTMdepth1_heatLoc = gotmkpp_heatLoc.z1; GTMnuh_heatLoc = gotmkpp_heatLoc.nuh; GTMu_heatLoc = gotmkpp_heatLoc.u; GTMv_heatLoc = gotmkpp_heatLoc.v; time_heatLoc = gotmkpp_heatLoc.time
GTMtemp_heatg05 = gotmkpp_heatg05.temp; GTMdepth_heatg05 = gotmkpp_heatg05.z; GTMnum_heatg05 = gotmkpp_heatg05.num; GTMnuh_heatg05 = gotmkpp_heatg05.nuh; GTMdepth1_heatg05 = gotmkpp_heatg05.z1; GTMu_heatg05 = gotmkpp_heatg05.u; GTMv_heatg05 = gotmkpp_heatg05.v; time_heatg05 = gotmkpp_heatg05.time
GTMtemp_heatg2 = gotmkpp_heatg2.temp; GTMdepth_heatg2 = gotmkpp_heatg2.z; GTMnum_heatg2 = gotmkpp_heatg2.num; GTMnuh_heatg2 = gotmkpp_heatg2.nuh; GTMdepth1_heatg2 = gotmkpp_heatg2.z1; GTMu_heatg2 = gotmkpp_heatg2.u; GTMv_heatg2 = gotmkpp_heatg2.v;  time_heatg2 = gotmkpp_heatg2.time
GTMtemp_heat05dlta = gotmkpp_heat05dlta.temp; GTMdepth_heat05dlta = gotmkpp_heat05dlta.z; GTMnum_heat05dlta = gotmkpp_heat05dlta.num; GTMnuh_heat05dlta = gotmkpp_heat05dlta.nuh; GTMdepth1_heat05dlta = gotmkpp_heat05dlta.z1; GTMu_heat05dlta = gotmkpp_heat05dlta.u; GTMv_heat05dlta = gotmkpp_heat05dlta.v;  time_heat05dlta = gotmkpp_heat05dlta.time
GTMtemp_heat025dlta = gotmkpp_heat025dlta.temp; GTMdepth_heat025dlta = gotmkpp_heat025dlta.z; GTMnum_heat025dlta = gotmkpp_heat025dlta.num; GTMnuh_heat025dlta = gotmkpp_heat025dlta.nuh; GTMdepth1_heat025dlta = gotmkpp_heat025dlta.z1; GTMu_heat025dlta = gotmkpp_heat025dlta.u; GTMv_heat025dlta = gotmkpp_heat025dlta.v;  time_heat025dlta = gotmkpp_heat025dlta.time
GTMtemp_heat125dlta = gotmkpp_heat125dlta.temp; GTMdepth_heat125dlta = gotmkpp_heat125dlta.z; GTMnum_heat125dlta = gotmkpp_heat125dlta.num; GTMnuh_heat125dlta = gotmkpp_heat125dlta.nuh; GTMdepth1_heat125dlta = gotmkpp_heat125dlta.z1; GTMu_heat125dlta = gotmkpp_heat125dlta.u; GTMv_heat125dlta = gotmkpp_heat125dlta.v;  time_heat125dlta = gotmkpp_heat125dlta.time
GTMtemp_heat0625dlta = gotmkpp_heat0625dlta.temp; GTMdepth_heat0625dlta = gotmkpp_heat0625dlta.z; GTMnum_heat0625dlta = gotmkpp_heat0625dlta.num; GTMnuh_heat0625dlta = gotmkpp_heat0625dlta.nuh; GTMdepth1_heat0625dlta = gotmkpp_heat0625dlta.z1; GTMu_heat0625dlta = gotmkpp_heat0625dlta.u; GTMv_heat0625dlta = gotmkpp_heat0625dlta.v;  time_heat0625dlta = gotmkpp_heat0625dlta.time
GTMtemp_mh = gotmkpp_m_h.temp; GTMdepth_mh = gotmkpp_m_h.z; GTMnum_mh = gotmkpp_m_h.num; GTMnuh_mh = gotmkpp_m_h.nuh; GTMdepth1_mh = gotmkpp_m_h.z1; GTMu_mh = gotmkpp_m_h.u; GTMv_mh = gotmkpp_m_h.v; time_mh = gotmkpp_m_h.time
GTMtemp_mhLoc = gotmkpp_m_hLoc.temp; GTMdepth_mhLoc = gotmkpp_m_hLoc.z; GTMnum_mhLoc = gotmkpp_m_hLoc.num; GTMnuh_mhLoc = gotmkpp_m_hLoc.nuh; GTMdepth1_mhLoc = gotmkpp_m_hLoc.z1; GTMu_mhLoc = gotmkpp_m_hLoc.u; GTMv_mhLoc = gotmkpp_m_hLoc.v; time_mhLoc = gotmkpp_m_hLoc.time
GTMtemp_mhg05 = gotmkpp_m_hg05.temp; GTMdepth_mhg05 = gotmkpp_m_hg05.z; GTMnum_mhg05 = gotmkpp_m_hg05.num; GTMnuh_mhg05 = gotmkpp_m_hg05.nuh; GTMdepth1_mhg05 = gotmkpp_m_hg05.z1; GTMu_mhg05 = gotmkpp_m_hg05.u; GTMv_mhg05 = gotmkpp_m_hg05.v; time_mhg05 = gotmkpp_m_hg05.time
GTMtemp_mhg2 = gotmkpp_m_hg2.temp; GTMdepth_mhg2 = gotmkpp_m_hg2.z; GTMnum_mhg2 = gotmkpp_m_hg2.num; GTMnuh_mhg2 = gotmkpp_m_hg2.nuh; GTMdepth1_mhg2 = gotmkpp_m_hg2.z1; GTMu_mhg2 = gotmkpp_m_hg2.u; GTMv_mhg2 = gotmkpp_m_hg2.v;  time_mhg2 = gotmkpp_m_hg2.time
GTMtemp_mh05dlta = gotmkpp_m_h05dlta.temp; GTMdepth_mh05dlta = gotmkpp_m_h05dlta.z; GTMnum_mh05dlta = gotmkpp_m_h05dlta.num; GTMnuh_mh05dlta = gotmkpp_m_h05dlta.nuh; GTMdepth1_mh05dlta = gotmkpp_m_h05dlta.z1; GTMu_mh05dlta = gotmkpp_m_h05dlta.u; GTMv_mh05dlta = gotmkpp_m_h05dlta.v;  time_mh05dlta = gotmkpp_m_h05dlta.time
GTMtemp_mh025dlta = gotmkpp_m_h025dlta.temp; GTMdepth_mh025dlta = gotmkpp_m_h025dlta.z; GTMnum_mh025dlta = gotmkpp_m_h025dlta.num; GTMnuh_mh025dlta = gotmkpp_m_h025dlta.nuh; GTMdepth1_mh025dlta = gotmkpp_m_h025dlta.z1; GTMu_mh025dlta = gotmkpp_m_h025dlta.u; GTMv_mh025dlta = gotmkpp_m_h025dlta.v;  time_mh025dlta = gotmkpp_m_h025dlta.time
GTMtemp_mh125dlta = gotmkpp_m_h125dlta.temp; GTMdepth_mh125dlta = gotmkpp_m_h125dlta.z; GTMnum_mh125dlta = gotmkpp_m_h125dlta.num; GTMnuh_mh125dlta = gotmkpp_m_h125dlta.nuh; GTMdepth1_mh125dlta = gotmkpp_m_h125dlta.z1; GTMu_mh125dlta = gotmkpp_m_h125dlta.u; GTMv_mh125dlta = gotmkpp_m_h125dlta.v;  time_mh125dlta = gotmkpp_m_h125dlta.time
GTMtemp_mh0625dlta = gotmkpp_m_h0625dlta.temp; GTMdepth_mh0625dlta = gotmkpp_m_h0625dlta.z; GTMnum_mh0625dlta = gotmkpp_m_h0625dlta.num; GTMnuh_mh0625dlta = gotmkpp_m_h0625dlta.nuh; GTMdepth1_mh0625dlta = gotmkpp_m_h0625dlta.z1; GTMu_mh0625dlta = gotmkpp_m_h0625dlta.u; GTMv_mh0625dlta = gotmkpp_m_h0625dlta.v;  time_mh0625dlta = gotmkpp_m_h0625dlta.time
GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkedepth1 = gotmkeps.z1; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time

## Convert time
time2 = 277.0 + (time[:]/86400); ketime2 = 277.0 + (ketime[:]/86400)
time2Loc = 277.0 + (timeLoc[:]/86400); 
time2La = 277.0 + (timeLa[:]/86400); time2LaLoc = 277.0 + (timeLaLoc[:]/86400); time2Lag05 = 277.0 + (timeLag05[:]/86400);  time2Lag2 = 277.0 + (timeLag2[:]/86400);  time2La05dlta = 277.0 + (timeLa05dlta[:]/86400);  time2La025dlta = 277.0 + (timeLa025dlta[:]/86400); time2La125dlta = 277.0 + (timeLa125dlta[:]/86400); time2La0625dlta = 277.0 + (timeLa0625dlta[:]/86400); 
time2_heat = 277.0 + (time_heat[:]/86400); time2_heatLoc = 277.0 + (time_heatLoc[:]/86400); time2_heatg05 = 277.0 + (time_heatg05[:]/86400); time2_heatg2 = 277.0 + (time_heatg2[:]/86400); time2_heat05dlta = 277.0 + (time_heat05dlta[:]/86400); time2_heat025dlta = 277.0 + (time_heat025dlta[:]/86400); time2_heat125dlta = 277.0 + (time_heat125dlta[:]/86400); time2_heat0625dlta = 277.0 + (time_heat0625dlta[:]/86400);
time2_mh = 277.0 + (time_mh[:]/86400); time2_mhLoc = 277.0 + (time_mhLoc[:]/86400);time2_mhg05 = 277.0 + (time_mhg05[:]/86400); time2_mhg2 = 277.0 + (time_mhg2[:]/86400); time2_mh05dlta = 277.0 + (time_heat05dlta[:]/86400); time2_mh025dlta = 277.0 + (time_mh025dlta[:]/86400); time2_mh125dlta = 277.0 + (time_mh125dlta[:]/86400); time2_mh0625dlta = 277.0 + (time_mh0625dlta[:]/86400);

# times = np.c_[time,time2]
# print('GOTM times for high La_t :', time2.values, time2.shape);
gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkppLoc = gotmkppLoc.assign(time2Loc=("time", time2Loc)); gotmkppLoc = gotmkppLoc.swap_dims({"time" : "time2Loc"})
gotmkppLa = gotmkppLa.assign(time2La=("time", time2La)); gotmkppLa = gotmkppLa.swap_dims({"time" : "time2La"})
gotmkppLaLoc = gotmkppLaLoc.assign(time2LaLoc=("time", time2LaLoc)); gotmkppLaLoc = gotmkppLaLoc.swap_dims({"time" : "time2LaLoc"})
gotmkppLag05 = gotmkppLag05.assign(time2Lag05=("time", time2Lag05)); gotmkppLag05 = gotmkppLag05.swap_dims({"time" : "time2Lag05"})
gotmkppLag2 = gotmkppLag2.assign(time2Lag2=("time", time2Lag2)); gotmkppLag2 = gotmkppLag2.swap_dims({"time" : "time2Lag2"})
gotmkppLa05dlta = gotmkppLa05dlta.assign(time2La05dlta=("time", time2La05dlta)); gotmkppLa05dlta = gotmkppLa05dlta.swap_dims({"time" : "time2La05dlta"})
gotmkppLa025dlta = gotmkppLa025dlta.assign(time2La025dlta=("time", time2La025dlta)); gotmkppLa025dlta = gotmkppLa025dlta.swap_dims({"time" : "time2La025dlta"})
gotmkppLa125dlta = gotmkppLa125dlta.assign(time2La125dlta=("time", time2La125dlta)); gotmkppLa125dlta = gotmkppLa125dlta.swap_dims({"time" : "time2La125dlta"})
gotmkppLa0625dlta = gotmkppLa0625dlta.assign(time2La0625dlta=("time", time2La0625dlta)); gotmkppLa0625dlta = gotmkppLa0625dlta.swap_dims({"time" : "time2La0625dlta"})
gotmkpp_heat = gotmkpp_heat.assign(time2_heat=("time", time2_heat)); gotmkpp_heat = gotmkpp_heat.swap_dims({"time" : "time2_heat"})
gotmkpp_heatLoc = gotmkpp_heatLoc.assign(time2_heatLoc=("time", time2_heatLoc)); gotmkpp_heatLoc = gotmkpp_heatLoc.swap_dims({"time" : "time2_heatLoc"})
gotmkpp_heatg05 = gotmkpp_heatg05.assign(time2_heatg05=("time", time2_heatg05)); gotmkpp_heatg05 = gotmkpp_heatg05.swap_dims({"time" : "time2_heatg05"})
gotmkpp_heatg2 = gotmkpp_heatg2.assign(time2_heatg2=("time", time2_heatg2)); gotmkpp_heatg2 = gotmkpp_heatg2.swap_dims({"time" : "time2_heatg2"})
gotmkpp_heat05dlta = gotmkpp_heat05dlta.assign(time2_heat05dlta=("time", time2_heat05dlta)); gotmkpp_heat05dlta = gotmkpp_heat05dlta.swap_dims({"time" : "time2_heat05dlta"})
gotmkpp_heat025dlta = gotmkpp_heat025dlta.assign(time2_heat025dlta=("time", time2_heat025dlta)); gotmkpp_heat025dlta = gotmkpp_heat025dlta.swap_dims({"time" : "time2_heat025dlta"})
gotmkpp_heat125dlta = gotmkpp_heat125dlta.assign(time2_heat125dlta=("time", time2_heat125dlta)); gotmkpp_heat125dlta = gotmkpp_heat125dlta.swap_dims({"time" : "time2_heat125dlta"})
gotmkpp_heat0625dlta = gotmkpp_heat0625dlta.assign(time2_heat0625dlta=("time", time2_heat0625dlta)); gotmkpp_heat0625dlta = gotmkpp_heat0625dlta.swap_dims({"time" : "time2_heat0625dlta"})
gotmkpp_m_h = gotmkpp_m_h.assign(time2_mh=("time", time2_mh)); gotmkpp_m_h = gotmkpp_m_h.swap_dims({"time" : "time2_mh"})
gotmkpp_m_hLoc = gotmkpp_m_hLoc.assign(time2_mhLoc=("time", time2_mhLoc)); gotmkpp_m_hLoc = gotmkpp_m_hLoc.swap_dims({"time" : "time2_mhLoc"})
gotmkpp_m_hg05 = gotmkpp_m_hg05.assign(time2_mhg05=("time", time2_mhg05)); gotmkpp_m_hg05 = gotmkpp_m_hg05.swap_dims({"time" : "time2_mhg05"})
gotmkpp_m_hg2 = gotmkpp_m_hg2.assign(time2_mhg2=("time", time2_mhg2)); gotmkpp_m_hg2 = gotmkpp_m_hg2.swap_dims({"time" : "time2_mhg2"})
gotmkpp_m_h05dlta = gotmkpp_m_h05dlta.assign(time2_mh05dlta=("time", time2_mh05dlta)); gotmkpp_m_h05dlta = gotmkpp_m_h05dlta.swap_dims({"time" : "time2_mh05dlta"})
gotmkpp_m_h025dlta = gotmkpp_m_h025dlta.assign(time2_mh025dlta=("time", time2_mh025dlta)); gotmkpp_m_h025dlta = gotmkpp_m_h025dlta.swap_dims({"time" : "time2_mh025dlta"})
gotmkpp_m_h125dlta = gotmkpp_m_h125dlta.assign(time2_mh125dlta=("time", time2_mh125dlta)); gotmkpp_m_h125dlta = gotmkpp_m_h125dlta.swap_dims({"time" : "time2_mh125dlta"})
gotmkpp_m_h0625dlta = gotmkpp_m_h0625dlta.assign(time2_mh0625dlta=("time", time2_mh0625dlta)); gotmkpp_m_h0625dlta = gotmkpp_m_h0625dlta.swap_dims({"time" : "time2_mh0625dlta"})
gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeps = gotmkeps.swap_dims({"time" : "ketime2"})

# print('GOTM time from 2012-09-17',time2.shape)
# print('GOTM time from 2012-09-17',ketime2.shape)
# print('Glider time from 2012-09-04',new_time[143:950])

# print('Glider time from high values Lat',new_time[308:383])
# print('GOTM time from high values Lat v1',time2.values, time2.shape)

# GliderTimeIndex = np.c_[np.arange(308,383),new_time[308:383]]
# print('2D GLider time with index', GliderTimeIndex)
# GTMTimeIndex = np.c_[np.arange(0,954),time2.values]
# print('2D GOTM time with index', GTMTimeIndex)

# ### ----------------------------- calculate mixed layer depth ----------------------------- ###
# ## Psedocode 
# # 1. extract temperature values for depth z_ref
# # 2. new_temp = temp(z_ref) - delta_t
# # 3. depth at new_temp

# # MLD_DT02 = depth where (T = T(z_ref) ± delta_t)
# # z_ref = -10 -- reference depth; 	T(z_ref) -- temperature at reference depth 		  	
# # delta_t = 0.2 -- ; 	

# ## using KPP model
# mld_tempKPP = GTMtemp.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp-mld_tempKPP).argmin('z'); # print(z_indexes)
# mldKPP = GTMtemp.z[z_indexes]; 
# mldKPP = mldKPP.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# # plt.plot(time2,mldKPP); plt.xlabel('time -julian day'); plt.ylabel('mixed layer depth -m')
# # plt.axis([260,341,-140,0]); plt.title('GOTM Mixed Layer Depth')
# # plt.savefig('GOTM_Mixed_Layer_Depth_test6.png'); plt.show()

# ## using KPP EV model
# mld_tempKPPLa = GTMtempLa.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLa-mld_tempKPPLa).argmin('z'); # print(z_indexes)
# mldKPPLa = GTMtempLa.z[z_indexes]; 
# mldKPPLa = mldKPPLa.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP HD model
# mld_tempKPPheat = GTMtemp_heat.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp_heat-mld_tempKPPheat).argmin('z'); # print(z_indexes)
# mldKPPheat = GTMtemp_heat.z[z_indexes]; 
# mldKPPheat = mldKPPheat.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP EV-HD model
# mld_tempKPPmh = GTMtemp_mh.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp_mh-mld_tempKPPmh).argmin('z'); # print(z_indexes)
# mldKPPmh = GTMtemp_mh.z[z_indexes]; 
# mldKPPmh = mldKPPmh.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using k-epsilon model
# mld_tempKeps = GTMketemp.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMketemp-mld_tempKeps).argmin('z'); # print(z_indexes)
# mldKeps = GTMketemp.z[z_indexes]; 
# mldKeps = mldKeps.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# ### ----------------------------- mixed layer depth .csv file ----------------------------- ###
# def read_lines():
#     with open('mld.csv', 'rU') as mldFile:
#         reader = csv.reader(mldFile)
#         for row in reader:
#             yield [ float(i) for i in row ]
# for i in read_lines():
#     print('mld :',i)

# # to get a list, instead of a generator, use
# # xy = list(read_lines())

# def read__lines():
#     with open('dats.csv', 'rU') as datFile:
#         reader_ = csv.reader(datFile)
#         for row in reader_:
#             yield [ float(j) for j in row ]
# for j in read__lines():
#     print('date :',j)
# j = np.array(j)
# j = j + 241

# # to get a list, instead of a generator, use
# # x_y = list(read__lines())
# ### Plot of MLD from OSMOSIS Observations
# # plt.plot(j,np.negative(i)); plt.xlabel('time -julian day'); plt.ylabel('mixed layer depth -m')
# # # plt.axis([260,341,-400,0]); 
# # plt.title('OSMOSIS Mixed Layer Depth')
# # plt.savefig('OSMOSIS_Mixed_Layer_Depth_.png'); plt.show()

# ### Plots-compare all MLDs ###
# plt.plot(j,np.negative(i), color='olivedrab',linestyle='dotted',label = 'OSMOSIS'); plt.legend(loc=(0.10,0.15)); 
# plt.plot(time2,mldKPP,color='gray',linestyle='solid',label = 'GOTM KPP'); plt.legend(loc=(0.10,0.15));
# plt.plot(time2La,mldKPPLa,color='steelblue',linestyle='solid',label = 'GOTM KPP EV');
# plt.plot(time2_heat,mldKPPheat,color='firebrick',linestyle='solid',label = 'GOTM KPP HD');
# plt.plot(time_mh,mldKPPmh,color='darkorchid',linestyle='solid',label = 'GOTM KPP EV-HD');
# plt.plot(time2,mldKeps,color='goldenrod',linestyle='solid',label = 'GOTM k-eps');
# plt.legend(loc='best',fontsize='small')
# plt.grid()
# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')
# plt.axis([277.0, 283.625,-150,0]);# plt.title('Compare OSMOSIS & GOTM-KPP,-k-epsilon MLDs')
# plt.savefig('GOTM_OSMOSIS_Mixed_Layer_Depth_ref11.png'); plt.show()

# ### --------------------------- quantitative analysis --------------------------- ###
# # # remove the nan values 
# GldrtempIdxNaN = np.argwhere(np.isnan(GldrTemp.values));# print('where are the nan values in glider MLD temp', tempIdxNaN)
# GldrtempIdxRowNaN = np.unique(GldrtempIdxNaN[:,0]);# print('row indices of Glider temp NaN values', tempIdxRowNaN)
# new_GldrTemp = np.delete(GldrTemp, GldrtempIdxRowNaN, axis=0);
# Glider_newtime = np.delete(new_time, GldrtempIdxRowNaN, axis=0);

# GldrtempIdxColNaN = np.unique(GldrtempIdxNaN[:,:]);# print('row indices of Glider temp NaN values', tempIdxRowNaN)
# newC_GldrTemp = np.delete(GldrTemp, GldrtempIdxColNaN, axis=1);
# GliderC_depth = np.delete(GldrDepth, GldrtempIdxColNaN, axis=0);
# # interpolate the matching times for OSMOSIS and GOTM temperature
# # GliderTempInterp = interpolate.interp2d(GldrDepth,Glider_newtime,new_GldrTemp)

# print('deleted nan in Temp col',newC_GldrTemp.values,newC_GldrTemp.shape, GliderC_depth.shape)
# print('deleted nan in Temp row',new_GldrTemp.values,new_GldrTemp.shape, Glider_newtime.shape)
# # print(new_time[309:381].shape, GldrDepth.shape, GldrTemp[309:381,:].shape )
# # print(time2.shape, GTMdepth.shape, GTMtemp.shape )
# print('Gldr depth',GldrDepth)
# print('GOTM depth',GTMdepth)

# # GliderTempInterp = np.interp(time2,new_time,GldrTemp[:,0])
# # print('interpolated Glider Temp', GliderTempInterp.shape)

# ################################################ 
# ## 1 ## interpolate Glider temperature values ##
# ################################################
# kpptemplist = []
# kpptempLalist = []
# kpptempLa125dltalist = []
# kpptemp_heatlist = []
# kpptemp_heat125dltalist = []
# kpptemp_mhlist = []
# kpptemp_mh125dltalist = []
# kepstemp1dlist =[]
# GldrTempIntlist = []

# correlation_KPP_Glider = []
# correlation_KPP_Glider = []
# correlation_KPPEV_Glider = []
# correlation_KPPEV125d_Glider = []
# correlation_KPPEV125d_KPPEV = []
# correlation_KPPHD_Glider = []
# correlation_KPPHD125d_Glider = []
# correlation_KPPEVHD_Glider = []
# correlation_KPPEVHD125d_Glider =[]
# correlation_KPPEV125d_KPPEVHD =[]
# correlation_keps_Glider = []

# mean_abs_error_KPP = []
# mean_abs_error_KPPEV = []
# mean_abs_error_KPPEV125d = []
# mean_abs_error_KPPHD = []
# mean_abs_error_KPPHD125d = []
# mean_abs_error_KPPEVHD = []
# mean_abs_error_KPPEVHD125d = []
# mean_abs_error_keps = []

# root_mean_square_error_KPP = []
# root_mean_square_error_KPPEV = []
# root_mean_square_error_KPPEV125d = []
# root_mean_square_error_KPPHD = []
# root_mean_square_error_KPPHD125d = []
# root_mean_square_error_KPPEVHD = []
# root_mean_square_error_KPPEVHD125d = []
# root_mean_square_error_keps = []

# Oz_sel = [75] #0:50:100:150:200
# Gz_sel = [424] #500:450:400:350:300:299
# # for depth in Oz_sel:
# # 	# GldrTempInterp = np.interp(time2,new_time,GldrTemp[:,depth])
# # 	GldrTempInterp = np.interp(time2,new_time,GldrTemp[:,depth])
# # 	Gldrdepth = GldrDepth[depth]; Gldrtime = time2[:]
# # 	for depthG in Gz_sel:
# # 		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[:,depthG]; 
# # 		kppdepth = GTMdepth[depthG]; kpptime = time2[:]
# # plt.plot(Gldrtime,GldrTempInterp, '-k', label='interpolated OSMOSIS Temperature at z=%0.2f'%Gldrdepth.values)
# # plt.plot(kpptime,kpptemp1d, '-b', label='GOTM KPP Temperature at z=%0.2f'%kppdepth.values)
# # plt.legend(loc='best',fontsize='small')
# # plt.xlabel('Time -Julian days'); plt.ylabel('Temperature -C')
# # plt.savefig('osmosis_autumn_variablePlots/Gldr_GOTM_Temperature_interpolated_timeseries_Aut_highLat.png');#plt.show()

# for depth in Oz_sel:
# 	# GldrTempInterp = np.interp(time2,new_time,GldrTemp[:,depth])
# 	GldrTemp = new_GldrTemp.isel()[:,depth]
# 	Gldrdepth = GldrDepth[depth]; Gldrtime = Glider_newtime[:]
# 	GldrTempIntlist = np.append(GldrTempIntlist, GldrTemp.values)
# 	for depthG in Gz_sel:
# 		kpptempinterp = np.interp(Glider_newtime,time2,GTMtemp.isel(lon=0,lat=0)[:,depthG]); kppdepth = GTMdepth[depthG]; kpptime = time2[:]
# 		kpptemplist = np.append(kpptemplist, kpptempinterp)
# 		kpptempLainterp = np.interp(Glider_newtime,time2La,GTMtempLa.isel(lon=0,lat=0)[:,depthG]); kppdepthLa = GTMdepthLa[depthG]; kpptimeLa = time2La[:]
# 		kpptempLalist = np.append(kpptempLalist,kpptempLainterp)
# 		kpptempLa125dltainterp = np.interp(Glider_newtime,time2La125dlta,GTMtempLa125dlta.isel(lon=0,lat=0)[:,depthG]); kppdepthLa125dlta = GTMdepthLa125dlta[depthG]; kpptimeLa125dlta = time2La125dlta[:]
# 		kpptempLa125dltalist = np.append(kpptempLa125dltalist,kpptempLa125dltainterp)
# 		kpptemp_heatinterp = np.interp(Glider_newtime,time2_heat,GTMtemp_heat.isel(lon=0,lat=0)[:,depthG]); kppdepth_heat = GTMdepth_heat[depthG]; kpptime_heat = time2_heat[:]
# 		kpptemp_heatlist = np.append(kpptemp_heatlist,kpptemp_heatinterp) 
# 		kpptemp_heat125dltainterp = np.interp(Glider_newtime,time2_heat125dlta,GTMtemp_heat125dlta.isel(lon=0,lat=0)[:,depthG]); kppdepth_heat125dlta = GTMdepth_heat125dlta[depthG]; kpptime_heat125dlta = time2_heat125dlta[:]
# 		kpptemp_heat125dltalist = np.append(kpptemp_heat125dltalist, kpptemp_heat125dltainterp)
# 		kpptemp_mhinterp = np.interp(Glider_newtime,time2_mh,GTMtemp_mh.isel(lon=0,lat=0)[:,depthG]); kppdepth_mh = GTMdepth_mh[depthG]; kpptime_mh = time2_mh[:]
# 		kpptemp_mhlist = np.append(kpptemp_mhlist,kpptemp_mhinterp)
# 		kpptemp_mh125dltainterp = np.interp(Glider_newtime,time2_mh125dlta,GTMtemp_mh125dlta.isel(lon=0,lat=0)[:,depthG]); kppdepth_mh125dlta = GTMdepth_mh125dlta[depthG]; kpptime_mh125dlta = time2_mh125dlta[:]
# 		kpptemp_mh125dltalist = np.append(kpptemp_mh125dltalist, kpptemp_mh125dltainterp)
# 		kepstemp1dinterp = np.interp(Glider_newtime,ketime2,GTMketemp.isel(lon=0,lat=0)[:,depthG]); kepsdepth = GTMkedepth[depthG]; kepstime = ketime2[:]
# 		kepstemp1dlist = np.append(kepstemp1dlist,kepstemp1dinterp)
# 		#analysis variables
# 		sumkppGldr = sum(kpptemplist*GldrTempIntlist);sumkpp = sum((kpptemplist)**(2)); sumGldr = sum((GldrTempIntlist)**(2))
# 	# # correlation coeef
# 	correlation_KPP_Glider = np.append(correlation_KPP_Glider, sum(kpptemplist*GldrTempIntlist)/(sum((kpptemplist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPEV_Glider = np.append(correlation_KPPEV_Glider, sum(kpptempLalist*GldrTempIntlist)/(sum((kpptempLalist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPEV125d_Glider = np.append(correlation_KPPEV125d_Glider, sum(kpptempLa125dltalist*GldrTempIntlist)/(sum((kpptempLa125dltalist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPEV125d_KPPEV = np.append(correlation_KPPEV125d_KPPEV, sum(kpptempLa125dltalist*kpptempLalist)/(sum((kpptempLa125dltalist)**(2))*sum((kpptempLalist)**(2)))**(0.5) )
# 	correlation_KPPHD_Glider = np.append(correlation_KPPHD_Glider, sum(kpptemp_heatlist*GldrTempIntlist)/(sum((kpptemp_heatlist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPHD125d_Glider = np.append(correlation_KPPHD125d_Glider, sum(kpptemp_heat125dltalist*GldrTempIntlist)/(sum((kpptemp_heat125dltalist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPEVHD_Glider = np.append(correlation_KPPEVHD_Glider, sum(kpptemp_mhlist*GldrTempIntlist)/(sum((kpptemp_mhlist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPEVHD125d_Glider = np.append(correlation_KPPEVHD125d_Glider, sum(kpptemp_mh125dltalist*GldrTempIntlist)/(sum((kpptemp_mh125dltalist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	correlation_KPPEV125d_KPPEVHD = np.append(correlation_KPPEV125d_KPPEVHD, sum(kpptemp_mh125dltalist*kpptemp_mhlist)/(sum((kpptemp_mh125dltalist)**(2))*sum((kpptemp_mhlist)**(2)))**(0.5) )
# 	correlation_keps_Glider = np.append(correlation_keps_Glider, sum(kepstemp1dlist*GldrTempIntlist)/(sum((kepstemp1dlist)**(2))*sum((GldrTempIntlist)**(2)))**(0.5) )
# 	# # mean absolute error
# 	mean_abs_error_KPP = np.append(mean_abs_error_KPP, 1/np.size(kpptemplist)*sum(np.absolute(kpptemplist-GldrTempIntlist)))
# 	mean_abs_error_KPPEV = np.append(mean_abs_error_KPPEV, 1/np.size(kpptempLalist)*sum(np.absolute(kpptempLalist-GldrTempIntlist)))
# 	mean_abs_error_KPPEV125d = np.append(mean_abs_error_KPPEV125d, 1/np.size(kpptempLa125dltalist)*sum(np.absolute(kpptempLa125dltalist-GldrTempIntlist)))
# 	mean_abs_error_KPPHD = np.append(mean_abs_error_KPPHD, 1/np.size(kpptemp_heatlist)*sum(np.absolute(kpptemp_heatlist-GldrTempIntlist)))
# 	mean_abs_error_KPPHD125d = np.append(mean_abs_error_KPPHD125d, 1/np.size(kpptemp_heat125dltalist)*sum(np.absolute(kpptemp_heat125dltalist-GldrTempIntlist)))
# 	mean_abs_error_KPPEVHD = np.append(mean_abs_error_KPPEVHD, 1/np.size(kpptemp_mhlist)*sum(np.absolute(kpptemp_mhlist-GldrTempIntlist)))
# 	mean_abs_error_KPPEVHD125d = np.append(mean_abs_error_KPPEVHD125d, 1/np.size(kpptemp_mh125dltalist)*sum(np.absolute(kpptemp_mh125dltalist-GldrTempIntlist)))
# 	mean_abs_error_keps = np.append(mean_abs_error_keps, 1/np.size(kepstemp1dlist)*sum(np.absolute(kepstemp1dlist-GldrTempIntlist)))

# 	# # rms
# 	root_mean_square_error_KPP = np.append(root_mean_square_error_KPP, (1/np.size(kpptemplist)*sum((kpptemplist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_KPPEV = np.append(root_mean_square_error_KPPEV, (1/np.size(kpptempLalist)*sum((kpptempLalist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_KPPEV125d = np.append(root_mean_square_error_KPPEV125d, (1/np.size(kpptempLa125dltalist)*sum((kpptempLa125dltalist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_KPPHD = np.append(root_mean_square_error_KPPHD, (1/np.size(kpptemp_heatlist)*sum((kpptemp_heatlist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_KPPHD125d = np.append(root_mean_square_error_KPPHD125d, (1/np.size(kpptemp_heat125dltalist)*sum((kpptemp_heat125dltalist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_KPPEVHD = np.append(root_mean_square_error_KPPEVHD, (1/np.size(kpptemp_mhlist)*sum((kpptemp_mhlist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_KPPEVHD125d = np.append(root_mean_square_error_KPPEVHD125d, (1/np.size(kpptemp_mh125dltalist)*sum((kpptemp_mh125dltalist-GldrTempIntlist)**(2)))**(0.5))
# 	root_mean_square_error_keps = np.append(root_mean_square_error_keps, (1/np.size(kepstemp1dlist)*sum((kepstemp1dlist-GldrTempIntlist)**(2)))**(0.5))
# print('KPP temp interpolated : ',kpptemplist)
# print('OSMOSIS temp',GldrTempIntlist)
# print('GOTMdepth=%0.2f'%kppdepth.values); print('Gliderdepth=%0.2f'%Gldrdepth.values)
# ##
# print('correlation coefficient between KPP and OSMOSIS Temperature :', correlation_KPP_Glider)
# print('correlation coefficient between KPP EV and OSMOSIS Temperature :', correlation_KPPEV_Glider)
# print('correlation coefficient between KPP EV 125dlta and OSMOSIS Temperature :', correlation_KPPEV125d_Glider)
# print('correlation coefficient between KPP EV 125dlta and KPP EV Temperature :', correlation_KPPEV125d_KPPEV)
# print('correlation coefficient between KPP HD and OSMOSIS Temperature :', correlation_KPPHD_Glider)
# print('correlation coefficient between KPP HD 125dlta and OSMOSIS Temperature :', correlation_KPPHD125d_Glider)
# print('correlation coefficient between KPP EV-HD and OSMOSIS Temperature :', correlation_KPPEVHD_Glider)
# print('correlation coefficient between KPP EV-HD 125dlta and OSMOSIS Temperature :', correlation_KPPEVHD125d_Glider)
# print('correlation coefficient between KPP EV-HD 125dlta and KPP EV-HD Temperature :', correlation_KPPEV125d_KPPEVHD)
# print('correlation coefficient between keps and OSMOSIS Temperature :', correlation_keps_Glider)
# ##
# print('mean absolute err between KPP and OSMOSIS Temperature :', mean_abs_error_KPP)
# print('mean absolute err between KPP EV and OSMOSIS and OSMOSIS Temperature :', mean_abs_error_KPPEV)
# print('mean absolute err between KPP EV 125dlta and OSMOSIS Temperature :', mean_abs_error_KPPEV125d)
# print('mean absolute err between KPP HD and OSMOSIS Temperature :', mean_abs_error_KPPHD)
# print('mean absolute err between KPP HD 125dlta and OSMOSIS Temperature :', mean_abs_error_KPPHD125d)
# print('mean absolute err between KPP EV-HD and OSMOSIS Temperature :', mean_abs_error_KPPEVHD)
# print('mean absolute err between KPP EV-HD 125dlta and OSMOSIS Temperature :', mean_abs_error_KPPEVHD125d)
# print('mean absolute err between keps and OSMOSIS Temperature :', mean_abs_error_keps)
# ##
# print('root mean sqare err between KPP and OSMOSIS Temperature :', root_mean_square_error_KPP)
# print('root mean sqare err between KPP EV and OSMOSIS Temperature :', root_mean_square_error_KPPEV)
# print('root mean sqare err between KPP EV 125dlta and OSMOSIS Temperature :', root_mean_square_error_KPPEV125d)
# print('root mean sqare err between KPP HD and OSMOSIS Temperature :', root_mean_square_error_KPPHD)
# print('root mean sqare err between KPP HD 125dlta and OSMOSIS Temperature :', root_mean_square_error_KPPHD125d)
# print('root mean sqare err between KPP EV-HD and OSMOSIS Temperature :', root_mean_square_error_KPPEVHD)
# print('root mean sqare err between KPP EV-HD 125dlta and OSMOSIS Temperature :', root_mean_square_error_KPPEVHD125d)
# print('root mean sqare err between keps and OSMOSIS Temperature :', root_mean_square_error_keps)
# # plt.plot(Gldrtime,GldrTemp, '-k', label='OSMOSIS Temperature at z=%0.2f'%Gldrdepth.values)
# # plt.plot(Glider_newtime,kpptempinterp, '-b', label=' interpolated GOTM KPP Temperature at z=%0.2f'%kppdepth.values)
# # plt.legend(loc='best',fontsize='small')
# # plt.xlabel('Time -Julian days'); plt.ylabel('Temperature -C')
# # plt.savefig('osmosis_autumn_variablePlots/Gldr_GOTM_Temperature_interpolated_timeseries_Aut_highLat.png');plt.show()

################################## 
## 2 ## correlation coefficient ##
##################################

# print('GOTM KPP temperature 1D : ', kpptemp1d.shape)
# print('OSMOSIS interpolated temp 1D : ', GldrTempInterp.shape)

# # kpptemp1d = np.array(list(kpptemp1d.item())) ; GldrTempInterp = np.array(list(GldrTempInterp.item()))
# # kpptemp1d = np.fromiter(kpptemp1d, np.float) ; GldrTempInterp = np.fromiter(GldrTempInterp, np.float)

# # correlation_KPP_Glider = []

# # for time in range(len(kpptemp1d)):
# 	# Gldrdepth = GldrDepth[:]; Gldrtime = time2[time]
# 	# GldrTempInterp = np.interp(time2,new_time,GldrTemp[time,:])
# 	# kpptemp1d = GTMtemp.isel(lon=0,lat=0)[time,:]; 
# 	# kppdepth = GTMdepth[:]; kpptime = time2[time]

# Ot_sel = [381] # [f:m:l]=[309:343:381] 
# Gt_sel = [476] # [f:m:l]=[0:215:476]

# for time in Ot_sel:
# 	# GldrTempInterp = np.interp(time2,new_time,GldrTemp[:,depth])
# 	GldrTemp = newC_GldrTemp.isel()[time,:]
# 	Gldrdepth = GliderC_depth[:]; Gldrtime = time2[time]
# 	GldrTempIntlist = np.append(GldrTempIntlist, GldrTemp.values)
# 	for timeG in Gt_sel:
# 		kpptempinterp = np.interp(Glider_newtime,time2,GTMtemp.isel(lon=0,lat=0)[:,timeG]); 
# 		kppdepth = GTMdepth[:]; kpptime = time2[timeG]
# 		kpptemplist = np.append(kpptemplist, kpptempinterp)
# sumkppGldr = sum(kpptemplist*GldrTempIntlist)
# sumkpp = sum((kpptemplist)**(2))
# sumGldr = sum((GldrTempIntlist)**(2))
# correlation_KPP_Glider = sumkppGldr/(sumkpp*sumGldr)**(0.5) 
# print('KPP temp interpolated : ',kpptemplist)
# print('OSMOSIS temp',GldrTempIntlist)
# print('correlation coefficient between KPP and OSMOSIS Temperature :', correlation_KPP_Glider)
# # plt.plot(correlation_KPP_Glider,Gldrdepth, '-k', label='interpolated OSMOSIS Temperature on d=%0.2f'%Gldrtime.values)
# # plt.plot(kpptemp1d,kppdepth, '-b', label='GOTM KPP Temperature on d=%0.2f'%kpptime.values)
# # plt.legend(loc='best',fontsize='small')
# # plt.ylabel('Depth -m'); plt.xlabel('Temperature -C')
# # plt.savefig('osmosis_autumn_variablePlots/Gldr_GOTM_Temperature_corrcoeff_Aut_highLat.png');plt.show()

### --------------------------- PLOTS --------------------------- ###

## -------------------------- Glider -------------------------- ##

# #choose the indices for each day
# a = [311, 314, 316, 319, 327, 329, 340, 349, 356, 358, 364, 365, 369, 378, 382] ;# interval 
# a = [311, 316, 327, 340, 356, 364, 369, 382]
# a = [309, 321, 333, 343, 353, 363, 375]
# for i in a:
# 	Gldrtemp1d = GldrTemp.isel()[i,:]; Gldrdepth = GldrDepth[:]
# 	Gldrtime = new_time[i] 
# 	plt.plot(Gldrtemp1d,Gldrdepth,label='Glider Jul_day: %s'%Gldrtime.values, linewidth=1.5);
# 	# plt.legend(loc=(0.55,0.05), fontsize = 'x-small')
# 	plt.legend(loc='best', fontsize = 'small')
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	plt.title('Glider Temperature Profiles for Autumn high La_t values')
# 	plt.axis([11,18,-150,0])
# plt.savefig('Glider_Temp_profile_highLat_v1.png')
# plt.show()

## -------------------------- k-eps -------------------------- ##
# #choose the indices for each day
# b = [0, 71, 143, 215, 287, 359, 431]
# for i in b:
# 	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
# 	kepstime = ketime2[i]
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps Day %s'%kepstime.values); 
# 	plt.legend(loc='best', fontsize='small')
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	plt.title('GOTM k-epsilon Temperature profiles for Autumn high La_t values')
# 	plt.axis([11,18,-150,0])
# plt.savefig('GOTM_keps_Temp_profile_highLat_v1.png')
# plt.show()

## -------------------------- KPP -------------------------- ##

# for i in b:
# 	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
# 	kpptime = time2[i]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_KPP Day %s'%kpptime.values);# plt.legend(loc=(0.05,0.75))
# 	plt.legend(loc='best', fontsize='small')
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	plt.title('GOTM KPP Temperature profiles for Autumn high La_t values')
# 	plt.axis([11,18,-150,0])
# plt.savefig('GOTM_kpp_Temp_profile_highLat_v1.png')
# plt.show()

## --------------------------- PLOTS- velocities --------------------------- ###
# c = [199]

# for i in c:
# 	kppu1d = GTMu.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
# 	kppuLa = GTMuLa.isel(lon=0,lat=0)[:,i]; kppdepthLa = GTMdepthLa[i]
# 	kppuLaLoc = GTMuLaLoc.isel(lon=0,lat=0)[:,i]; kppdepthLaLoc = GTMdepthLaLoc[i]
# 	kppu_heat = GTMu_heat.isel(lon=0,lat=0)[:,i]; kppdepth_heat = GTMdepth_heat[i]
# 	kppu_mh = GTMu_mh.isel(lon=0,lat=0)[:,i]; kppdepth_mh = GTMdepth_mh[i]
# 	kepsu1d = GTMkeu.isel(lon=0,lat=0)[:,i]; kepsdepth = GTMkedepth[i]
# 	kpptime = time2[:] 
# 	kpptimeLa = time2La[:]
# 	kpptimeLaLoc = time2LaLoc[:]
# 	kpptime_heat = time2_heat[:]
# 	kpptime_mh = time2_mh[:]
# 	kepstime = ketime2[:]
# 	plt.plot(kpptime,kppu1d,'k-',label='KPP');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptimeLa,kppuLa,'b-',label='KPP EV');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptimeLaLoc,kppuLaLoc,color='purple',linestyle='dashed',label='KPP EV');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_heat,kppu_heat,'g-',label='KPP HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_mh,kppu_mh,'r--',label='KPP EV-HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kepstime,kepsu1d,'k:',label='k-eps');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.legend(loc='best',fontsize ='small')
# 	plt.grid()
# 	plt.xlabel('Time -Julian days');plt.ylabel('Surface streamwise velocity -m/s'); 
# 	# plt.axis([11.0,15.5,-150,0])
# plt.savefig('osmosis_autumn_variablePlots/GOTM_streamwise_velocity_timeseries_AutumnhighLat.png')
# plt.show()

# for i in c:
# 	kppv1d = GTMv.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
# 	kppvLa = GTMvLa.isel(lon=0,lat=0)[:,i]; kppdepthLa = GTMdepthLa[i]
# 	kppvLaLoc = GTMvLaLoc.isel(lon=0,lat=0)[:,i]; kppdepthLaLoc = GTMdepthLaLoc[i]
# 	kppv_heat = GTMv_heat.isel(lon=0,lat=0)[:,i]; kppdepth_heat = GTMdepth_heat[i]
# 	kppv_mh = GTMv_mh.isel(lon=0,lat=0)[:,i]; kppdepth_mh = GTMdepth_mh[i]
# 	kepsv1d = GTMkev.isel(lon=0,lat=0)[:,i]; kepsdepth = GTMkedepth[i]
# 	kpptime = time2[:] 
# 	kpptimeLa = time2La[:]
# 	kpptimeLaLoc = time2LaLoc[:]
# 	kpptime_heat = time2_heat[:]
# 	kpptime_mh = time2_mh[:]
# 	kepstime = ketime2[:]
# 	plt.plot(kpptime,kppv1d,'k-',label='KPP');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptimeLa,kppvLa,'b-',label='KPP EV');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptimeLaLoc,kppvLaLoc,color='purple',linestyle='dashed',label='KPP EV');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_heat,kppv_heat,'g-',label='KPP HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kpptime_mh,kppv_mh,'r--',label='KPP EV-HD');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.plot(kepstime,kepsv1d,'k:',label='k-eps');# plt.legend(loc=(0.66,0.02),fontsize ='xx-small')
# 	plt.legend(loc='best',fontsize ='small')
# 	plt.grid()
# 	plt.xlabel('Time -Julian days');plt.ylabel('Surface spanwise velocity -m/s'); 
# 	# plt.axis([11.0,15.5,-150,0])
# plt.savefig('osmosis_autumn_variablePlots/GOTM_spanwise_velocity_timeseries_AutumnhighLat.png')
# plt.show()

# # phase plot == u v
# for i in c:
# 	kppu1d = GTMu.isel(lon=0,lat=0)[:,i]; kppv1d = GTMv.isel(lon=0,lat=0)[:,i];
# 	kppuLa = GTMuLa.isel(lon=0,lat=0)[:,i]; kppvLa = GTMvLa.isel(lon=0,lat=0)[:,i];
# 	kppu_heat = GTMu_heat.isel(lon=0,lat=0)[:,i]; kppv_heat = GTMv_heat.isel(lon=0,lat=0)[:,i];
# 	kppu_mh = GTMu_mh.isel(lon=0,lat=0)[:,i]; kppv_mh = GTMv_mh.isel(lon=0,lat=0)[:,i];
# 	# kppnuh_heatg5 = GTMnuh_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
# 	# kppnuh_heatg15 = GTMnuh_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
# 	kepsu1d = GTMkeu.isel(lon=0,lat=0)[:,i]; kepsv1d = GTMkeu.isel(lon=0,lat=0)[:,i];
# 	kpptime = time2[i]
# 	kpptimeLa = time2La[i]
# 	kpptime_heat = time2_heat[i]
# 	kpptime_mh = time2_mh[i]
# 	# kpptime_heatg5 = time2_heatg5[i]
# 	# kpptime_heatg15 = time2_heatg15[i]
# 	kepstime = ketime2[i]
# plt.plot(kppu1d,kppv1d,'k-.',label = 'KPP'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppuLa,kppvLa,'b-' ,label = 'KPP EV'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_heat,kppv_heat,'g-',label = 'KPP HD'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppu_mh,kppv_mh,'r--',label = 'KPP EV-HD'); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_heatg5,kppdepth_heatg5,'r-',label = 'KPP diff gam=0.5 day %0.2f'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_heatg15,kppdepth_heatg15,'r--',label = 'KPP diff gam=1.5 day %0.2f'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kepsu1d,kepsv1d,'k:' ,label = 'k-eps'); plt.legend(loc='best', fontsize='small')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.grid()
# plt.xlabel('Streamwise velocity -m/s');plt.ylabel('Spanwise velocity -m/s'); 
# # plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# # plt.axis([-0.01,0.15,-100,0])
# # plt.ylim(-150,0)
# plt.savefig('autumn_highLatKPP/GOTM_kpp_kppLa_kppheat_kppmh_keps_velocity_phase_plots.png')
# plt.show()

## -------------------------- Comparisons ------------------------------ ##
# select desired rows for the loop to plot day after change in wave direction
# print('Glider time from high values Lat',new_time[286:380])
# print('Glider high values Lat wind change 1',new_time[339])
# print('GOTM high values Lat wind change 1',time2[188])
# print('Glider high values Lat wind change 2',new_time[315])
# print('GOTM high values Lat wind change 2',time2[35])

Gldr_sel = [309] # range [f:m:l]=[309:343:381] 
GTM_sel = [0] # range [f:m:l]=[0:431:953]

# for i in GTM_sel:
# 	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
# 	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
# 	kpptime = time2[i]
# 	kepstime = ketime2[i]
# 	plt.plot(kpptemp1d,kppdepth,label = 'GOTM_KPP Day %s'%kpptime.values); plt.legend(loc=(0.05,0.75))
# 	plt.plot(kepstemp1d,kepsdepth,label = 'GOTM_k-eps Day %s'%kepstime.values); plt.legend(loc=(0.05,0.75))
# 	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
# 	plt.title('Compared GOTM Temperature Profiles : KPP and k-eps models')
# 	plt.axis([11,18,-100,0])
# plt.savefig('GOTM_kppkeps_Temp_profile_highLat_test.png')
# plt.show()

# Compared plot GLider and GOTM_KPP
for j in Gldr_sel:
	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
	Gldrtime = new_time[j] 
	for i in GTM_sel:
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
		# kpptempLoc = GTMtempLoc.isel(lon=0,lat=0)[i,:]; kppdepthLoc = GTMdepthLoc[:]
		kpptempLa = GTMtempLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
		# kpptempLaLoc = GTMtempLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
		kpptempLag05 = GTMtempLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
		kpptempLag2 = GTMtempLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
		kpptempLa05dlta = GTMtempLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
		kpptempLa025dlta = GTMtempLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
		kpptempLa125dlta = GTMtempLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
		kpptempLa0625dlta = GTMtempLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
		kpptemp_heat = GTMtemp_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
		# kpptemp_heatLoc = GTMtemp_heatLoc.isel(lon=0,lat=0)[i,:]; kppdepth_heatLoc = GTMdepth_heatLoc[:]
		kpptemp_heatg05 = GTMtemp_heatg05.isel(lon=0,lat=0)[i,:]; kppdepth_heatg05 = GTMdepth_heatg05[:]
		kpptemp_heatg2 = GTMtemp_heatg2.isel(lon=0,lat=0)[i,:]; kppdepth_heatg2 = GTMdepth_heatg2[:]
		kpptemp_heat05dlta = GTMtemp_heat05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat05dlta = GTMdepth_heat05dlta[:]
		kpptemp_heat025dlta = GTMtemp_heat025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat025dlta = GTMdepth_heat025dlta[:]
		kpptemp_heat125dlta = GTMtemp_heat125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat125dlta = GTMdepth_heat125dlta[:]
		kpptemp_heat0625dlta = GTMtemp_heat0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_heat0625dlta = GTMdepth_heat0625dlta[:]
		kpptemp_mh = GTMtemp_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
		# kpptemp_mhLoc = GTMtemp_mhLoc.isel(lon=0,lat=0)[i,:]; kppdepth_mhLoc = GTMdepth_mhLoc[:]
		kpptemp_mhg05 = GTMtemp_mhg05.isel(lon=0,lat=0)[i,:]; kppdepth_mhg05 = GTMdepth_mhg05[:]
		kpptemp_mhg2 = GTMtemp_mhg2.isel(lon=0,lat=0)[i,:]; kppdepth_mhg2 = GTMdepth_mhg2[:]
		kpptemp_mh05dlta = GTMtemp_mh05dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh05dlta = GTMdepth_mh05dlta[:]
		kpptemp_mh025dlta = GTMtemp_mh025dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh025dlta = GTMdepth_mh025dlta[:]
		kpptemp_mh125dlta = GTMtemp_mh125dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh125dlta = GTMdepth_mh125dlta[:]
		kpptemp_mh0625dlta = GTMtemp_mh0625dlta.isel(lon=0,lat=0)[i,:]; kppdepth_mh0625dlta = GTMdepth_mh0625dlta[:]
		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
		kpptime = time2[i]
		# kpptimeLoc = time2Loc[i]
		kpptimeLa = time2La[i]
		# kpptimeLaLoc = time2LaLoc[i]
		kpptimeLag05 = time2Lag05[i]
		kpptimeLag2 = time2Lag2[i]
		kpptimeLa05dlta = time2La05dlta[i]
		kpptimeLa025dlta = time2La025dlta[i]
		kpptimeLa125dlta = time2La125dlta[i]
		kpptimeLa0625dlta = time2La0625dlta[i]
		kpptime_heat = time2_heat[i]
		# kpptime_heatLoc = time2_heatLoc[i]
		kpptime_heatg05 = time2_heatg05[i]
		kpptime_heatg2 = time2_heatg2[i]
		kpptime_heat05dlta = time2_heat05dlta[i]
		kpptime_heat025dlta = time2_heat025dlta[i]
		kpptime_heat125dlta = time2_heat125dlta[i]		
		kpptime_heat0625dlta = time2_heat0625dlta[i]		
		kpptime_mh = time2_mh[i]
		# kpptime_mhLoc = time2_mhLoc[i]
		kpptime_mhg05 = time2_mhg05[i]
		kpptime_mhg2 = time2_mhg2[i]
		kpptime_mh05dlta = time2_mh05dlta[i]
		kpptime_mh025dlta = time2_mh025dlta[i]
		kpptime_mh125dlta = time2_mh125dlta[i]
		kpptime_mh0625dlta = time2_mh0625dlta[i]
		kepstime = ketime2[i]
	plt.plot(kpptemp1d,kppdepth,color='gray',linestyle='solid',label = 'KPP d=%0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLoc,kppdepthLoc,color='gray',linestyle='dashed',label = 'KPP local day %0.2f'%kpptime.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
	# # plt.plot(kpptempLaLoc,kppdepthLaLoc,color='steelblue',linestyle='dashed',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat,kppdepth_heat,color='firebrick',linestyle='solid',label = 'KPP HD d=%0.2f'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptemp_heatLoc,kppdepth_heatLoc,color='firebrick',linestyle='dashed',label = 'KPP HD loc d=%0.2f'%kpptime_heatLoc.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heatg05,kppdepth_heatg05,color='firebrick',linestyle=(0, (1, 10)),label = 'KPP HD g=0.5 d=%0.2f'%kpptime_heatg05.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heatg2,kppdepth_heatg2,color='firebrick',linestyle='dashdot',label = 'KPP HD g=2 d=%0.2f'%kpptime_heatg2.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat05dlta,kppdepth_heat05dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1)),label = 'KPP HD 2*dlt d=%0.2f'%kpptime_heat05dlta.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat025dlta,kppdepth_heat025dlta,color='firebrick',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP HD 4*dlt d=%0.2f'%kpptime_heat025dlta.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat125dlta,kppdepth_heat125dlta,color='firebrick',linestyle=(0, (1, 1)),label = 'KPP HD 8*dlt d=%0.2f'%kpptime_heat125dlta.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat0625dlta,kppdepth_heat0625dlta,color='firebrick',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP HD 16*dlt d=%0.2f'%kpptime_heat0625dlta.values); plt.legend(loc='best', fontsize='small')
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
	plt.axis([11,18,-150,0])
	# plt.ylim(-150,0)
plt.savefig('osmosis_autumn_variablePlots/Gldr_GOTM_Temperature_profile_Aut_highLat.png')
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
# 	plt.axis([11,18,-100,0])
# plt.savefig('GOTM_kppkeps_Temp_profile_lowLat_test.png')
# plt.show()

## EDDY VISCOSITY
for i in GTM_sel:
	kppnum1d = GTMnum.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kppnumLa = GTMnumLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppnumLaLoc = GTMnumLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	kppnumLag05 = GTMnumLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	kppnumLag2 = GTMnumLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	kppnumLa05dlta = GTMnumLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	kppnumLa025dlta = GTMnumLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	kppnumLa125dlta = GTMnumLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	kppnumLa0625dlta = GTMnumLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
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
	kpptimeLa = time2La[i]
	# kpptimeLaLoc = time2LaLoc[i]
	kpptimeLag05 = time2Lag05[i]
	kpptimeLag2 = time2Lag2[i]
	kpptimeLa05dlta = time2La05dlta[i]
	kpptimeLa025dlta = time2La025dlta[i]
	kpptimeLa125dlta = time2La125dlta[i]
	kpptimeLa0625dlta = time2La0625dlta[i]
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
plt.plot(kppnumLa,kppdepthLa,color='steelblue',linestyle='solid' ,label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
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
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.grid()
plt.xlabel('Eddy Viscosity -m2/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-150,0)
plt.savefig('osmosis_autumn_variablePlots/GOTM_eddy_viscosity_profile_Aut_highLat.png')
plt.show()

## HEAT DIFFUSIVITY
for i in GTM_sel:
	kppnuh1d = GTMnuh.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	# kppnuh1 = GTMnuh1.isel(lon=0,lat=0)[i,:]; kppdepth1 = GTMdepth1[:]
	kppnuhLa = GTMnuhLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppnuhLaLoc = GTMnuhLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	kppnuhLag05 = GTMnuhLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	kppnuhLag2 = GTMnuhLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	kppnuhLa05dlta = GTMnuhLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	kppnuhLa025dlta = GTMnuhLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	kppnuhLa125dlta = GTMnuhLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	kppnuhLa0625dlta = GTMnuhLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
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
	kpptimeLa = time2La[i]
	kpptimeLag05 = time2Lag05[i]
	kpptimeLag2 = time2Lag2[i]
	kpptimeLa05dlta = time2La05dlta[i]
	kpptimeLa025dlta = time2La025dlta[i]
	kpptimeLa125dlta = time2La125dlta[i]
	kpptimeLa0625dlta = time2La0625dlta[i]
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
plt.plot(kppnuhLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuhLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
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
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-150,0)
plt.savefig('osmosis_autumn_variablePlots/GOTM_heat_diffusivity_profile_Aut_highLat.png')
plt.show()

## Velocity profiles

# U
for i in GTM_sel:
	kppu1d = GTMu.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kppuLa = GTMuLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppuLaLoc = GTMuLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	kppuLag05 = GTMuLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	kppuLag2 = GTMuLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	kppuLa05dlta = GTMuLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	kppuLa025dlta = GTMuLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	kppuLa125dlta = GTMuLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	kppuLa0625dlta = GTMuLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
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
	kpptimeLa = time2La[i]
	kpptimeLag05 = time2Lag05[i]
	kpptimeLag2 = time2Lag2[i]
	kpptimeLa05dlta = time2La05dlta[i]
	kpptimeLa025dlta = time2La025dlta[i]
	kpptimeLa125dlta = time2La125dlta[i]
	kpptimeLa0625dlta = time2La0625dlta[i]
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
plt.plot(kppuLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppuLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppuLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppuLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppuLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppuLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppuLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppuLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
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
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-150,0)
plt.savefig('osmosis_autumn_variablePlots/GOTM_streamwise_velocity_profile_Aut_highLat.png')
plt.show()

# V
for i in GTM_sel:
	kppv1d = GTMv.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kppvLa = GTMvLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
	# kppuLaLoc = GTMvLaLoc.isel(lon=0,lat=0)[i,:]; kppdepthLaLoc = GTMdepthLaLoc[:]
	kppvLag05 = GTMvLag05.isel(lon=0,lat=0)[i,:]; kppdepthLag05 = GTMdepthLag05[:]
	kppvLag2 = GTMvLag2.isel(lon=0,lat=0)[i,:]; kppdepthLag2 = GTMdepthLag2[:]
	kppvLa05dlta = GTMvLa05dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa05dlta = GTMdepthLa05dlta[:]
	kppvLa025dlta = GTMvLa025dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa025dlta = GTMdepthLa025dlta[:]
	kppvLa125dlta = GTMvLa125dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa125dlta = GTMdepthLa125dlta[:]
	kppvLa0625dlta = GTMvLa0625dlta.isel(lon=0,lat=0)[i,:]; kppdepthLa0625dlta = GTMdepthLa0625dlta[:]
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
	kpptimeLa = time2La[i]
	kpptimeLag05 = time2Lag05[i]
	kpptimeLag2 = time2Lag2[i]
	kpptimeLa05dlta = time2La05dlta[i]
	kpptimeLa025dlta = time2La025dlta[i]
	kpptimeLa125dlta = time2La125dlta[i]
	kpptimeLa0625dlta = time2La0625dlta[i]
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
plt.plot(kppvLa,kppdepthLa,color='steelblue',linestyle='solid',label = 'KPP EV d=%0.2f'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppvLaLoc,kppdepthLaLoc,color='purple',label = 'KPP EV loc d=%0.2f'%kpptimeLaLoc.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppvLag05,kppdepthLag05,color='steelblue',linestyle=(0, (1, 10)),label = 'KPP EV g=0.5 d=%0.2f'%kpptimeLag05.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppvLag2,kppdepthLag2,color='steelblue',linestyle='dashdot',label = 'KPP EV g=2 d=%0.2f'%kpptimeLag2.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppvLa05dlta,kppdepthLa05dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1)),label = 'KPP EV 2*dlt d=%0.2f'%kpptimeLa05dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppvLa025dlta,kppdepthLa025dlta,color='steelblue',linestyle=(0, (3, 5, 1, 5, 1, 5)),label = 'KPP EV 4*dlt d=%0.2f'%kpptimeLa025dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppvLa125dlta,kppdepthLa125dlta,color='steelblue',linestyle=(0, (1, 1)),label = 'KPP EV 8*dlt d=%0.2f'%kpptimeLa125dlta.values); plt.legend(loc='best', fontsize='small')
plt.plot(kppvLa0625dlta,kppdepthLa0625dlta,color='steelblue',linestyle=(0, (3, 1, 1, 1, 1, 1)),label = 'KPP EV 16*dlt d=%0.2f'%kpptimeLa0625dlta.values); plt.legend(loc='best', fontsize='small')
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
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-150,0)
plt.savefig('osmosis_autumn_variablePlots/GOTM_spanwise_velocity_profile_Aut_highLat.png')
plt.show()


exit()
# Compared plot GLider and GOTM_KPP
for j in Gldr_sel:
	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
	Gldrtime = new_time[j] 
	for i in GTM_sel:
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
		# kpptempLoc = GTMtempLoc.isel(lon=0,lat=0)[i,:]; kppdepthLoc = GTMdepthLoc[:]
		kpptempLa = GTMtempLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
		kpptemp_heat = GTMtemp_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
		kpptemp_mh = GTMtemp_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
		kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
		kpptime = time2[i]
		kpptimeLa = time2La[i]
		kpptime_heat = time2_heat[i]
		kpptime_mh = time2_mh[i]
		kepstime = ketime2[i]
	plt.plot(kpptemp1d,kppdepth,'r-',label = 'GOTM_KPP Day %s'%kpptime.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kpptempLoc,kppdepthLoc,'k--' ,label = 'GOTM_KPP local Day %s'%kpptime.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptempLa,kppdepthLa,'b--' ,label = 'GOTM_KPP visc Day %s'%kpptimeLa.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_heat,kppdepth_heat,'b-',label = 'GOTM_KPP diff day %s'%kpptime_heat.values); plt.legend(loc='best', fontsize='small')
	plt.plot(kpptemp_mh,kppdepth_mh,'k:',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
	# plt.plot(kepstemp1d,kepsdepth,'k:' ,label = 'GOTM_k-eps Day %s'%kepstime.values); plt.legend(loc='best', fontsize='small')
	plt.plot(Gldrtemp1d,Gldrdepth,'g-' ,label='Glider Day %s'%Gldrtime.values); plt.legend(loc='best', fontsize='x-small')
	plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
	# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
	plt.axis([11,18,-150,0])
plt.savefig('autumn_highLatKPP/Gldr_GOTM_kpp_kpploc_kppLa_kppmh_keps_Temp_profile_highLat_v1.png')
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
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-100,0)
plt.savefig('autumn_highLatKPP/GOTM_kpp_kpploc_kppvisc_kppheat_kppmh_keps_eddy_viscosity_profile.png')
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
plt.xlabel('Eddy Diffusivity -m2/s');plt.ylabel('Depth -m'); 
# plt.title('Compared Temperature Profiles for high La_t values : Glider & GOTM-KPP,-k-eps')
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-100,0)
plt.savefig('autumn_highLatKPP/GOTM_kpp_kpploc_kppLa_kppheat_kppmh_keps_heat_diffusivity_profile.png')
plt.show()
exit()

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
#print('Converted Days: ', convertedDays)

## CALCULATIONS

hflux = np.genfromtxt('heat_flux.dat')
mflux = np.genfromtxt('momentum_flux.dat')
# print('Heat Flux Data: ', hflux); print(hflux.shape)


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
# print('Stoke drift veloctiy :',u_s); print(u_s.shape)
# print('friction velocity :',fric_u); print(fric_u.shape)

## calculate la_t

la_t = np.empty((634,0), float)
for i in range(len(u_s)):
	la_t = np.append(la_t, math.sqrt(fric_u[i]/u_s[i]))
# print('turbulent Langmuir number :',la_t); print(la_t.shape)

## mean value of la_t

la_tMean = np.mean(la_t)
# print('mean value of La_t', la_tMean)

##  values of la_t
la_tt = np.c_[convertedDays,la_t]; print('added day of year to la_t', la_tt)
print(la_tt.shape)
# la_tMax = np.amax(la_tt); print('maximum value of la_t',la_tMax)
# la_tMin = np.amin(la_tt); print('minimum value of la_t',la_tMin)

exit()
## ------------ period with low La_t values ------------ ##
rows, cols = np.where(la_tt < 0.6); print('values below la_t=0.6',la_tt[rows])
la_tL = la_tt[rows][:,0]; print('times of low la_t values: ', la_tL);
rowvalL = np.in1d(time2, la_tL).nonzero()[0]; rowval_low = rowvalL.tolist(); print(rowval_low)
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

## Individual plots

# rvL_kpp = list(range(0,395))

for i in range(len(GTMtemp)):
	kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
	kpptime = time2[i]; 
	plt.plot(kpptemp1d,kppdepth)#,label = 'GOTM_kpp day %s'%kpptime.values); plt.legend(loc=(0.45,0.35))
plt.xlabel('Temp -C'); plt.ylabel('Depth -m'); plt.title('GOTM_kpp_test Low La_t values')
plt.axis([15,18,-100,0])
plt.savefig('GOTM_kpp_Temp_lowLat_test.png')
plt.show()

exit()
for i in range(len(GTMketemp)):
	kepstemp1d = GTMketemp.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
	kepstime = ketime2[i]
	plt.plot(kepstemp1d,kepsdepth)#,label = 'GOTM_k-eps day %s'%kepstime.values); plt.legend(loc=(0.45,0.35))
plt.xlabel('Temp -C'); plt.ylabel('Depth -m'); plt.title('GOTM_keps_test Low La_t values')
plt.axis([11,18,-400,0])
plt.savefig('GOTM_keps_Temp_lowLat_test.png')
plt.show()


Gldr_sel = list(range(144, 202)) # range [143:950]

for j in Gldr_sel:
	Gldrtemp1d = GldrTemp.isel()[j,:]; Gldrdepth = GldrDepth[:]
	Gldrtime = new_time[j] 
	plt.plot(Gldrtemp1d,Gldrdepth)#,label='Glider Day %s'%Gldrtime.values); plt.legend(loc=(0.45,0.35))
plt.xlabel('Temperature -C');plt.ylabel('Depth -m'); 
plt.title('Glider Temperature Profiles')
plt.axis([11,18,-400,0])
plt.savefig('Gldr__Temp_profile_lowLat_test.png')
plt.show()