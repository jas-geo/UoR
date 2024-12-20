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
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
# plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

### ----------------------------- GLIDER ----------------------------- ##
# file1 = 'C:/home/users/mc837749/Documents/gotm-4.0.0/simulations/annualOSMOSISforcing_2502/glider_timeseries.nc'
# gliderdata= Dataset('../glider_timeseries.nc'); 
gliderdata= xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc'); 
glider = xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen
# print('original Gldr temp', GldrTemp.shape)
# GldrTempNaN = GldrTemp.isnull(); print('number of missing values in GldrTemp', GldrTempNaN)

# GldrTemp = GldrTemp.to_dataframe().dropna(how='all'); print('Gldr temp with dropped NaN values?', GldrTemp, GldrTemp.shape)
# GldrTemp = GldrTemp.to_dataframe();

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)

new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days

glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})
print(gliderdata.variables.items())
### ----------------------------- GOTM ----------------------------- ###

gotmkpp = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_placebo1.nc", decode_times=False)
gotmkppL = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_placebo1_LOCAL.nc", decode_times=False)
# gotmkeps = xr.open_dataset("osmosis_annual_surface_forcingkeps/OSMOSIS_glider_comparison.nc", decode_times=False)
# changed the definition of eddy viscosity (num(k)) in kpp.F90
# gotmkppLa = xr.open_dataset("../../gotm-4.0.0-edit/simulations/sim485_633/OSMOSIS_KPPEV_PAP.nc",decode_times=False)
gotmkppLaRib = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib.nc",decode_times=False)
gotmkppLaRibL = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_LOCAL.nc",decode_times=False)
# gotmkppLaRibnum2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_num2.nc",decode_times=False)
# gotmkppLaRibnuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_nuh2.nc",decode_times=False)
# gotmkppLaRibnum2nuh2 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_num2_nuh2.nc",decode_times=False)
# gotmkppLaRibhL = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_hLtest.nc",decode_times=False)
# gotmkpphLPve = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Rib_orig_hL_alpha_Pve.nc",decode_times=False)
# gotmkpphLPveN100 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Rib_orig_hL_alpha_PveN100.nc",decode_times=False)
# gotmkpphLPveP100 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Rib_orig_hL_alpha_PveP100.nc",decode_times=False)

# gotmkppLaRibg1a0 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Rib_LaT_g1a0.nc",decode_times=False)
# gotmkppLaRibg0a1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Rib_LaT_g0a1.nc",decode_times=False)

# gotmkppLaRibnum15 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_num15.nc",decode_times=False)
# gotmkppLaRibnuh15 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_nuh15.nc",decode_times=False)
# gotmkppLaRibnum15nuh15 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_num15_nuh15.nc",decode_times=False)

# gotmkppLaRibhLNve = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLNve.nc",decode_times=False)
# gotmkppLaRibhLPve = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPve.nc",decode_times=False)
# gotmkppLaRibLaTPve03 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_LaTPve03.nc",decode_times=False)
# gotmkppLaRibhLPveN100 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveN100.nc",decode_times=False)
# gotmkppLaRibhLPveP100 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP100.nc",decode_times=False)
# gotmkppLaRibhLPveP200 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP200.nc",decode_times=False)
# gotmkppLaRibhLPveP500 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP500.nc",decode_times=False)
# gotmkppLaRibhLPveP1000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP1000.nc",decode_times=False)
# gotmkppLaRibhLPveP2000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP2000.nc",decode_times=False)
gotmkppLaRibhLPveP3000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP3000.nc",decode_times=False)
gotmkppLaRibhLPveP3000L = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP3000_LOCAL.nc",decode_times=False)
gotmkppLaRibhLPveP3500 = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP3500.nc",decode_times=False)
gotmkppLaRibhLPveP3500L = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP3500_LOCAL.nc",decode_times=False)
gotmkppLaRibhLPveP4000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP4000.nc",decode_times=False)
gotmkppLaRibhLPveP4000L = xr.open_dataset("../gotm-4.0.0-edit/simulations/sim480_633/OSMOSIS_nu_Lat_Rib_hLPveP4000_LOCAL.nc",decode_times=False)
# gotmkppLaRibhLPveP5000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP5000.nc",decode_times=False)
# gotmkppLaRibhLPveP10000 = xr.open_dataset("../gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat_Rib_hLPveP10000.nc",decode_times=False)


GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMsal = gotmkpp.salt; GTMo2 = gotmkpp.o2_obs; GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMxflx = gotmkpp.gamu; GTMyflx = gotmkpp.gamv; time = gotmkpp.time
GTMtempL = gotmkppL.temp; GTMdepthL = gotmkppL.z; GTMnumL = gotmkppL.num; GTMnuhL = gotmkppL.nuh; GTMsalL = gotmkppL.salt; GTMo2L = gotmkppL.o2_obs; GTMuL = gotmkppL.u; GTMvL = gotmkppL.v; GTMxflxL = gotmkppL.gamu; GTMyflxL = gotmkppL.gamv; timeL = gotmkppL.time
# GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time
# GTMtempLa = gotmkppLa.temp; GTMdepthLa = gotmkppLa.z; GTMnumLa = gotmkppLa.num; GTMnuhLa = gotmkppLa.nuh; GTMuLa = gotmkppLa.u; GTMvLa = gotmkppLa.v; timeLa = gotmkppLa.time
GTMtempLaRib = gotmkppLaRib.temp; GTMdepthLaRib = gotmkppLaRib.z; GTMnumLaRib = gotmkppLaRib.num; GTMnuhLaRib = gotmkppLaRib.nuh; GTMuLaRib = gotmkppLaRib.u; GTMvLaRib = gotmkppLaRib.v; timeLaRib = gotmkppLaRib.time
GTMtempLaRibL = gotmkppLaRibL.temp; GTMdepthLaRibL = gotmkppLaRibL.z; GTMnumLaRibL = gotmkppLaRibL.num; GTMnuhLaRibL = gotmkppLaRibL.nuh; GTMuLaRibL = gotmkppLaRibL.u; GTMvLaRibL = gotmkppLaRibL.v; timeLaRibL = gotmkppLaRibL.time
# GTMtempLaRibnum2 = gotmkppLaRibnum2.temp; GTMdepthLaRibnum2 = gotmkppLaRibnum2.z; GTMnumLaRibnum2 = gotmkppLaRibnum2.num; GTMnuhLaRibnum2 = gotmkppLaRibnum2.nuh; GTMuLaRibnum2 = gotmkppLaRibnum2.u; GTMvLaRibnum2 = gotmkppLaRibnum2.v; timeLaRibnum2 = gotmkppLaRibnum2.time
# GTMtempLaRibnuh2 = gotmkppLaRibnuh2.temp; GTMdepthLaRibnuh2 = gotmkppLaRibnuh2.z; GTMnumLaRibnuh2 = gotmkppLaRibnuh2.num; GTMnuhLaRibnuh2 = gotmkppLaRibnuh2.nuh; GTMuLaRibnuh2 = gotmkppLaRibnuh2.u; GTMvLaRibnuh2 = gotmkppLaRibnuh2.v; timeLaRibnuh2 = gotmkppLaRibnuh2.time
# GTMtempLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.temp; GTMdepthLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.z; GTMnumLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.num; GTMnuhLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.nuh; GTMuLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.u; GTMvLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.v; timeLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.time
# GTMtempLaRibhL = gotmkppLaRibhL.temp; GTMdepthLaRibhL = gotmkppLaRibhL.z; GTMnumLaRibhL = gotmkppLaRibhL.num; GTMnuhLaRibhL = gotmkppLaRibhL.nuh; GTMuLaRibhL = gotmkppLaRibhL.u; GTMvLaRibhL = gotmkppLaRibhL.v; timeLaRibhL = gotmkppLaRibhL.time

# GTMtemphLPve = gotmkpphLPve.temp; GTMdepthhLPve = gotmkpphLPve.z; GTMnumhLPve = gotmkpphLPve.num; GTMnuhhLPve = gotmkpphLPve.nuh; GTMuhLPve = gotmkpphLPve.u; GTMvhLPve = gotmkpphLPve.v; timehLPve = gotmkpphLPve.time
# GTMtemphLPveN100 = gotmkpphLPveN100.temp; GTMdepthhLPveN100 = gotmkpphLPveN100.z; GTMnumPveN100 = gotmkpphLPveN100.num; GTMnuhPveN100 = gotmkpphLPveN100.nuh; GTMuhLPveN100 = gotmkpphLPveN100.u; GTMvhLPveN100 = gotmkpphLPveN100.v; timehLPveN100 = gotmkpphLPveN100.time
# GTMtemphLPveP100 = gotmkpphLPveP100.temp; GTMdepthhLPveP100 = gotmkpphLPveP100.z; GTMnumhLPveP100 = gotmkpphLPveP100.num; GTMnuhhLPveP100 = gotmkpphLPveP100.nuh; GTMuhLPveP100 = gotmkpphLPveP100.u; GTMvhLPveP100 = gotmkpphLPveP100.v; timehLPveP100 = gotmkpphLPveP100.time

# GTMtempLaRibg0a1 = gotmkppLaRibg0a1.temp; GTMdepthLaRibg0a1 = gotmkppLaRibg0a1.z; GTMnumLaRibg0a1 = gotmkppLaRibg0a1.num; GTMnuhLaRibg0a1 = gotmkppLaRibg0a1.nuh; GTMuLaRibg0a1 = gotmkppLaRibg0a1.u; GTMvLaRibg0a1 = gotmkppLaRibg0a1.v; timeLaRibg0a1 = gotmkppLaRibg0a1.time
# GTMtempLaRibg1a0 = gotmkppLaRibg1a0.temp; GTMdepthLaRibg1a0 = gotmkppLaRibg1a0.z; GTMnumLaRibg1a0 = gotmkppLaRibg1a0.num; GTMnuhLaRibg1a0 = gotmkppLaRibg1a0.nuh; GTMuLaRibg1a0 = gotmkppLaRibg1a0.u; GTMvLaRibg1a0 = gotmkppLaRibg1a0.v; timeLaRibg1a0 = gotmkppLaRibg1a0.time

# GTMtempLaRibnum15 = gotmkppLaRibnum15.temp; GTMdepthLaRibnum15 = gotmkppLaRibnum15.z; GTMnumLaRibnum15 = gotmkppLaRibnum15.num; GTMnuhLaRibnum15 = gotmkppLaRibnum15.nuh; GTMuLaRibnum15 = gotmkppLaRibnum15.u; GTMvLaRibnum15 = gotmkppLaRibnum15.v; timeLaRibnum15 = gotmkppLaRibnum15.time
# GTMtempLaRibnuh15 = gotmkppLaRibnuh15.temp; GTMdepthLaRibnuh15 = gotmkppLaRibnuh15.z; GTMnumLaRibnuh15 = gotmkppLaRibnuh15.num; GTMnuhLaRibnuh15 = gotmkppLaRibnuh15.nuh; GTMuLaRibnuh15 = gotmkppLaRibnuh15.u; GTMvLaRibnuh15 = gotmkppLaRibnuh15.v; timeLaRibnuh15 = gotmkppLaRibnuh15.time
# GTMtempLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.temp; GTMdepthLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.z; GTMnumLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.num; GTMnuhLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.nuh; GTMuLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.u; GTMvLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.v; timeLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.time

# GTMtempLaRibhLNve = gotmkppLaRibhLNve.temp; GTMdepthLaRibhLNve = gotmkppLaRibhLNve.z; GTMnumLaRibhLNve = gotmkppLaRibhLNve.num; GTMnuhLaRibhLNve = gotmkppLaRibhLNve.nuh; GTMuLaRibhLNve = gotmkppLaRibhLNve.u; GTMvLaRibhLNve = gotmkppLaRibhLNve.v; timeLaRibhLNve = gotmkppLaRibhLNve.time
# GTMtempLaRibhLPve = gotmkppLaRibhLPve.temp; GTMdepthLaRibhLPve = gotmkppLaRibhLPve.z; GTMnumLaRibhLPve = gotmkppLaRibhLPve.num; GTMnuhLaRibhLPve = gotmkppLaRibhLPve.nuh; GTMuLaRibhLPve = gotmkppLaRibhLPve.u; GTMvLaRibhLPve = gotmkppLaRibhLPve.v; timeLaRibhLPve = gotmkppLaRibhLPve.time
# GTMtempLaRibLaTPve03 = gotmkppLaRibLaTPve03.temp; GTMdepthLaRibLaTPve03 = gotmkppLaRibLaTPve03.z; GTMnumLaRibLaTPve03 = gotmkppLaRibLaTPve03.num; GTMnuhLaRibLaTPve03 = gotmkppLaRibLaTPve03.nuh; GTMuLaRibLaTPve03 = gotmkppLaRibLaTPve03.u; GTMvLaRibLaTPve03 = gotmkppLaRibLaTPve03.v; timeLaRibLaTPve03 = gotmkppLaRibLaTPve03.time
# GTMtempLaRibhLPveN100 = gotmkppLaRibhLPveN100.temp; GTMdepthLaRibhLPveN100 = gotmkppLaRibhLPveN100.z; GTMnumLaRibhLPveN100 = gotmkppLaRibhLPveN100.num; GTMnuhLaRibhLPveN100 = gotmkppLaRibhLPveN100.nuh; GTMuLaRibhLPveN100 = gotmkppLaRibhLPveN100.u; GTMvLaRibhLPveN100 = gotmkppLaRibhLPveN100.v; timeLaRibhLPveN100 = gotmkppLaRibhLPveN100.time
# GTMtempLaRibhLPveP100 = gotmkppLaRibhLPveP100.temp; GTMdepthLaRibhLPveP100 = gotmkppLaRibhLPveP100.z; GTMnumLaRibhLPveP100 = gotmkppLaRibhLPveP100.num; GTMnuhLaRibhLPveP100 = gotmkppLaRibhLPveP100.nuh; GTMuLaRibhLPveP100 = gotmkppLaRibhLPveP100.u; GTMvLaRibhLPveP100 = gotmkppLaRibhLPveP100.v; timeLaRibhLPveP100 = gotmkppLaRibhLPveP100.time
# GTMtempLaRibhLPveP200 = gotmkppLaRibhLPveP200.temp; GTMdepthLaRibhLPveP200 = gotmkppLaRibhLPveP200.z; GTMnumLaRibhLPveP200 = gotmkppLaRibhLPveP200.num; GTMnuhLaRibhLPveP200 = gotmkppLaRibhLPveP200.nuh; GTMuLaRibhLPveP200 = gotmkppLaRibhLPveP200.u; GTMvLaRibhLPveP200 = gotmkppLaRibhLPveP200.v; timeLaRibhLPveP200 = gotmkppLaRibhLPveP200.time
# GTMtempLaRibhLPveP500 = gotmkppLaRibhLPveP500.temp; GTMdepthLaRibhLPveP500 = gotmkppLaRibhLPveP500.z; GTMnumLaRibhLPveP500 = gotmkppLaRibhLPveP500.num; GTMnuhLaRibhLPveP500 = gotmkppLaRibhLPveP500.nuh; GTMuLaRibhLPveP500 = gotmkppLaRibhLPveP500.u; GTMvLaRibhLPveP500 = gotmkppLaRibhLPveP500.v; timeLaRibhLPveP500 = gotmkppLaRibhLPveP500.time
# GTMtempLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.temp; GTMdepthLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.z; GTMnumLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.num; GTMnuhLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.nuh; GTMuLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.u; GTMvLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.v; timeLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.time
# GTMtempLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.temp; GTMdepthLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.z; GTMnumLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.num; GTMnuhLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.nuh; GTMuLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.u; GTMvLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.v; timeLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.time
GTMtempLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.temp; GTMdepthLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.z; GTMnumLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.num; GTMnuhLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.nuh; GTMuLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.u; GTMvLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.v; timeLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.time
GTMtempLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.temp; GTMdepthLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.z; GTMnumLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.num; GTMnuhLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.nuh; GTMuLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.u; GTMvLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.v; timeLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.time
GTMtempLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.temp; GTMdepthLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.z; GTMnumLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.num; GTMnuhLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.nuh; GTMuLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.u; GTMvLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.v; timeLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.time
GTMtempLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.temp; GTMdepthLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.z; GTMnumLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.num; GTMnuhLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.nuh; GTMuLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.u; GTMvLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.v; timeLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.time
GTMtempLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.temp; GTMdepthLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.z; GTMnumLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.num; GTMnuhLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.nuh; GTMuLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.u; GTMvLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.v; timeLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.time
GTMtempLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.temp; GTMdepthLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.z; GTMnumLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.num; GTMnuhLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.nuh; GTMuLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.u; GTMvLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.v; timeLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.time
# GTMtempLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.temp; GTMdepthLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.z; GTMnumLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.num; GTMnuhLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.nuh; GTMuLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.u; GTMvLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.v; timeLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.time
# GTMtempLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.temp; GTMdepthLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.z; GTMnumLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.num; GTMnuhLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.nuh; GTMuLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.u; GTMvLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.v; timeLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.time


## Convert time
time2 = 480 + (time[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
time2L = 480 + (timeL[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
# time2La = 445.0 + (timeLa[:]/86400); 
time2LaRib = 480 + (timeLaRib[:]/86400); 
time2LaRibL = 480 + (timeLaRibL[:]/86400); 
# time2LaRibg0a1 = 267.75 + (timeLaRibg0a1[:]/86400); 
# time2LaRibg1a0 = 267.75 + (timeLaRibg1a0[:]/86400); 
# time2LaRibnum2 = 267.75 + (timeLaRibnum2[:]/86400); 
# time2LaRibnuh2 = 267.75 + (timeLaRibnuh2[:]/86400); 
# time2LaRibnum2nuh2 = 267.75 + (timeLaRibnum2nuh2[:]/86400); 

# time2LaRibnum15 = 267.75 + (timeLaRibnum15[:]/86400); 
# time2LaRibnuh15 = 267.75 + (timeLaRibnuh15[:]/86400); 
# time2LaRibnum15nuh15 = 267.75 + (timeLaRibnum15nuh15[:]/86400); 

# time2LaRibhL = 267.75 + (timeLaRibhL[:]/86400); 
# time2hLPve = 267.75 + (timehLPve[:]/86400); 
# time2hLPveN100 = 267.75 + (timehLPveN100[:]/86400); 
# time2hLPveP100 = 267.75 + (timehLPveP100[:]/86400); 
# time2LaRibhLPve = 267.75 + (timeLaRibhLPve[:]/86400); 
# # time2LaRibhLNve = 267.75 + (timeLaRibhLNve[:]/86400); 
# time2LaRibLaTPve03 = 267.75 + (timeLaRibLaTPve03[:]/86400); 
# time2LaRibhLPveN100 = 267.75 + (timeLaRibhLPveN100[:]/86400); 
# time2LaRibhLPveP100 = 267.75 + (timeLaRibhLPveP100[:]/86400); 
# time2LaRibhLPveP200 = 267.75 + (timeLaRibhLPveP200[:]/86400); 
# time2LaRibhLPveP500 = 267.75 + (timeLaRibhLPveP500[:]/86400); 
# time2LaRibhLPveP1000 = 267.75 + (timeLaRibhLPveP1000[:]/86400); 
# time2LaRibhLPveP2000 = 267.75 + (timeLaRibhLPveP2000[:]/86400); 
time2LaRibhLPveP3000 = 480 + (timeLaRibhLPveP3000[:]/86400); 
time2LaRibhLPveP3000L = 480 + (timeLaRibhLPveP3000L[:]/86400); 
time2LaRibhLPveP3500 = 480 + (timeLaRibhLPveP3500[:]/86400); 
time2LaRibhLPveP3500L = 480 + (timeLaRibhLPveP3500L[:]/86400); 
time2LaRibhLPveP4000 = 480 + (timeLaRibhLPveP4000[:]/86400); 
time2LaRibhLPveP4000L = 480 + (timeLaRibhLPveP4000L[:]/86400); 
# time2LaRibhLPveP5000 = 267.75 + (timeLaRibhLPveP5000[:]/86400); 
# time2LaRibhLPveP10000 = 267.75 + (timeLaRibhLPveP10000[:]/86400); 

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkppL = gotmkppL.assign(time2L=("time", time2L)); gotmkppL = gotmkppL.swap_dims({"time" : "time2L"})
# gotmkppLa = gotmkppLa.assign(time2La=("time", time2La)); gotmkppLa = gotmkppLa.swap_dims({"time" : "time2La"})
gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})
gotmkppLaRibL = gotmkppLaRibL.assign(time2LaRibL=("time", time2LaRibL)); gotmkppLaRibL = gotmkppLaRibL.swap_dims({"time" : "time2LaRibL"})
# gotmkppLaRibnum2 = gotmkppLaRibnum2.assign(time2LaRibnum2=("time", time2LaRibnum2)); gotmkppLaRibnum2 = gotmkppLaRibnum2.swap_dims({"time" : "time2LaRibnum2"})
# gotmkppLaRibnuh2 = gotmkppLaRibnuh2.assign(time2LaRibnuh2=("time", time2LaRibnuh2)); gotmkppLaRibnuh2 = gotmkppLaRibnuh2.swap_dims({"time" : "time2LaRibnuh2"})
# gotmkppLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.assign(time2LaRibnum2nuh2=("time", time2LaRibnum2nuh2)); gotmkppLaRibnum2nuh2 = gotmkppLaRibnum2nuh2.swap_dims({"time" : "time2LaRibnum2nuh2"})
# gotmkppLaRibhL = gotmkppLaRibhL.assign(time2LaRibhL=("time", time2LaRibhL)); gotmkppLaRibhL = gotmkppLaRibhL.swap_dims({"time" : "time2LaRibhL"})
# gotmkpphLPve = gotmkpphLPve.assign(time2hLPve=("time", time2hLPve)); gotmkpphLPve = gotmkpphLPve.swap_dims({"time" : "time2hLPve"})
# gotmkpphLPveP100 = gotmkpphLPveP100.assign(time2hLPveP100=("time", time2hLPveP100)); gotmkpphLPveP100 = gotmkpphLPveP100.swap_dims({"time" : "time2hLPveP100"})
# gotmkpphLPveN100 = gotmkpphLPveN100.assign(time2hLPveN100=("time", time2hLPveN100)); gotmkpphLPveN100 = gotmkpphLPveN100.swap_dims({"time" : "time2hLPveN100"})

# gotmkppLaRibg0a1 = gotmkppLaRibg0a1.assign(time2LaRibg0a1=("time", time2LaRibg0a1)); gotmkppLaRibg0a1 = gotmkppLaRibg0a1.swap_dims({"time" : "time2LaRibg0a1"})
# gotmkppLaRibg1a0 = gotmkppLaRibg1a0.assign(time2LaRibg1a0=("time", time2LaRibg1a0)); gotmkppLaRibg1a0 = gotmkppLaRibg1a0.swap_dims({"time" : "time2LaRibg1a0"})

# gotmkppLaRibnum15 = gotmkppLaRibnum15.assign(time2LaRibnum15=("time", time2LaRibnum15)); gotmkppLaRibnum15 = gotmkppLaRibnum15.swap_dims({"time" : "time2LaRibnum15"})
# gotmkppLaRibnuh15 = gotmkppLaRibnuh15.assign(time2LaRibnuh15=("time", time2LaRibnuh15)); gotmkppLaRibnuh15 = gotmkppLaRibnuh15.swap_dims({"time" : "time2LaRibnuh15"})
# gotmkppLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.assign(time2LaRibnum15nuh15=("time", time2LaRibnum15nuh15)); gotmkppLaRibnum15nuh15 = gotmkppLaRibnum15nuh15.swap_dims({"time" : "time2LaRibnum15nuh15"})

# gotmkppLaRibhLNve = gotmkppLaRibhLNve.assign(time2LaRibhLNve=("time", time2LaRibhLNve)); gotmkppLaRibhLNve = gotmkppLaRibhLNve.swap_dims({"time" : "time2LaRibhLNve"})
# gotmkppLaRibhLPve = gotmkppLaRibhLPve.assign(time2LaRibhLPve=("time", time2LaRibhLPve)); gotmkppLaRibhLPve = gotmkppLaRibhLPve.swap_dims({"time" : "time2LaRibhLPve"})
# gotmkppLaRibLaTPve03 = gotmkppLaRibLaTPve03.assign(time2LaRibLaTPve03=("time", time2LaRibLaTPve03)); gotmkppLaRibLaTPve03 = gotmkppLaRibLaTPve03.swap_dims({"time" : "time2LaRibLaTPve03"})
# gotmkppLaRibhLPveN100 = gotmkppLaRibhLPveN100.assign(time2LaRibhLPveN100=("time", time2LaRibhLPveN100)); gotmkppLaRibhLPveN100 = gotmkppLaRibhLPveN100.swap_dims({"time" : "time2LaRibhLPveN100"})
# gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.assign(time2LaRibhLPveP100=("time", time2LaRibhLPveP100)); gotmkppLaRibhLPveP100 = gotmkppLaRibhLPveP100.swap_dims({"time" : "time2LaRibhLPveP100"})
# gotmkppLaRibhLPveP200 = gotmkppLaRibhLPveP200.assign(time2LaRibhLPveP200=("time", time2LaRibhLPveP200)); gotmkppLaRibhLPveP200 = gotmkppLaRibhLPveP200.swap_dims({"time" : "time2LaRibhLPveP200"})
# gotmkppLaRibhLPveP500 = gotmkppLaRibhLPveP500.assign(time2LaRibhLPveP500=("time", time2LaRibhLPveP500)); gotmkppLaRibhLPveP500 = gotmkppLaRibhLPveP500.swap_dims({"time" : "time2LaRibhLPveP500"})
# gotmkppLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.assign(time2LaRibhLPveP1000=("time", time2LaRibhLPveP1000)); gotmkppLaRibhLPveP1000 = gotmkppLaRibhLPveP1000.swap_dims({"time" : "time2LaRibhLPveP1000"})
# gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.assign(time2LaRibhLPveP2000=("time", time2LaRibhLPveP2000)); gotmkppLaRibhLPveP2000 = gotmkppLaRibhLPveP2000.swap_dims({"time" : "time2LaRibhLPveP2000"})
gotmkppLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.assign(time2LaRibhLPveP3000=("time", time2LaRibhLPveP3000)); gotmkppLaRibhLPveP3000 = gotmkppLaRibhLPveP3000.swap_dims({"time" : "time2LaRibhLPveP3000"})
gotmkppLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.assign(time2LaRibhLPveP3000L=("time", time2LaRibhLPveP3000L)); gotmkppLaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.swap_dims({"time" : "time2LaRibhLPveP3000L"})
gotmkppLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.assign(time2LaRibhLPveP3500=("time", time2LaRibhLPveP3500)); gotmkppLaRibhLPveP3500 = gotmkppLaRibhLPveP3500.swap_dims({"time" : "time2LaRibhLPveP3500"})
gotmkppLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.assign(time2LaRibhLPveP3500L=("time", time2LaRibhLPveP3500L)); gotmkppLaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.swap_dims({"time" : "time2LaRibhLPveP3500L"})
gotmkppLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.assign(time2LaRibhLPveP4000=("time", time2LaRibhLPveP4000)); gotmkppLaRibhLPveP4000 = gotmkppLaRibhLPveP4000.swap_dims({"time" : "time2LaRibhLPveP4000"})
gotmkppLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.assign(time2LaRibhLPveP4000L=("time", time2LaRibhLPveP4000L)); gotmkppLaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.swap_dims({"time" : "time2LaRibhLPveP4000L"})
# gotmkppLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.assign(time2LaRibhLPveP5000=("time", time2LaRibhLPveP5000)); gotmkppLaRibhLPveP5000 = gotmkppLaRibhLPveP5000.swap_dims({"time" : "time2LaRibhLPveP5000"})
# gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.assign(time2LaRibhLPveP10000=("time", time2LaRibhLPveP10000)); gotmkppLaRibhLPveP10000 = gotmkppLaRibhLPveP10000.swap_dims({"time" : "time2LaRibhLPveP10000"})

# print(gotmkppLaRibhLPveP3000.variables.items())

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
### ----------------------------- calculate mixed layer depth ----------------------------- ###
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
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtemp-mld_tempKPP).argmin('z'); # print(z_indexes)
print('GOTM temp MLD indexes', z_indexes)
mldKPP = GTMtemp.z[z_indexes]; 
mldKPP = mldKPP.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model LOCAL
mld_tempKPP_Local = GTMtempL.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempL-mld_tempKPP_Local).argmin('z'); # print(z_indexes)
print('GOTM temp MLD indexes', z_indexes)
mldKPP_Local = GTMtempL.z[z_indexes]; 
mldKPP_Local = mldKPP_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# plt.plot(time2,mldKPP); plt.xlabel('time -julian day'); plt.ylabel('mixed layer depth -m')
# plt.axis([260,341,-140,0]); plt.title('GOTM Mixed Layer Depth')
# plt.savefig('GOTM_Mixed_Layer_Depth_test6.png'); plt.show()

# ## using KPP EV model
# mld_tempKPPLa = GTMtempLa.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLa-mld_tempKPPLa).argmin('z'); # print(z_indexes)
# mldKPPLa = GTMtempLa.z[z_indexes]; 
# mldKPPLa = mldKPPLa.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original
mld_tempKPPLaRib = GTMtempLaRib.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRib-mld_tempKPPLaRib).argmin('z'); # print(z_indexes)
mldKPPLaRib = GTMtempLaRib.z[z_indexes]; 
mldKPPLaRib = mldKPPLaRib.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original
mld_tempKPPLaRib_Local = GTMtempLaRibL.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibL-mld_tempKPPLaRib_Local).argmin('z'); # print(z_indexes)
mldKPPLaRib_Local = GTMtempLaRib.z[z_indexes]; 
mldKPPLaRib_Local = mldKPPLaRib_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original
# mld_tempKPPLaRibg0a1 = GTMtempLaRibg0a1.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg0a1-mld_tempKPPLaRibg0a1).argmin('z'); # print(z_indexes)
# mldKPPLaRibg0a1 = GTMtempLaRibg0a1.z[z_indexes]; 
# mldKPPLaRibg0a1 = mldKPPLaRibg0a1.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original
# mld_tempKPPLaRibg1a0 = GTMtempLaRibg1a0.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg1a0-mld_tempKPPLaRibg1a0).argmin('z'); # print(z_indexes)
# mldKPPLaRibg1a0 = GTMtempLaRibg1a0.z[z_indexes]; 
# mldKPPLaRibg1a0 = mldKPPLaRibg1a0.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib num*2
# mld_tempKPPLaRibnum2 = GTMtempLaRibnum2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibnum2-mld_tempKPPLaRibnum2).argmin('z'); # print(z_indexes)
# mldKPPLaRibnum2 = GTMtempLaRibnum2.z[z_indexes]; 
# mldKPPLaRibnum2 = mldKPPLaRibnum2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib nuh*2
# mld_tempKPPLaRibnuh2 = GTMtempLaRibnuh2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibnuh2-mld_tempKPPLaRibnuh2).argmin('z'); # print(z_indexes)
# mldKPPLaRibnuh2 = GTMtempLaRibnuh2.z[z_indexes]; 
# mldKPPLaRibnuh2 = mldKPPLaRibnuh2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original
# mld_tempKPPLaRibnum2nuh2 = GTMtempLaRibnum2nuh2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibnum2nuh2-mld_tempKPPLaRibnum2nuh2).argmin('z'); # print(z_indexes)
# mldKPPLaRibnum2nuh2 = GTMtempLaRibnum2nuh2.z[z_indexes]; 
# mldKPPLaRibnum2nuh2 = mldKPPLaRibnum2nuh2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib num*1.5
# mld_tempKPPLaRibnum15 = GTMtempLaRibnum15.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibnum15-mld_tempKPPLaRibnum15).argmin('z'); # print(z_indexes)
# mldKPPLaRibnum15 = GTMtempLaRibnum15.z[z_indexes]; 
# mldKPPLaRibnum15 = mldKPPLaRibnum15.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib nuh*1.5
# mld_tempKPPLaRibnuh15 = GTMtempLaRibnuh15.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibnuh15-mld_tempKPPLaRibnuh15).argmin('z'); # print(z_indexes)
# mldKPPLaRibnuh15 = GTMtempLaRibnuh15.z[z_indexes]; 
# mldKPPLaRibnuh15 = mldKPPLaRibnuh15.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original num*1.5 nuh*1.5
# mld_tempKPPLaRibnum15nuh15 = GTMtempLaRibnum15nuh15.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibnum15nuh15-mld_tempKPPLaRibnum15nuh15).argmin('z'); # print(z_indexes)
# mldKPPLaRibnum15nuh15 = GTMtempLaRibnum15nuh15.z[z_indexes]; 
# mldKPPLaRibnum15nuh15 = mldKPPLaRibnum15nuh15.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition
# mld_tempKPPLaRibhL = GTMtempLaRibhL.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhL-mld_tempKPPLaRibhL).argmin('z'); # print(z_indexes)
# mldKPPLaRibhL = GTMtempLaRibhL.z[z_indexes]; 
# mldKPPLaRibhL = mldKPPLaRibhL.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=0 for hL>0)
# mld_tempKPPLaRibhLPve = GTMtempLaRibhLPve.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPve-mld_tempKPPLaRibhLPve).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPve = GTMtempLaRibhLPve.z[z_indexes]; 
# mldKPPLaRibhLPve = mldKPPLaRibhLPve.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=0 for hL>0)
# mld_tempKPPhLPve = GTMtemphLPve.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemphLPve-mld_tempKPPhLPve).argmin('z'); # print(z_indexes)
# mldKPPhLPve = GTMtemphLPve.z[z_indexes]; 
# mldKPPhLPve = mldKPPhLPve.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=0 for hL>0) for gamma
# mld_tempKPPLaRibhLgamPve = GTMtempLaRibhLgamPve.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLgamPve-mld_tempKPPLaRibhLgamPve).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLgamPve = GTMtempLaRibhLgamPve.z[z_indexes]; 
# mldKPPLaRibhLgamPve = mldKPPLaRibhLgamPve.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=0 for hL<0)
# mld_tempKPPLaRibhLNve = GTMtempLaRibhLNve.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLNve-mld_tempKPPLaRibhLNve).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLNve = GTMtempLaRibhLNve.z[z_indexes]; 
# mldKPPLaRibhLNve = mldKPPLaRibhLNve.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for LaT>0.3)
# mld_tempKPPLaRibLaTPve03 = GTMtempLaRibLaTPve03.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibLaTPve03-mld_tempKPPLaRibLaTPve03).argmin('z'); # print(z_indexes)
# mldKPPLaRibLaTPve03 = GTMtempLaRibLaTPve03.z[z_indexes]; 
# mldKPPLaRibLaTPve03 = mldKPPLaRibLaTPve03.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>-100)
# mld_tempKPPLaRibhLPveN100 = GTMtempLaRibhLPveN100.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveN100-mld_tempKPPLaRibhLPveN100).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveN100 = GTMtempLaRibhLPveN100.z[z_indexes]; 
# mldKPPLaRibhLPveN100 = mldKPPLaRibhLPveN100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>-100)
# mld_tempKPPhLPveN100 = GTMtemphLPveN100.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemphLPveN100-mld_tempKPPhLPveN100).argmin('z'); # print(z_indexes)
# mldKPPhLPveN100 = GTMtemphLPveN100.z[z_indexes]; 
# mldKPPhLPveN100 = mldKPPhLPveN100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>-100) for gamma
# mld_tempKPPLaRibhLgamPveN100 = GTMtempLaRibhLgamPveN100.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLgamPveN100-mld_tempKPPLaRibhLgamPveN100).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLgamPveN100 = GTMtempLaRibhLgamPveN100.z[z_indexes]; 
# mldKPPLaRibhLgamPveN100 = mldKPPLaRibhLgamPveN100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100)
# mld_tempKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP100-mld_tempKPPLaRibhLPveP100).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP100 = GTMtempLaRibhLPveP100.z[z_indexes]; 
# mldKPPLaRibhLPveP100 = mldKPPLaRibhLPveP100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>200) for gamma
# mld_tempKPPLaRibhLPveP200 = GTMtempLaRibhLPveP200.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP200-mld_tempKPPLaRibhLPveP200).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP200 = GTMtempLaRibhLPveP200.z[z_indexes]; 
# mldKPPLaRibhLPveP200 = mldKPPLaRibhLPveP200.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>500)
# mld_tempKPPLaRibhLPveP500 = GTMtempLaRibhLPveP500.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP500-mld_tempKPPLaRibhLPveP500).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP500 = GTMtempLaRibhLPveP500.z[z_indexes]; 
# mldKPPLaRibhLPveP500 = mldKPPLaRibhLPveP500.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100)
# mld_tempKPPhLPveP100 = GTMtemphLPveP100.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemphLPveP100-mld_tempKPPhLPveP100).argmin('z'); # print(z_indexes)
# mldKPPhLPveP100 = GTMtemphLPveP100.z[z_indexes]; 
# mldKPPhLPveP100 = mldKPPhLPveP100.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>1000)
# mld_tempKPPLaRibhLPveP1000 = GTMtempLaRibhLPveP1000.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP1000-mld_tempKPPLaRibhLPveP1000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP1000 = GTMtempLaRibhLPveP1000.z[z_indexes]; 
# mldKPPLaRibhLPveP1000 = mldKPPLaRibhLPveP1000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP2000 = GTMtempLaRibhLPveP2000.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP2000-mld_tempKPPLaRibhLPveP2000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP2000 = GTMtempLaRibhLPveP2000.z[z_indexes]; 
# mldKPPLaRibhLPveP2000 = mldKPPLaRibhLPveP2000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP3000 = GTMtempLaRibhLPveP3000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP3000-mld_tempKPPLaRibhLPveP3000).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP3000 = GTMtempLaRibhLPveP3000.z[z_indexes]; 
mldKPPLaRibhLPveP3000 = mldKPPLaRibhLPveP3000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)
print('MLD for h/L > 3000', mldKPPLaRibhLPveP3000)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP3000_Local = GTMtempLaRibhLPveP3000L.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP3000L-mld_tempKPPLaRibhLPveP3000_Local).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP3000_Local = GTMtempLaRibhLPveP3000L.z[z_indexes]; 
mldKPPLaRibhLPveP3000_Local = mldKPPLaRibhLPveP3000_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP3500 = GTMtempLaRibhLPveP3500.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP3500-mld_tempKPPLaRibhLPveP3500).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP3500 = GTMtempLaRibhLPveP3500.z[z_indexes]; 
mldKPPLaRibhLPveP3500 = mldKPPLaRibhLPveP3500.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP3500_Local = GTMtempLaRibhLPveP3500L.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP3500L-mld_tempKPPLaRibhLPveP3500_Local).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP3500_Local = GTMtempLaRibhLPveP3500L.z[z_indexes]; 
mldKPPLaRibhLPveP3500_Local = mldKPPLaRibhLPveP3500_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP4000 = GTMtempLaRibhLPveP4000.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP4000-mld_tempKPPLaRibhLPveP4000).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP4000 = GTMtempLaRibhLPveP4000.z[z_indexes]; 
mldKPPLaRibhLPveP4000 = mldKPPLaRibhLPveP4000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
mld_tempKPPLaRibhLPveP4000_Local = GTMtempLaRibhLPveP4000L.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLaRibhLPveP4000L-mld_tempKPPLaRibhLPveP4000_Local).argmin('z'); # print(z_indexes)
mldKPPLaRibhLPveP4000_Local = GTMtempLaRibhLPveP4000L.z[z_indexes]; 
mldKPPLaRibhLPveP4000_Local = mldKPPLaRibhLPveP4000_Local.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>500)
# mld_tempKPPLaRibhLPveP5000 = GTMtempLaRibhLPveP5000.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP5000-mld_tempKPPLaRibhLPveP5000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP5000 = GTMtempLaRibhLPveP5000.z[z_indexes]; 
# mldKPPLaRibhLPveP5000 = mldKPPLaRibhLPveP5000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original with hL condition (Ch=1 for hL>100) for gamma
# mld_tempKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibhLPveP10000-mld_tempKPPLaRibhLPveP10000).argmin('z'); # print(z_indexes)
# mldKPPLaRibhLPveP10000 = GTMtempLaRibhLPveP10000.z[z_indexes]; 
# mldKPPLaRibhLPveP10000 = mldKPPLaRibhLPveP10000.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

### ----------------------------- moving averages h/L ----------------------------- ###

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

GTM_hL_LaRib = gotmkppLaRib.hL; GTM_hLR_LaRib = gotmkppLaRib.hLR; GTM_Lat_LaRib = gotmkppLaRib.La_t; GTM_tdotus_LaRib = gotmkppLaRib.tdotus; 
GTM_hL_LaRibL = gotmkppLaRibL.hL; GTM_hLR_LaRibL = gotmkppLaRibL.hLR; GTM_Lat_LaRibL = gotmkppLaRibL.La_t; GTM_tdotus_LaRibL = gotmkppLaRibL.tdotus; 
# GTM_hL_LaRibLaTPve03 = gotmkppLaRibLaTPve03.hL; GTM_hLR_LaRibLaTPve03 = gotmkppLaRibLaTPve03.hLR; GTM_Lat_LaRibLaTPve03 = gotmkppLaRibLaTPve03.La_t; GTM_tdotus_LaRibLaTPve03 = gotmkppLaRibLaTPve03.tdotus; 
# GTM_hL_LaRibhLPveLaTPve03 = gotmkppLaRibhLPveLaTPve03.hL; GTM_hLR_LaRibhLPveLaTPve03 = gotmkppLaRibhLPveLaTPve03.hLR; GTM_Lat_LaRibhLPveLaTPve03 = gotmkppLaRibhLPveLaTPve03.La_t; GTM_tdotus_LaRibhLPveLaTPve03 = gotmkppLaRibhLPveLaTPve03.tdotus; 
# GTM_hL_LaRibhLPve = gotmkppLaRibhLPve.hL; GTM_hLR_LaRibhLPve = gotmkppLaRibhLPve.hLR; GTM_Lat_LaRibhLPve = gotmkppLaRibhLPve.La_t; GTM_tdotus_LaRibhLPve = gotmkppLaRibhLPve.tdotus; 
# GTM_hL_hLPve = gotmkpphLPve.hL; GTM_hLR_hLPve = gotmkpphLPve.hLR; GTM_Lat_hLPve = gotmkpphLPve.La_t; GTM_tdotus_hLPve = gotmkpphLPve.tdotus; 
# GTM_hL_hLPveP100 = gotmkpphLPveP100.hL; GTM_hLR_hLPveP100 = gotmkpphLPveP100.hLR; GTM_Lat_hLPveP100 = gotmkpphLPveP100.La_t; GTM_tdotus_hLPveP100 = gotmkppLaRibhLPveP100.tdotus; 
# GTM_hL_hLPveN100 = gotmkpphLPveN100.hL; GTM_hLR_hLPveN100 = gotmkpphLPveN100.hLR; GTM_Lat_hLPveN100 = gotmkpphLPveN100.La_t; GTM_tdotus_hLPveN100 = gotmkppLaRibhLPveN100.tdotus; 
# GTM_hL_LaRibg0a1 = gotmkppLaRibg0a1.hL; GTM_hLR_LaRibg0a1 = gotmkppLaRibg0a1.hLR; GTM_Lat_LaRibg0a1 = gotmkppLaRibg0a1.La_t; GTM_tdotus_LaRibg0a1 = gotmkppLaRibg0a1.tdotus; 
# GTM_hL_LaRibg1a0 = gotmkppLaRibg1a0.hL; GTM_hLR_LaRibg1a0 = gotmkppLaRibg1a0.hLR; GTM_Lat_LaRibg1a0 = gotmkppLaRibg1a0.La_t; GTM_tdotus_LaRibg1a0 = gotmkppLaRibg1a0.tdotus; 

# GTM_hL_LaRibhLgamPve = gotmkppLaRibhLgamPve.hL; GTM_hLR_LaRibhLgamPve = gotmkppLaRibhLgamPve.hLR; GTM_Lat_LaRibhLgamPve = gotmkppLaRibhLgamPve.La_t; GTM_tdotus_LaRibhLgamPve = gotmkppLaRibhLgamPve.tdotus; 
# GTM_hL_LaRibhLgamPveP100 = gotmkppLaRibhLgamPveP100.hL; GTM_hLR_LaRibhLgamPveP100 = gotmkppLaRibhLgamPveP100.hLR; GTM_Lat_LaRibhLgamPveP100 = gotmkppLaRibhLgamPveP100.La_t; GTM_tdotus_LaRibhLgamPveP100 = gotmkppLaRibhLgamPveP100.tdotus; 
# GTM_hL_LaRibhLgamPveN100 = gotmkppLaRibhLgamPveN100.hL; GTM_hLR_LaRibhLgamPveN100 = gotmkppLaRibhLgamPveN100.hLR; GTM_Lat_LaRibhLgamPveN100 = gotmkppLaRibhLgamPveN100.La_t; GTM_tdotus_LaRibhLgamPveN100 = gotmkppLaRibhLgamPveN100.tdotus; 
# GTM_hL_LaRibhLNve = gotmkppLaRibhLNve.hL; GTM_hLR_LaRibhLNve = gotmkppLaRibhLNve.hLR; GTM_Lat_LaRibhLNve = gotmkppLaRibhLNve.La_t; GTM_tdotus_LaRibhLNve = gotmkppLaRibhLNve.tdotus; 
# GTM_hL_LaRibhLPveP100 = gotmkppLaRibhLPveP100.hL; GTM_hLR_LaRibhLPveP100 = gotmkppLaRibhLPveP100.hLR; GTM_Lat_LaRibhLPveP100 = gotmkppLaRibhLPveP100.La_t; GTM_tdotus_LaRibhLPveP100 = gotmkppLaRibhLPveP100.tdotus; 
# GTM_hL_LaRibhLPveN100 = gotmkppLaRibhLPveN100.hL; GTM_hLR_LaRibhLPveN100 = gotmkppLaRibhLPveN100.hLR; GTM_Lat_LaRibhLPveN100 = gotmkppLaRibhLPveN100.La_t; GTM_tdotus_LaRibhLPveN100 = gotmkppLaRibhLPveN100.tdotus; 
# GTM_hL_LaRibhLPveN200 = gotmkppLaRibhLPveN200.hL; GTM_hLR_LaRibhLPveN200 = gotmkppLaRibhLPveN200.hLR; GTM_Lat_LaRibhLPveN200 = gotmkppLaRibhLPveN200.La_t; GTM_tdotus_LaRibhLPveN200 = gotmkppLaRibhLPveN200.tdotus; 
# GTM_hL_LaRibhLPveN500 = gotmkppLaRibhLPveN500.hL; GTM_hLR_LaRibhLPveN500 = gotmkppLaRibhLPveN500.hLR; GTM_Lat_LaRibhLPveN500 = gotmkppLaRibhLPveN500.La_t; GTM_tdotus_LaRibhLPveN500 = gotmkppLaRibhLPveN500.tdotus; 
# GTM_hL_LaRibhLPveP500 = gotmkppLaRibhLPveP500.hL; GTM_hLR_LaRibhLPveP500 = gotmkppLaRibhLPveP500.hLR; GTM_Lat_LaRibhLPveP500 = gotmkppLaRibhLPveP500.La_t; GTM_tdotus_LaRibhLPveP500 = gotmkppLaRibhLPveP500.tdotus; 
# GTM_hL_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.hL; GTM_hLR_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.hLR; GTM_Lat_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.La_t; GTM_tdotus_LaRibhLPveN1000 = gotmkppLaRibhLPveN1000.tdotus; 
# GTM_hL_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.hL; GTM_hLR_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.hLR; GTM_Lat_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.La_t; GTM_tdotus_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.tdotus; GTM_Theatflux_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.total; GTM_fric_LaRibhLPveP1000 = gotmkppLaRibhLPveP1000.u_taus;
# GTM_hL_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.hL; GTM_hLR_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.hLR; GTM_Lat_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.La_t; GTM_tdotus_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.tdotus; GTM_Theatflux_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.total; GTM_fric_LaRibhLPveP2000 = gotmkppLaRibhLPveP2000.u_taus; 
GTM_hL_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.hL; GTM_hLR_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.hLR; GTM_Lat_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.La_t; GTM_tdotus_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.tdotus; GTM_Theatflux_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.total; GTM_fric_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.u_taus; GTM_us_LaRibhLPveP3000 = gotmkppLaRibhLPveP3000.u_s; 
GTM_hL_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.hL; GTM_hLR_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.hLR; GTM_Lat_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.La_t; GTM_tdotus_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.tdotus; GTM_Theatflux_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.total; GTM_fric_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.u_taus; GTM_us_LaRibhLPveP3000L = gotmkppLaRibhLPveP3000L.u_s; 
GTM_hL_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.hL; GTM_hLR_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.hLR; GTM_Lat_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.La_t; GTM_tdotus_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.tdotus; GTM_Theatflux_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.total; GTM_fric_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.u_taus; GTM_us_LaRibhLPveP3500 = gotmkppLaRibhLPveP3500.u_s; 
GTM_hL_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.hL; GTM_hLR_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.hLR; GTM_Lat_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.La_t; GTM_tdotus_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.tdotus; GTM_Theatflux_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.total; GTM_fric_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.u_taus; GTM_us_LaRibhLPveP3500L = gotmkppLaRibhLPveP3500L.u_s; 
GTM_hL_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.hL; GTM_hLR_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.hLR; GTM_Lat_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.La_t; GTM_tdotus_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.tdotus; GTM_Theatflux_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.total; GTM_fric_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.u_taus; GTM_us_LaRibhLPveP4000 = gotmkppLaRibhLPveP4000.u_s; 
GTM_hL_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.hL; GTM_hLR_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.hLR; GTM_Lat_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.La_t; GTM_tdotus_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.tdotus; GTM_Theatflux_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.total; GTM_fric_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.u_taus; GTM_us_LaRibhLPveP4000L = gotmkppLaRibhLPveP4000L.u_s; 
# GTM_hL_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.hL; GTM_hLR_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.hLR; GTM_Lat_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.La_t; GTM_tdotus_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.tdotus; GTM_Theatflux_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.total; GTM_fric_LaRibhLPveP5000 = gotmkppLaRibhLPveP5000.u_taus; 


# extract data for hL analysis
# GTM_hLR_LaRibhL2 = GTM_hLR_LaRibhL2.isel(lon=0,lat=0)
GTM_hLR_LaRib = GTM_hLR_LaRib.isel(lon=0,lat=0)
GTM_hLR_LaRibL = GTM_hLR_LaRibL.isel(lon=0,lat=0)
# GTM_hLR_LaRibg0a1 = GTM_hLR_LaRibg0a1.isel(lon=0,lat=0)
# GTM_hLR_LaRibg1a0 = GTM_hLR_LaRibg1a0.isel(lon=0,lat=0)
# GTM_hLR_LaRibLaTPve03 = GTM_hLR_LaRibLaTPve03.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveLaTPve03 = GTM_hLR_LaRibhLPveLaTPve03.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPve = GTM_hLR_LaRibhLPve.isel(lon=0,lat=0)
# GTM_hLR_hLPve = GTM_hLR_hLPve.isel(lon=0,lat=0)
# GTM_hLR_hLPveP100 = GTM_hLR_hLPveP100.isel(lon=0,lat=0)
# GTM_hLR_hLPveN100 = GTM_hLR_hLPveN100.isel(lon=0,lat=0)

# GTM_hLR_LaRibhLgamPve = GTM_hLR_LaRibhLgamPve.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLgamPveP100 = GTM_hLR_LaRibhLgamPveP100.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLgamPveN100 = GTM_hLR_LaRibhLgamPveN100.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLNve = GTM_hLR_LaRibhLNve.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP100 = GTM_hLR_LaRibhLPveP100.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveN100 = GTM_hLR_LaRibhLPveN100.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveN200 = GTM_hLR_LaRibhLPveN200.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveN500 = GTM_hLR_LaRibhLPveN500.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP500 = GTM_hLR_LaRibhLPveP500.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveN1000 = GTM_hLR_LaRibhLPveN1000.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP1000 = GTM_hLR_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP2000 = GTM_hLR_LaRibhLPveP2000.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP3000 = GTM_hLR_LaRibhLPveP3000.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP3000L = GTM_hLR_LaRibhLPveP3000L.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP3500 = GTM_hLR_LaRibhLPveP3500.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP3500L = GTM_hLR_LaRibhLPveP3500L.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP4000 = GTM_hLR_LaRibhLPveP4000L.isel(lon=0,lat=0)
GTM_hLR_LaRibhLPveP4000L = GTM_hLR_LaRibhLPveP4000L.isel(lon=0,lat=0)
# GTM_hLR_LaRibhLPveP5000 = GTM_hLR_LaRibhLPveP5000.isel(lon=0,lat=0)

# La_t - extract data for hL analysis
# GTM_Lat_LaRibhLPveP1000 = GTM_Lat_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP2000 = GTM_Lat_LaRibhLPveP2000.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP3000 = GTM_Lat_LaRibhLPveP3000.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP3000L = GTM_Lat_LaRibhLPveP3000L.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP3500 = GTM_Lat_LaRibhLPveP3500.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP3500L = GTM_Lat_LaRibhLPveP3500L.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP4000 = GTM_Lat_LaRibhLPveP4000.isel(lon=0,lat=0)
GTM_Lat_LaRibhLPveP4000L = GTM_Lat_LaRibhLPveP4000L.isel(lon=0,lat=0)
# GTM_Lat_LaRibhLPveP5000 = GTM_Lat_LaRibhLPveP5000.isel(lon=0,lat=0)

# Total heat flux with radiation - extract data for hL analysis
# GTM_Theatflux_LaRibhLPveP1000 = GTM_Theatflux_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP2000 = GTM_Theatflux_LaRibhLPveP2000.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP3000 = GTM_Theatflux_LaRibhLPveP3000.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP3000L = GTM_Theatflux_LaRibhLPveP3000L.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP3500 = GTM_Theatflux_LaRibhLPveP3500.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP3500L = GTM_Theatflux_LaRibhLPveP3500L.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP4000 = GTM_Theatflux_LaRibhLPveP4000.isel(lon=0,lat=0)
GTM_Theatflux_LaRibhLPveP4000L = GTM_Theatflux_LaRibhLPveP4000L.isel(lon=0,lat=0)
# GTM_Theatflux_LaRibhLPveP5000 = GTM_Theatflux_LaRibhLPveP5000.isel(lon=0,lat=0)

# windStress to Stokes drfit flux density - extract data for hL analysis
# GTM_tdotus_LaRibhLPveP1000 = GTM_tdotus_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP2000 = GTM_tdotus_LaRibhLPveP2000.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP3000 = GTM_tdotus_LaRibhLPveP3000.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP3000L = GTM_tdotus_LaRibhLPveP3000L.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP3500 = GTM_tdotus_LaRibhLPveP3500.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP3500L = GTM_tdotus_LaRibhLPveP3500L.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP4000 = GTM_tdotus_LaRibhLPveP4000.isel(lon=0,lat=0)
GTM_tdotus_LaRibhLPveP4000L = GTM_tdotus_LaRibhLPveP4000L.isel(lon=0,lat=0)
# GTM_tdotus_LaRibhLPveP5000 = GTM_tdotus_LaRibhLPveP5000.isel(lon=0,lat=0)

# friction velocity - extract data for hL analysis
# GTM_fric_LaRibhLPveP1000 = GTM_fric_LaRibhLPveP1000.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP2000 = GTM_fric_LaRibhLPveP2000.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP3000 = GTM_fric_LaRibhLPveP3000.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP3000L = GTM_fric_LaRibhLPveP3000L.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP3500 = GTM_fric_LaRibhLPveP3500.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP3500L = GTM_fric_LaRibhLPveP3500L.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP4000 = GTM_fric_LaRibhLPveP4000.isel(lon=0,lat=0)
GTM_fric_LaRibhLPveP4000L = GTM_fric_LaRibhLPveP4000L.isel(lon=0,lat=0)
# GTM_fric_LaRibhLPveP5000 = GTM_fric_LaRibhLPveP5000.isel(lon=0,lat=0)

# Stokes drift velocity - extract data for hL analysis
GTM_us_LaRibhLPveP3000 = GTM_us_LaRibhLPveP3000.isel(lon=0,lat=0)
GTM_us_LaRibhLPveP3000L = GTM_us_LaRibhLPveP3000L.isel(lon=0,lat=0)
GTM_us_LaRibhLPveP3500 = GTM_us_LaRibhLPveP3500.isel(lon=0,lat=0)
GTM_us_LaRibhLPveP3500L = GTM_us_LaRibhLPveP3500L.isel(lon=0,lat=0)
GTM_us_LaRibhLPveP4000 = GTM_us_LaRibhLPveP4000.isel(lon=0,lat=0)
GTM_us_LaRibhLPveP4000L = GTM_us_LaRibhLPveP4000L.isel(lon=0,lat=0)



# # apply moving average for 1 week
# time2LaRibhL2_1015 = moving_average(time2LaRibhL2,1015)
# GTM_hL_LaRibhL2_1015 = moving_average(GTM_hL_LaRibhL2.isel(lon=0,lat=0),1015)
# GTM_hLR_LaRibhL2_1015 = moving_average(GTM_hLR_LaRibhL2,1015)

#---------------------------------------------#
#          hL parameter -- SMA=2wks           #
#---------------------------------------------#
# hL2_test
# time2LaRibhL2_2030 = moving_average(time2LaRibhL2,2030)
# GTM_hL_LaRibhL2_2030 = moving_average(GTM_hL_LaRibhL2.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhL2_2030 = moving_average(GTM_hLR_LaRibhL2,2030)
# nu_LaT_Rib Ch=Cm=1
# time2LaRib_2030 = moving_average(time2LaRib,2030)
# GTM_hL_LaRib_2030 = moving_average(GTM_hL_LaRib.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRib_2030 = moving_average(GTM_hLR_LaRib,2030)
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
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>500
# time2LaRibhLPveP500_2030 = moving_average(time2LaRibhLPveP500,2030)
# GTM_hL_LaRibhLPveP500_2030 = moving_average(GTM_hL_LaRibhLPveP500.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP500_2030 = moving_average(GTM_hLR_LaRibhLPveP500,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>-1000
# time2LaRibhLPveN1000_2030 = moving_average(time2LaRibhLPveN1000,2030)
# GTM_hL_LaRibhLPveN1000_2030 = moving_average(GTM_hL_LaRibhLPveN1000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveN1000_2030 = moving_average(GTM_hLR_LaRibhLPveN1000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>1000
# time2LaRibhLPveP1000_2030 = moving_average(time2LaRibhLPveP1000,2030)
# GTM_hL_LaRibhLPveP1000_2030 = moving_average(GTM_hL_LaRibhLPveP1000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP1000_2030 = moving_average(GTM_hLR_LaRibhLPveP1000,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>2000
# time2LaRibhLPveP2000_2030 = moving_average(time2LaRibhLPveP2000,2030)
# GTM_hL_LaRibhLPveP2000_2030 = moving_average(GTM_hL_LaRibhLPveP2000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP2000_2030 = moving_average(GTM_hLR_LaRibhLPveP2000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
time2LaRibhLPveP3000_2030 = moving_average(time2LaRibhLPveP3000,2030)
GTM_hL_LaRibhLPveP3000_2030 = moving_average(GTM_hL_LaRibhLPveP3000.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP3000_2030 = moving_average(GTM_hLR_LaRibhLPveP3000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000 LOCAL
time2LaRibhLPveP3000L_2030 = moving_average(time2LaRibhLPveP3000L,2030)
GTM_hL_LaRibhLPveP3000L_2030 = moving_average(GTM_hL_LaRibhLPveP3000L.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP3000L_2030 = moving_average(GTM_hLR_LaRibhLPveP3000L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
time2LaRibhLPveP3500_2030 = moving_average(time2LaRibhLPveP3500,2030)
GTM_hL_LaRibhLPveP3500_2030 = moving_average(GTM_hL_LaRibhLPveP3500.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP3500_2030 = moving_average(GTM_hLR_LaRibhLPveP3500,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500 LOCAL
time2LaRibhLPveP3500L_2030 = moving_average(time2LaRibhLPveP3500L,2030)
GTM_hL_LaRibhLPveP3500L_2030 = moving_average(GTM_hL_LaRibhLPveP3500L.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP3500L_2030 = moving_average(GTM_hLR_LaRibhLPveP3500L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
time2LaRibhLPveP4000_2030 = moving_average(time2LaRibhLPveP4000,2030)
GTM_hL_LaRibhLPveP4000_2030 = moving_average(GTM_hL_LaRibhLPveP4000.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP4000_2030 = moving_average(GTM_hLR_LaRibhLPveP4000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000 LOCAL
time2LaRibhLPveP4000L_2030 = moving_average(time2LaRibhLPveP4000L,2030)
GTM_hL_LaRibhLPveP4000L_2030 = moving_average(GTM_hL_LaRibhLPveP4000L.isel(lon=0,lat=0),2030)
GTM_hLR_LaRibhLPveP4000L_2030 = moving_average(GTM_hLR_LaRibhLPveP4000L,2030)
# # nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>5000
# time2LaRibhLPveP5000_2030 = moving_average(time2LaRibhLPveP5000,2030)
# GTM_hL_LaRibhLPveP5000_2030 = moving_average(GTM_hL_LaRibhLPveP5000.isel(lon=0,lat=0),2030)
# GTM_hLR_LaRibhLPveP5000_2030 = moving_average(GTM_hLR_LaRibhLPveP5000,2030)

#---------------------------------------------#
#          LaT parameter -- SMA=2wks          #
#---------------------------------------------#
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
GTM_Lat_LaRibhLPveP3000_2030 = moving_average(GTM_Lat_LaRibhLPveP3000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000 LOCAL
GTM_Lat_LaRibhLPveP3000L_2030 = moving_average(GTM_Lat_LaRibhLPveP3000L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
GTM_Lat_LaRibhLPveP3500_2030 = moving_average(GTM_Lat_LaRibhLPveP3500,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500 LOCAL
GTM_Lat_LaRibhLPveP3500L_2030 = moving_average(GTM_Lat_LaRibhLPveP3500L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
GTM_Lat_LaRibhLPveP4000_2030 = moving_average(GTM_Lat_LaRibhLPveP4000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000 LOCAL
GTM_Lat_LaRibhLPveP4000L_2030 = moving_average(GTM_Lat_LaRibhLPveP4000L,2030)

#--------------------------------------------------------#
#          Heat flux with radiation -- SMA=2wks          #
#--------------------------------------------------------#
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
GTM_Theatflux_LaRibhLPveP3000_2030 = moving_average(GTM_Theatflux_LaRibhLPveP3000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000 LOCAL
GTM_Theatflux_LaRibhLPveP3000L_2030 = moving_average(GTM_Theatflux_LaRibhLPveP3000L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
GTM_Theatflux_LaRibhLPveP3500_2030 = moving_average(GTM_Theatflux_LaRibhLPveP3500,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500 LOCAL
GTM_Theatflux_LaRibhLPveP3500L_2030 = moving_average(GTM_Theatflux_LaRibhLPveP3500L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
GTM_Theatflux_LaRibhLPveP4000_2030 = moving_average(GTM_Theatflux_LaRibhLPveP4000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000 LOCAL
GTM_Theatflux_LaRibhLPveP4000L_2030 = moving_average(GTM_Theatflux_LaRibhLPveP4000L,2030)

#-------------------------------------------------#
#          friction velocity -- SMA=2wks          #
#-------------------------------------------------#
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
GTM_fric_LaRibhLPveP3000_2030 = moving_average(GTM_fric_LaRibhLPveP3000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000 LOCAL
GTM_fric_LaRibhLPveP3000L_2030 = moving_average(GTM_fric_LaRibhLPveP3000L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
GTM_fric_LaRibhLPveP3500_2030 = moving_average(GTM_fric_LaRibhLPveP3500,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500 LOCAL
GTM_fric_LaRibhLPveP3500L_2030 = moving_average(GTM_fric_LaRibhLPveP3500L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
GTM_fric_LaRibhLPveP4000_2030 = moving_average(GTM_fric_LaRibhLPveP4000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000 LOCAL
GTM_fric_LaRibhLPveP4000L_2030 = moving_average(GTM_fric_LaRibhLPveP4000L,2030)

#-----------------------------------------------------#
#          Stokes drift velocity -- SMA=2wks          #
#-----------------------------------------------------#
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000
GTM_us_LaRibhLPveP3000_2030 = moving_average(GTM_us_LaRibhLPveP3000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3000 LOCAL
GTM_us_LaRibhLPveP3000L_2030 = moving_average(GTM_us_LaRibhLPveP3000L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500
GTM_us_LaRibhLPveP3500_2030 = moving_average(GTM_us_LaRibhLPveP3500,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>3500 LOCAL
GTM_us_LaRibhLPveP3500L_2030 = moving_average(GTM_us_LaRibhLPveP3500L,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000
GTM_us_LaRibhLPveP4000_2030 = moving_average(GTM_us_LaRibhLPveP4000,2030)
# nu_LaT_Rib Ch=Cm=1 w Ch=0 for hL>4000 LOCAL
GTM_us_LaRibhLPveP4000L_2030 = moving_average(GTM_us_LaRibhLPveP4000L,2030)


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

print(j)

OSMOSIS_MLD_time = np.c_[j,i]

print('OSMOSIS MLD with time',OSMOSIS_MLD_time)

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


print('manualy moving average :', df_mld)
# df_mld = pd.read_csv('mld.csv')
# print(df_mld.info())
# print(df_mld)
# df_mld = pd.DataFrame(i)
# df_mld = i.rolling(3,min_periods=1).mean()
# print('moving average MLD', df_mld)


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

# ########################################################################
# ## Comparison of MLDs using hL thresholds ranging [100,200,500,1000]  ##
# ########################################################################


# fig, ax = plt.subplots(4,1)

# ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[0].plot(time2LaRibhLPveP100,mldKPPLaRibhLPveP100,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ hL>100$')
# ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[1].plot(time2LaRibhLPveP200,mldKPPLaRibhLPveP200,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ hL>200$')
# ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
# ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

# ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[2].plot(time2LaRibhLPveP500,mldKPPLaRibhLPveP500,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ hL>500$')
# ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
# ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

# ax[3].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[3].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[3].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[3].plot(time2LaRibhLPveP1000,mldKPPLaRibhLPveP1000,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ hL>1000$')
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
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [100,200,500,1000]')
# plt.show()

# ############################################################################
# ## Comparison of MLDs using hL thresholds ranging [2000,3000,3500,4000]  ##
# ############################################################################


# fig, ax = plt.subplots(4,1)

# ax[0].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[0].plot(time2LaRibhLPveP1000,mldKPPLaRibhLPveP2000,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>2000$')
# ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
# # ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

# ax[1].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[1].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[1].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[1].plot(time2LaRibhLPveP2000,mldKPPLaRibhLPveP3000,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
# ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
# # ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

# ax[2].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[2].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[2].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[2].plot(time2LaRibhLPveP3500,mldKPPLaRibhLPveP3500,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
# ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
# # ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

# ax[3].plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5')
# ax[3].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
# ax[3].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1')
# ax[3].plot(time2LaRibhLPveP4000,mldKPPLaRibhLPveP4000,color='darkorange',linestyle='dashed',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[3].grid(); ax[3].legend(loc='best',fontsize='small')
# # ax[3].set_xlim(490,635); ax[3].set_ylim(-150,0)

# fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
# ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
# ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
# ax[3].text(0.06, 0.08, '(d)', horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)

# # plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# # plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [2000,3000,3500,4000]')
# plt.show()

#################################################################################
## Evolution of MLDs, hL and Lat using hL thresholds ranging [3000,3500,4000]  ##
#################################################################################

fig, ax = plt.subplots(2,1)

ax[0].plot(df_dats5,np.negative(df_mld5), color='darkgoldenrod',linestyle='solid',label = 'OSMOSIS SMA n=5')
ax[0].plot(time2L,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_ylabel('Mixed layer depth -m');
# ax[0].set_xlim(480,650)
ax[0].set_xlim(470,635); ax[0].set_ylim(-500,50)

ax[1].plot(df_dats5,np.negative(df_mld5), color='darkgoldenrod',linestyle='solid',label = 'OSMOSIS SMA n=5')
ax[1].plot(time2LaRibhLPveP3500,mldKPPLaRibhLPveP3500,color='maroon',linestyle='solid',label = r'KPP original $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
# ax[1].set_ylabel(r'$h/L$'); ax[1].set_yscale('symlog',nonposy='clip', linthreshy=0.01)
ax[1].set_xlim(470,635); ax[1].set_ylim(-200,50)

# ax[2].plot(time2LaRibhLPveP4000L,GTM_Lat_LaRibhLPveP4000L,color='black',linestyle='dashed',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[2].plot(time2LaRibhLPveP4000L_2030,GTM_Lat_LaRibhLPveP4000L_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
# ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
# ax[2].set_ylabel(r'$La_T$');
# ax[2].set_ylim(0,1)
# # ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.02, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')



######################################################
##               Temperature profiles               ##
######################################################
GTM_sel = [0];
Gldr_sel = [2543]


for j in Gldr_sel:
	GldrTemp = GldrTemp.isel()[j,:]
	Gldrdepth = GldrDepth[:]; Gldrtime = new_time[j]
	for i in GTM_sel:
		#
		kpptemp1d = GTMtemp.isel(lon=0,lat=0)[i,:]; 
		kppdepth = GTMdepth[:]
		kpptime = time2[i]
		#
		kpptempLaRibhLPveP35001d = GTMtempLaRibhLPveP3500.isel(lon=0,lat=0)[i,:]; 
		kppdepthLaRibhLPveP3500 = GTMdepthLaRibhLPveP3500[:]
		kpptimeLaRibhLPveP3500 = time2LaRibhLPveP3500[i]

fig, ax = plt.subplots(1,2)

ax[0].plot(kpptemp1d,kppdepth,'seagreen',label = r'GOTM KPP %0.2f'%kpptime.values);
ax[0].plot(GldrTemp,Gldrdepth,'gold',label = r'OSMOSIS %0.2f'%Gldrtime.values);
ax[0].grid(); ax[0].legend(loc='best', fontsize='small')
# ax[0].set_ylabel('Depth -m');ax[0].set_xlabel(r'Temperature - $\degree C$'); 
ax[0].set_ylim(-200,0);

ax[1].plot(kpptempLaRibhLPveP35001d,kppdepthLaRibhLPveP3500,'seagreen',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ h/L>3500$ %0.2f'%kpptimeLaRibhLPveP3500.values);
ax[1].plot(GldrTemp,Gldrdepth,'gold',label = r'OSMOSIS %0.2f'%Gldrtime.values);
ax[1].grid(); ax[1].legend(loc='best', fontsize='small')
# ax[1].set_ylabel('Depth -m');ax[1].set_xlabel(r'Temperature - $\degree C$'); 
ax[1].set_ylim(-200,0);

fig.text(0.5, 0.04, r'Temperature - $\degree C$', ha='center', va='center')
fig.text(0.08, 0.5, r'Depth - $m$', ha='center', va='center', rotation='vertical')
ax[0].text(0.02, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.02, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
# ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
plt.show()

# plt.title('KPP tracer variable profiles for %s'%kpptime.values)
# plt.axis([-0.01,0.15,-100,0])

exit()
###################################################################################################################
## Compare evolution of MLDs with NON-LOCAL MIXING vs LOCAL MIXING using hL thresholds ranging [3000,3500,4000]  ## 
###################################################################################################################

fig, ax = plt.subplots(2,1)

ax[0].plot(df_dats,np.negative(df_mld), color='maroon',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[0].plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[0].plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
ax[0].plot(time2LaRibhLPveP3000,mldKPPLaRibhLPveP3000,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_ylabel('MLD (NON-LOCAL enabled) -m');
ax[0].set_xlim(480,650)
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(df_dats,np.negative(df_mld), color='maroon',linestyle='solid',label = 'OSMOSIS moving average n=5')
ax[1].plot(time2L,mldKPP_Local,color='olivedrab',linestyle='solid',label = 'GOTM KPP')
ax[1].plot(time2LaRibL,mldKPPLaRib_Local,color='steelblue',linestyle='solid',label = 'KPP Lat Rib Cm=Ch=1')
ax[1].plot(time2LaRibhLPveP3000L,mldKPPLaRibhLPveP3000_Local,color='darkorange',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_ylabel('MLD (LOCAL) -m');
ax[1].set_xlim(480,650)
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
########################################################################################
## Evolution of tdotus, u* and heatflux using hL thresholds ranging [3000,3500,4000]  ##
########################################################################################

fig, ax = plt.subplots(3,1)

ax[0].plot(time2LaRibhLPveP4000L,GTM_fric_LaRibhLPveP4000L,color='black',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[0].plot(time2LaRibhLPveP4000L_2030,GTM_fric_LaRibhLPveP4000L_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_ylabel(r'$u_{*}^2$ - $ms^{-1}$');
# ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(time2LaRibhLPveP4000L,GTM_Theatflux_LaRibhLPveP4000L,color='black',linestyle='solid',label = r'KPP original  $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[1].plot(time2LaRibhLPveP4000L_2030,GTM_Theatflux_LaRibhLPveP4000L_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_ylabel(r'$[H + \lambda E + R]$ - $Wm^{-2}$');
# ax[1].set_xlim(490, 635); ax[1].set_ylim(-150,0)

ax[2].plot(time2LaRibhLPveP4000L,GTM_us_LaRibhLPveP4000L,color='black',linestyle='solid',label = r'KPP $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[2].plot(time2LaRibhLPveP4000L_2030,GTM_us_LaRibhLPveP4000L_2030,color='darkorange',linestyle='solid',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_ylabel(r'$u_s$ - $ms^{-1}$');
# ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
# fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
# plt.suptitle('Comparison of MLDs using hL thresholds ranging [3000,3500,4000]')
plt.show()

exit()
##############################################################################
## Comparison of hL time series w & wout SMA ranging [2000,3000,3500,4000]  ##
##############################################################################


fig, ax = plt.subplots(4,1)

ax[0].plot(time2LaRibhLPveP2000,GTM_hLR_LaRibhLPveP2000,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>2000$')
ax[0].plot(time2LaRibhLPveP2000_2030,GTM_hLR_LaRibhLPveP2000_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>2000$')
ax[0].grid(); ax[0].legend(loc='best',fontsize='small')
ax[0].set_xlim(490,635); ax[0].set_ylim(-150,0)

ax[1].plot(time2LaRibhLPveP3000,GTM_hLR_LaRibhLPveP3000,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[1].plot(time2LaRibhLPveP3000_2030,GTM_hLR_LaRibhLPveP3000_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>3000$')
ax[1].grid(); ax[1].legend(loc='best',fontsize='small')
ax[1].set_xlim(490,635); ax[1].set_ylim(-150,0)

ax[2].plot(time2LaRibhLPveP3500,GTM_hLR_LaRibhLPveP3500,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
ax[2].plot(time2LaRibhLPveP3500_2030,GTM_hLR_LaRibhLPveP3500_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>3500$')
ax[2].grid(); ax[2].legend(loc='best',fontsize='small')
ax[2].set_xlim(490,635); ax[2].set_ylim(-150,0)

ax[3].plot(time2LaRibhLPveP4000,GTM_hLR_LaRibhLPveP4000,color='steelblue',linestyle='dashed',label = r'KPP original \ $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[3].plot(time2LaRibhLPveP4000_2030,GTM_hLR_LaRibhLPveP4000_2030,color='darkorange',linestyle='dashed',label = r'KPP SMA=2 wks $\gamma=1 \ \alpha=0 \ for \ hL>4000$')
ax[3].grid(); ax[3].legend(loc='best',fontsize='small')
ax[3].set_xlim(490,635); ax[3].set_ylim(-150,0)

fig.text(0.5, 0.04, 'Time -Julian days', ha='center', va='center')
fig.text(0.08, 0.5, 'Mixed Layer Depth -m', ha='center', va='center', rotation='vertical')
ax[0].text(0.06, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
ax[1].text(0.06, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)
ax[2].text(0.06, 0.08, '(c)', horizontalalignment='center', verticalalignment='center', transform=ax[2].transAxes)
ax[3].text(0.06, 0.08, '(d)', horizontalalignment='center', verticalalignment='center', transform=ax[3].transAxes)

# plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')

# plt.legend(loc='best',fontsize='small')
plt.suptitle('Comparison of hL time series for simualtions with & without SMA ranging [2000,3000,3500,4000]')
plt.show()

exit()
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
plt.show()

exit()
### Plots-compare all MLDs ###

# plt.plot(j,np.negative(i), color='gold',linestyle='solid',label = 'OSMOSIS');
# plt.plot(df_dats,np.negative(df_mld), color='green',linestyle='solid',label = 'OSMOSIS moving average n=3'); 
plt.plot(df_dats5,np.negative(df_mld5), color='gold',linestyle='solid',label = 'OSMOSIS moving average n=5'); 
# plt.plot(df_dats7,np.negative(df_mld7), color='red',linestyle='solid',label = 'OSMOSIS moving average n=7'); 
# plt.plot(new_time,mldOS, color='salmon',linestyle='dotted',label = 'OSMOSIS-calculated'); plt.legend(loc=(0.10,0.15)); 
plt.plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = 'GOTM KPP');
# plt.plot(time2La,mldKPPLa,color='gray',linestyle='dotted',label = 'KPP EV');
plt.plot(time2LaRib,mldKPPLaRib,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1');
# plt.plot(time2LaRibhL,mldKPPLaRibhL,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=1');
# plt.plot(time2LaRibhLcr0,mldKPPLaRibhLcr0,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=0');
# plt.plot(time2LaRibhLcrN100,mldKPPLaRibhLcrN100,color='violet',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=-100');
# plt.plot(time2LaRibhLcr100,mldKPPLaRibhLcr100,color='maroon',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=100');
# plt.plot(time2LaRibhLcr200a0,mldKPPLaRibhLcr200a0,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=200 alpha=0');
# plt.plot(time2LaRibhLcr200,mldKPPLaRibhLcr200,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=200');
# plt.plot(time2LaRibhLPve,mldKPPLaRibhLPve,color='black',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr > 0');
# plt.plot(time2LaRibhLNve,mldKPPLaRibhLNve,color='maroon',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < 0');
# plt.plot(time2LaRibhLPveN100,mldKPPLaRibhLPveN100,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -100');
# plt.plot(time2LaRibhLPveN200,mldKPPLaRibhLPveN200,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -200');
# plt.plot(time2LaRibhLPveN500,mldKPPLaRibhLPveN500,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -500');
# plt.plot(time2LaRibhLa1g2PveN500,mldKPPLaRibhLa1g2PveN500,color='violet',linestyle='dashed',label = 'KPP Lat Rib Cm=2, Ch=1 w hL cr < -500');
plt.plot(time2LaRibhLa2g1PveN500,mldKPPLaRibhLa2g1PveN500,color='black',linestyle='dashed',label = 'KPP Lat Rib Cm=1, Ch=2 w hL cr < -500');
# plt.plot(time2LaRibhLPveN1000,mldKPPLaRibhLPveN1000,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr < -1000');
# plt.plot(time2LaRibhLa1g2PveN1000,mldKPPLaRibhLa1g2PveN1000,color='black',linestyle='dashed',label = 'KPP Lat Rib Cm=2, Ch=1 w hL cr < -1000');
# plt.plot(time2LaRibhLa2g1PveN1000,mldKPPLaRibhLa2g1PveN1000,color='black',linestyle='dashed',label = 'KPP Lat Rib Cm=1, Ch=2 w hL cr < -1000');
# plt.plot(time2LaRibhLN,mldKPPLaRibhLN,color='maroon',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w -hL cr=1');
# plt.plot(time2LaRibhLN05,mldKPPLaRibhLN05,color='limegreen',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w -hL cr=0.5');
# plt.plot(time2LaRibhL05,mldKPPLaRibhL05,color='brown',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w hL cr=0.5');
# plt.plot(time2LaRibhLN2,mldKPPLaRibhLN2,color='violet',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1 w -hL cr=2');
# plt.plot(time2LaRiba2,mldKPPLaRiba2,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=1;Ch=2');
# plt.plot(time2LaRiba025,mldKPPLaRiba025,color='lawngreen',linestyle='dashed',label = 'KPP Lat Rib Cm=1;Ch=0.25');
# plt.plot(time2LaRiba05,mldKPPLaRiba05,color='salmon',linestyle='dashed',label = 'KPP Lat Rib Cm=1;Ch=0.5');
# plt.plot(time2LaRiba075,mldKPPLaRiba075,color='mediumpurple',linestyle='dashed',label = 'KPP Lat Rib Cm=1;Ch=0.75');
# plt.plot(time2LaRiba125,mldKPPLaRiba125,color='lawngreen',linestyle='dashed',label = 'KPP Lat Rib Cm=1;Ch=1.25');
# plt.plot(time2LaRibg2,mldKPPLaRibg2,color='steelblue',linestyle='dashed',label = 'KPP Lat Rib Cm=2;Ch=1');
# plt.plot(time2LaRibg025,mldKPPLaRibg025,color='darkorange',linestyle='dashed',label = 'KPP Lat Rib Cm=0.25;Ch=1');
# plt.plot(time2LaRibg05,mldKPPLaRibg05,color='maroon',linestyle='dashed',label = 'KPP Lat Rib Cm=0.5;Ch=1');
# plt.plot(time2LaRibg075,mldKPPLaRibg075,color='teal',linestyle='dashed',label = 'KPP Lat Rib Cm=0.75;Ch=1');
# plt.plot(time2LaRibg125,mldKPPLaRibg125,color='violet',linestyle='dashed',label = 'KPP Lat Rib Cm=1.25;Ch=1');
# plt.plot(time2_heat,mldKPPheat,color='gray',linestyle='dashed',label = 'KPP HD');
# plt.plot(time2_heat125dlta,mldKPPheat125dlta,color='gray',linestyle=(0, (3,1,1,1)),label = 'KPP HD 8*dlt');
# plt.plot(time_11,mldKPP11,color='salmon',linestyle='solid',label = 'KPP num=1; nuh=1');
# plt.plot(time_12,mldKPP12,color='steelblue',linestyle='solid',label = 'KPP num=1; nuh=2 oriG');
# plt.plot(time_R12,mldKPPR12,color='slateblue',linestyle='dashed',label = 'KPP num=1; nuh=2 mod RiG');
# plt.plot(time_RO12,mldKPPRO12,color='violet',linestyle='dotted',label = 'KPP num=1; nuh=2 mod Ric=0.6');
# plt.plot(time_RR12,mldKPPRR12,color='lawngreen',linestyle='dotted',label = 'KPP num=1; nuh=2 mod Ric=0.15');
# plt.plot(time_21,mldKPP21,color='gray',linestyle='solid',label = 'KPP num=2; nuh=1');
# plt.plot(time_R21,mldKPPR21,color='gray',linestyle='dashed',label = 'KPP num=2; nuh=1 mod RiG');
# plt.plot(time_RO21,mldKPPRO21,color='chocolate',linestyle='dotted',label = 'KPP num=2; nuh=1 mod Ric=0.15');
# plt.plot(time_RR21,mldKPPRR21,color='lawngreen',linestyle='dotted',label = 'KPP num=2; nuh=1 mod Ric=0.6');
# plt.plot(time2,mldKeps,color='goldenrod',linestyle='solid',label = 'GOTM k-eps');
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize='small');
plt.legend(loc='best',fontsize='small')
plt.grid()
plt.xlabel('Time -Julian days'); plt.ylabel('Mixed Layer Depth -m')
# plt.xlim(265,635)
plt.xlim(490,635)
plt.ylim(-150,0)
# plt.ylim(-100,0)
# plt.xlim(535,635)
# plt.xlim(440,540) spring
# plt.xlim(480,540)
# plt.axis([265,635,-1000,0]);
# plt.axis([440,540,400,0]);
# plt.xlim(445.0,537.0)
# plt.axis([445.0,537.0,-150,0]);# plt.title('Compare OSMOSIS & GOTM-KPP,-k-epsilon MLDs')
# plt.tight_layout(rect=[0,0,0.75,1])
# plt.title(' z_ref=-11.0; DeltaT=0.2 ')
plt.savefig('GOTM_OSMOSIS__Mixed_Layer_Depth_ref11.png'); 
plt.show(); #,bbox_inches='tight'
exit()

### --------------------------- quantitative analysis --------------------------- ###
# # remove the nan values 
# GldrtempIdxNaN = np.argwhere(np.isnan(GldrTemp.values));# print('where are the nan values in glider MLD temp', tempIdxNaN)
# GldrtempIdxRowNaN = np.unique(GldrtempIdxNaN[:,0]);# print('row indices of Glider temp NaN values', tempIdxRowNaN)
# new_GldrTemp = np.delete(GldrTemp, GldrtempIdxRowNaN, axis=0);
# Glider_newtime = np.delete(new_time, GldrtempIdxRowNaN, axis=0);

# GldrtempIdxColNaN = np.unique(GldrtempIdxNaN[:,:]);# print('row indices of Glider temp NaN values', tempIdxRowNaN)
# newC_GldrTemp = np.delete(GldrTemp, GldrtempIdxColNaN, axis=1);
# GliderC_depth = np.delete(GldrDepth, GldrtempIdxColNaN, axis=0);
# interpolate the matching times for OSMOSIS and GOTM temperature
# GliderTempInterp = interpolate.interp2d(GldrDepth,Glider_newtime,new_GldrTemp)

print('OSMOSIS MLD SMA n-5 :',np.negative(df_mld5).shape, df_dats5.shape)
print('GOTM KPP MLD :',mldKPP.shape, time2.shape)
print('GOTM KPP MLD Lat Rib Cm=Ch=1 :',mldKPPLaRib.shape, time2LaRib.shape)
# print('deleted nan in Temp col',newC_GldrTemp.values,newC_GldrTemp.shape, GliderC_depth.shape)
# print('deleted nan in Temp row',new_GldrTemp.shape, Glider_newtime.shape)
# print(new_time[309:381].shape, GldrDepth.shape, GldrTemp[309:381,:].shape )
# print(time2.shape, GTMdepth.shape, GTMtemp.shape )

# interpolate GOTM time and 
GTMmldInterp = np.interp(df_dats5,time2,mldKPP)
# print('interpolated GOTM MLD', GTMmldInterp.shape, GTMmldInterp)
# print('MLD time : ', df_dats5)
OStimeSpr = np.where((df_dats5 <= 538.0) & (df_dats5 >= 445.0))
print('OSMOSIS spring time elements :', OStimeSpr)


exit()
GTMmldLaRigInterp = np.interp(df_dats5,time2LaRib,mldKPPLaRib)
# print('interpolated GOTM MLD', GTMmldLaRigInterp.shape, GTMmldLaRigInterp)

correlation_OS_GTM_MLD = sum(np.multiply(np.negative(df_mld5), GTMmldInterp))/np.sqrt(np.sum(np.square(np.negative(df_mld5)))*np.sum(np.square(GTMmldInterp)))
print('correlation coefficient between OSMOSIS and GOTM KPP mixed layer depth :', correlation_OS_GTM_MLD)

correlation_OS_GTMLaRig_MLD = sum(np.multiply(np.negative(df_mld5), GTMmldLaRigInterp))/np.sqrt(np.sum(np.square(np.negative(df_mld5)))*np.sum(np.square(GTMmldLaRigInterp)))
print('correlation coefficient between OSMOSIS and GOTM KPP La Rig mixed layer depth :', correlation_OS_GTMLaRig_MLD)
exit()

### ----------------- PLOTS -compare GOTM KPP and k-eps velocities ---------------------- ###

# a = [499]

# for i in a:
#   kppu1d = GTMu.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
#   kpptime = time2[:] 
#   plt.plot(kpptime,kppu1d,label='GTM kpp depth= %s'%kppdepth.values);# plt.legend(loc=(0.66,0.02),fontsize ='xsmall')
#   plt.legend(loc='best',fontsize ='small')
#   plt.xlabel('Time -Julian days');plt.ylabel('Streamwise velocity -m/s'); 
#   plt.title('Annual streamwise velocity time series')
#   # plt.axis([11.0,15.5,-150,0])
# plt.savefig('osmosis_annual_variablePlots/GOTM_streamwise_velocity_timeseries.png')
# # plt.show()

# for i in a:
#   kppv1d = GTMv.isel(lon=0,lat=0)[:,i]; kppdepth = GTMdepth[i]
#   kpptime = time2[:] 
#   plt.plot(kpptime,kppv1d,label='GTM kpp depth= %s'%kppdepth.values);# plt.legend(loc=(0.66,0.02),fontsize ='xsmall')
#   plt.legend(loc='best',fontsize ='small')
#   plt.xlabel('Time -Julian days');plt.ylabel('Spanwise velocity -m/s'); 
#   plt.title('Annual spanwise velocity time series')
#   # plt.axis([11.0,15.5,-150,0])
# plt.savefig('osmosis_annual_variablePlots/GOTM_spanwise_velocity_timeseries.png')
# plt.show()


### ----------------- PLOTS - GOTM KPP and k-eps other variables ---------------------- ###
print(time2[0:500].values, time2.shape)
GTM_sel = [0];# [0,13515]
Gldr_sel = [2127] #[2127, ,3255]
# SALINITY
for j in Gldr_sel:
	GldrSal = GldrSalt.isel()[j,:]
	Gldrdepth = GldrDepth[:]; Gldrtime = new_time[j]
	for i in GTM_sel:
		kppsal1d = GTMsal.isel(lon=0,lat=0)[i,:];
		kppdepth = GTMdepth[:]
		kpptime = time2[i]
plt.plot(kppsal1d,kppdepth,'seagreen',label = 'GOTM salinity %0.2f'%kpptime.values);
plt.plot(GldrSal,Gldrdepth,'gold',label = 'OSMOSIS salinity %0.2f'%Gldrtime.values);
plt.legend(loc='best', fontsize='small')
plt.ylabel('Depth -m');plt.xlabel('Salinity -psu'); 
plt.ylim(-200,0)
plt.grid(); plt.show()

# OXYGEN
for j in Gldr_sel:
	GldrO2 = Gldro2.isel()[j,:]
	Gldrdepth = GldrDepth[:]; Gldrtime = new_time[j]
	for i in GTM_sel:
		kppyo21d = GTMo2.isel(lon=0,lat=0)[i,:]; 
		kppdepth = GTMdepth[:]
		kpptime = time2[i]
plt.plot(kppyo21d,kppdepth,'seagreen',label = 'GOTM oxygen %0.2f'%kpptime.values);
plt.plot(GldrO2,Gldrdepth,'gold',label = 'OSMOSIS oxygen %0.2f'%Gldrtime.values);
plt.legend(loc='best', fontsize='small')
plt.ylabel('Depth -m');plt.xlabel('Oxygen -mmol/m3'); 
plt.grid()
# plt.title('KPP tracer variable profiles for %s'%kpptime.values)
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-200,0)
# plt.savefig('GOTM_kpp_keps_eddy_viscosity_num_025_profiles.png')
plt.show()
exit()

### ----------------- PLOTS -compare GOTM KPP and k-eps viscosity ---------------------- ###
print(time2[0:500].values, time2.shape)
GTM_sel = [0];# [0,52500]

# EDDY VISCOSITY
for i in GTM_sel:
    kppnum1d = GTMnum.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
    kppnumLa = GTMnumLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
    kppnum_heat = GTMnum_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
    # kppnum_mh = GTMnum_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
    # kppnum_heatg5 = GTMnum_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
    # kppnum_heatg15 = GTMnum_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
    kepsnum1d = GTMkenum.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
    kppnum025025 = GTMnum025025.isel(lon=0,lat=0)[i,:]; kppdepth025025 = GTMdepth025025[:]
    kppnum02505 = GTMnum02505.isel(lon=0,lat=0)[i,:]; kppdepth02505 = GTMdepth02505[:]
    kppnum0251 = GTMnum0251.isel(lon=0,lat=0)[i,:]; kppdepth0251 = GTMdepth0251[:]
    kppnum0252 = GTMnum0252.isel(lon=0,lat=0)[i,:]; kppdepth0252 = GTMdepth0252[:]
    kppnum0254 = GTMnum0254.isel(lon=0,lat=0)[i,:]; kppdepth0254 = GTMdepth0254[:]
    kppnum05025 = GTMnum05025.isel(lon=0,lat=0)[i,:]; kppdepth05025 = GTMdepth05025[:]
    kppnum0505 = GTMnum0505.isel(lon=0,lat=0)[i,:]; kppdepth0505 = GTMdepth0505[:]
    kppnum051 = GTMnum051.isel(lon=0,lat=0)[i,:]; kppdepth051 = GTMdepth051[:]
    kppnum052 = GTMnum052.isel(lon=0,lat=0)[i,:]; kppdepth052 = GTMdepth052[:]
    kppnum054 = GTMnum054.isel(lon=0,lat=0)[i,:]; kppdepth054 = GTMdepth054[:]
    kppnum1025 = GTMnum1025.isel(lon=0,lat=0)[i,:]; kppdepth1025 = GTMdepth1025[:]
    kppnum105 = GTMnum105.isel(lon=0,lat=0)[i,:]; kppdepth105 = GTMdepth105[:]
    kppnum11 = GTMnum11.isel(lon=0,lat=0)[i,:]; kppdepth11 = GTMdepth11[:]
    kppnum12 = GTMnum12.isel(lon=0,lat=0)[i,:]; kppdepth12 = GTMdepth12[:]
    kppnum14 = GTMnum14.isel(lon=0,lat=0)[i,:]; kppdepth14 = GTMdepth14[:]
    kppnum2025 = GTMnum2025.isel(lon=0,lat=0)[i,:]; kppdepth2025 = GTMdepth2025[:]
    kppnum205 = GTMnum205.isel(lon=0,lat=0)[i,:]; kppdepth205 = GTMdepth205[:]
    kppnum21 = GTMnum21.isel(lon=0,lat=0)[i,:]; kppdepth21 = GTMdepth21[:]
    kppnum22 = GTMnum22.isel(lon=0,lat=0)[i,:]; kppdepth22 = GTMdepth22[:]
    kppnum24 = GTMnum24.isel(lon=0,lat=0)[i,:]; kppdepth24 = GTMdepth24[:]
    kppnum4025 = GTMnum4025.isel(lon=0,lat=0)[i,:]; kppdepth4025 = GTMdepth4025[:]
    kppnum405 = GTMnum405.isel(lon=0,lat=0)[i,:]; kppdepth405 = GTMdepth405[:]
    kppnum41 = GTMnum41.isel(lon=0,lat=0)[i,:]; kppdepth41 = GTMdepth41[:]
    kppnum42 = GTMnum42.isel(lon=0,lat=0)[i,:]; kppdepth42 = GTMdepth42[:]
    kppnum44 = GTMnum44.isel(lon=0,lat=0)[i,:]; kppdepth44 = GTMdepth44[:]
    kppnum58025 = GTMnum58025.isel(lon=0,lat=0)[i,:]; kppdepth58025 = GTMdepth58025[:]
    kppnum5805 = GTMnum5805.isel(lon=0,lat=0)[i,:]; kppdepth5805 = GTMdepth5805[:]
    kppnum581 = GTMnum581.isel(lon=0,lat=0)[i,:]; kppdepth581 = GTMdepth581[:]
    # kppnum582 = GTMnum582.isel(lon=0,lat=0)[i,:]; kppdepth582 = GTMdepth582[:]
    kppnum584 = GTMnum584.isel(lon=0,lat=0)[i,:]; kppdepth584 = GTMdepth584[:]
    kpptime = time2[i]
    kpptimeLa = time2La[i]
    kpptime_heat = time2_heat[i]
    # kpptime_mh = time2_mh[i]
    # kpptime_heatg5 = time2_heatg5[i]
    # kpptime_heatg15 = time2_heatg15[i]
    kepstime = ketime2[i]
    kpptime025025 = time_025025[i]
    kpptime02505 = time_02505[i]
    kpptime0251 = time_0251[i]
    kpptime0252 = time_0252[i]
    kpptime0254 = time_0254[i]
    kpptime05025 = time_05025[i]
    kpptime0505 = time_0505[i]
    kpptime051 = time_051[i]
    kpptime052 = time_052[i]
    kpptime054 = time_054[i]
    kpptime1025 = time_1025[i]
    kpptime105 = time_105[i]
    kpptime11 = time_11[i]
    kpptime12 = time_12[i]
    kpptime14 = time_14[i]
    kpptime2025 = time_2025[i]
    kpptime205 = time_205[i]
    kpptime21 = time_21[i]
    kpptime22 = time_22[i]
    kpptime24 = time_24[i]
    kpptime4025 = time_4025[i]
    kpptime405 = time_405[i]
    kpptime41 = time_41[i]
    kpptime42 = time_42[i]
    kpptime44 = time_44[i]
    kpptime58025 = time_58025[i]
    kpptime5805 = time_5805[i]
    kpptime581 = time_581[i]
    kpptime582 = time_582[i]
    kpptime584 = time_584[i]
plt.plot(kppnum1d,kppdepth,'k-',label = 'KPP'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnumLa,kppdepthLa,'b-' ,label = 'KPP EV'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum_heat,kppdepth_heat,'g--',label = 'KPP HD'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_mh,kppdepth_mh,'r--',label = 'KPP EV-HD'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heatg5,kppdepth_heatg5,'r-',label = 'GOTM_KPP diff gam=0.5 day %s'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnum_heatg15,kppdepth_heatg15,'r--',label = 'GOTM_KPP diff gam=1.5 day %s'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnum1d,kepsdepth,'k:' ,label = 'k-eps'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum025025,kppdepth025025,color='steelblue',linestyle=(0,(1,1)),label = 'KPP num=0.25, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum02505,kppdepth02505,color='steelblue',linestyle=(0,(5,1)),label = 'KPP num=0.25, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum0251,kppdepth0251,color='steelblue', linestyle='solid',label = 'KPP num=0.25, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum0252,kppdepth0252,color='steelblue',linestyle=(0,(3,1,1,1)),label = 'KPP num=0.25, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum0254,kppdepth0254,color='steelblue',linestyle=(0,(3,1,1,1,1,1)),label = 'KPP num=0.25, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum05025,kppdepth05025,color='olivedrab',linestyle=(0,(1,1)),label = 'KPP num=0.5, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum0505,kppdepth0505,color='olivedrab',linestyle=(0,(5,1)),label = 'KPP num=0.5, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum051,kppdepth051,color='olivedrab',linestyle='solid',label = 'KPP num=0.5, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum052,kppdepth052,color='olivedrab',linestyle=(0,(3,1,1,1)),label = 'KPP num=0.5, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum054,kppdepth054,color='olivedrab',linestyle=(0,(3,1,1,1,1,1)),label = 'KPP num=0.5, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum1025,kppdepth1025,color='salmon',linestyle=(0,(1,1)),label = 'KPP num=1, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum105,kppdepth105,color='salmon',linestyle=(0,(5,1)),label = 'KPP num=1, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum11,kppdepth11,color='salmon',linestyle='solid',label = 'KPP num=1, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum12,kppdepth12,color='salmon',linestyle=(0,(3,1,1,1)),label = 'KPP num=1, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum14,kppdepth14,color='salmon',linestyle=(0,(3,1,1,1,1,1)),label = 'KPP num=1, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum2025,kppdepth2025,color='darkorange',linestyle=(0,(1,1)),label = 'KPP num=2, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum205,kppdepth205,color='darkorange',linestyle=(0,(5,1)),label = 'KPP num=2, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum21,kppdepth21,color='darkorange',linestyle='solid',label = 'KPP num=2, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum22,kppdepth22,color='darkorange',linestyle=(0,(3,1,1,1)),label = 'KPP num=2, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum24,kppdepth24,color='darkorange',linestyle=(0,(3,1,1,1,1,1)),label = 'KPP num=2, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum4025,kppdepth4025,color='rosybrown',linestyle=(0,(1,1)),label = 'KPP num=4, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum405,kppdepth405,color='rosybrown',linestyle=(0,(5,1)),label = 'KPP num=4, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum41,kppdepth41,color='rosybrown',linestyle='solid',label = 'KPP num=4, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum42,kppdepth42,color='rosybrown',linestyle=(0,(3,1,1,1)),label = 'KPP num=4, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum44,kppdepth44,color='rosybrown',linestyle=(0,(3,1,1,1,1,1)),label = 'KPP num=4, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum58025,kppdepth58025,color='seagreen',linestyle=(0,(1,1)),label = 'KPP num=5.8, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum5805,kppdepth5805,color='seagreen',linestyle=(0,(5,1)),label = 'KPP num=5.8, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum581,kppdepth581,color='seagreen',linestyle='solid',label = 'KPP num=5.8, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum582,kppdepth582,color='seagreen',linestyle=(0,(3,1,1,1)),label = 'KPP num=5.8, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnum584,kppdepth584,color='seagreen',linestyle=(0,(3,1,1,1,1,1)),label = 'KPP num=5.8, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.xlabel('Eddy Viscosity -m2/s');plt.ylabel('Depth -m'); 
plt.grid()
plt.title('KPP,-k-eps eddy viscosity profiles for %s'%kpptime.values)
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-100,0)
# plt.savefig('GOTM_kpp_keps_eddy_viscosity_num_025_profiles.png')
plt.show()
exit()

# HEAT DIFFUSIVITY
# for i in GTM_sel:
#     kppnuh1d = GTMnuh.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
#     kppnuhLa = GTMnuhLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
#     kppnuh_heat = GTMnuh_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
#     # kppnuh_mh = GTMnuh_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
#     # kppnuh_heatg5 = GTMnuh_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
#     # kppnuh_heatg15 = GTMnuh_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
#     kepsnuh1d = GTMkenuh.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
#     kppnuh025025 = GTMnuh025025.isel(lon=0,lat=0)[i,:]; kppdepth025025 = GTMdepth025025[:]
#     kppnuh02505 = GTMnuh02505.isel(lon=0,lat=0)[i,:]; kppdepth02505 = GTMdepth02505[:]
#     kppnuh0251 = GTMnuh0251.isel(lon=0,lat=0)[i,:]; kppdepth0251 = GTMdepth0251[:]
#     kppnuh0252 = GTMnuh0252.isel(lon=0,lat=0)[i,:]; kppdepth0252 = GTMdepth0252[:]
#     kppnuh0254 = GTMnuh0254.isel(lon=0,lat=0)[i,:]; kppdepth0254 = GTMdepth0254[:]
#     kpptime = time2[i]
#     kpptimeLa = time2La[i]
#     kpptime_heat = time2_heat[i]
#     # kpptime_mh = time2_mh[i]
#     # kpptime_heatg5 = time2_heatg5[i]
#     # kpptime_heatg15 = time2_heatg15[i]
#     kepstime = ketime2[i]
#     kpptime025025 = time_025025[i]
#     kpptime02505 = time_02505[i]
#     kpptime0251 = time_0251[i]
#     kpptime0252 = time_0252[i]
#     kpptime0254 = time_0254[i]
# plt.plot(kppnuh1d,kppdepth,'k-',label = 'KPP'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuhLa,kppdepthLa,'b-' ,label = 'KPP EV'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heat,kppdepth_heat,'g--',label = 'KPP HD'); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_mh,kppdepth_mh,'r--',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_heatg5,kppdepth_heatg5,'r-',label = 'GOTM_KPP diff gam=0.5 day %s'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# # plt.plot(kppnuh_heatg15,kppdepth_heatg15,'r--',label = 'GOTM_KPP diff gam=1.5 day %s'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kepsnuh1d,kepsdepth,'k:' ,label = 'k-eps'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh025025,kppdepth025025,color='steelblue',label = 'KPP num=0.25, nuh=0.25'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh02505,kppdepth02505,color='olivedrab',label = 'KPP num=0.25, nuh=0.5'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh0251,kppdepth0251,color='salmon',label = 'KPP num=0.25, nuh=1'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh0252,kppdepth0252,color='darkorange',label = 'KPP num=0.25, nuh=2'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh0254,kppdepth0254,color='rosybrown',label = 'KPP num=0.25, nuh=4'); plt.legend(loc='best', fontsize='small')
# plt.xlabel('Heat Diffusivity -m2/s');plt.ylabel('Depth -m'); 
# plt.grid()
# plt.title('KPP,-k-eps heat diffusivity for %s'%kpptime.values)
# # plt.axis([-0.01,0.15,-100,0])
# plt.ylim(-100,0)
# plt.savefig('GOTM_kpp_keps_heat_diffusivity_nuh_025_profiles.png')
# plt.show()

for i in GTM_sel:
    kppnuh1d = GTMnuh.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
    kppnuhLa = GTMnuhLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
    kppnuh_heat = GTMnuh_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
    # kppnuh_mh = GTMnuh_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
    # kppnuh_heatg5 = GTMnuh_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
    # kppnuh_heatg15 = GTMnuh_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
    kepsnuh1d = GTMkenuh.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
    kppnuh05025 = GTMnuh05025.isel(lon=0,lat=0)[i,:]; kppdepth05025 = GTMdepth05025[:]
    kppnuh0505 = GTMnuh0505.isel(lon=0,lat=0)[i,:]; kppdepth0505 = GTMdepth0505[:]
    kppnuh051 = GTMnuh051.isel(lon=0,lat=0)[i,:]; kppdepth051 = GTMdepth051[:]
    kppnuh052 = GTMnuh052.isel(lon=0,lat=0)[i,:]; kppdepth052 = GTMdepth052[:]
    kppnuh054 = GTMnuh054.isel(lon=0,lat=0)[i,:]; kppdepth054 = GTMdepth054[:]
    kpptime = time2[i]
    kpptimeLa = time2La[i]
    kpptime_heat = time2_heat[i]
    # kpptime_mh = time2_mh[i]
    # kpptime_heatg5 = time2_heatg5[i]
    # kpptime_heatg15 = time2_heatg15[i]
    kepstime = ketime2[i]
    kpptime05025 = time_05025[i]
    kpptime0505 = time_0505[i]
    kpptime051 = time_051[i]
    kpptime052 = time_052[i]
    kpptime054 = time_054[i]
plt.plot(kppnuh1d,kppdepth,'k-',label = 'KPP'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa,kppdepthLa,'b-' ,label = 'KPP EV'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh_heat,kppdepth_heat,'g--',label = 'KPP HD'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh,kppdepth_mh,'r--',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg5,kppdepth_heatg5,'r-',label = 'GOTM_KPP diff gam=0.5 day %s'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg15,kppdepth_heatg15,'r--',label = 'GOTM_KPP diff gam=1.5 day %s'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnuh1d,kepsdepth,'k:' ,label = 'k-eps'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh05025,kppdepth05025,color='steelblue',label = 'KPP num=0.5, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh0505,kppdepth0505,color='olivedrab',label = 'KPP num=0.5, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh051,kppdepth051,color='salmon',label = 'KPP num=0.5, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh052,kppdepth052,color='darkorange',label = 'KPP num=0.5, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh054,kppdepth054,color='rosybrown',label = 'KPP num=0.5, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.xlabel('Heat Diffusivity -m2/s');plt.ylabel('Depth -m'); 
plt.grid()
plt.title('KPP,-k-eps heat diffusivity for %s'%kpptime.values)
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-100,0)
plt.savefig('GOTM_kpp_keps_heat_diffusivity_nuh_05_profiles.png')
plt.show()

for i in GTM_sel:
    kppnuh1d = GTMnuh.isel(lon=0,lat=0)[i,:]; kppdepth = GTMdepth[:]
    kppnuhLa = GTMnuhLa.isel(lon=0,lat=0)[i,:]; kppdepthLa = GTMdepthLa[:]
    kppnuh_heat = GTMnuh_heat.isel(lon=0,lat=0)[i,:]; kppdepth_heat = GTMdepth_heat[:]
    # kppnuh_mh = GTMnuh_mh.isel(lon=0,lat=0)[i,:]; kppdepth_mh = GTMdepth_mh[:]
    # kppnuh_heatg5 = GTMnuh_heatg5.isel(lon=0,lat=0)[i,:]; kppdepth_heatg5 = GTMdepth_heatg5[:]
    # kppnuh_heatg15 = GTMnuh_heatg15.isel(lon=0,lat=0)[i,:]; kppdepth_heatg15 = GTMdepth_heatg15[:]
    kepsnuh1d = GTMkenuh.isel(lon=0,lat=0)[i,:]; kepsdepth = GTMkedepth[:]
    kppnuh1025 = GTMnuh1025.isel(lon=0,lat=0)[i,:]; kppdepth1025 = GTMdepth1025[:]
    kppnuh105 = GTMnuh105.isel(lon=0,lat=0)[i,:]; kppdepth105 = GTMdepth105[:]
    kppnuh11 = GTMnuh11.isel(lon=0,lat=0)[i,:]; kppdepth11 = GTMdepth11[:]
    kppnuh12 = GTMnuh12.isel(lon=0,lat=0)[i,:]; kppdepth12 = GTMdepth12[:]
    kppnuh14 = GTMnuh14.isel(lon=0,lat=0)[i,:]; kppdepth14 = GTMdepth14[:]
    kpptime = time2[i]
    kpptimeLa = time2La[i]
    kpptime_heat = time2_heat[i]
    # kpptime_mh = time2_mh[i]
    # kpptime_heatg5 = time2_heatg5[i]
    # kpptime_heatg15 = time2_heatg15[i]
    kepstime = ketime2[i]
    kpptime1025 = time_1025[i]
    kpptime105 = time_105[i]
    kpptime11 = time_11[i]
    kpptime12 = time_12[i]
    kpptime14 = time_14[i]
plt.plot(kppnuh1d,kppdepth,'k-',label = 'KPP'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuhLa,kppdepthLa,'b-' ,label = 'KPP EV'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh_heat,kppdepth_heat,'g--',label = 'KPP HD'); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_mh,kppdepth_mh,'r--',label = 'GOTM_KPP visc-diff day %s'%kpptime_mh.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg5,kppdepth_heatg5,'r-',label = 'GOTM_KPP diff gam=0.5 day %s'%kpptime_heatg5.values); plt.legend(loc='best', fontsize='small')
# plt.plot(kppnuh_heatg15,kppdepth_heatg15,'r--',label = 'GOTM_KPP diff gam=1.5 day %s'%kpptime_heatg15.values); plt.legend(loc='best', fontsize='small')
plt.plot(kepsnuh1d,kepsdepth,'k:' ,label = 'k-eps'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh1025,kppdepth1025,color='steelblue',label = 'KPP num=1, nuh=0.25'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh105,kppdepth105,color='olivedrab',label = 'KPP num=1, nuh=0.5'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh11,kppdepth11,color='salmon',label = 'KPP num=1, nuh=1'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh12,kppdepth12,color='darkorange',label = 'KPP num=1, nuh=2'); plt.legend(loc='best', fontsize='small')
plt.plot(kppnuh14,kppdepth14,color='rosybrown',label = 'KPP num=1, nuh=4'); plt.legend(loc='best', fontsize='small')
plt.xlabel('Heat Diffusivity -m2/s');plt.ylabel('Depth -m'); 
plt.grid()
plt.title('KPP,-k-eps heat diffusivity for %s'%kpptime.values)
# plt.axis([-0.01,0.15,-100,0])
plt.ylim(-100,0)
plt.savefig('GOTM_kpp_keps_heat_diffusivity_nuh_1_profiles.png')
plt.show()

exit()
### ----------------------------- flux .dat files ----------------------------- ###

hdate = []; htime = []; solarRad = []; surfaceHflux = []

with open('osmosis_annual_surface_forcingKPP/heat_flux.dat', 'r') as f:
    for row in f:
        a, b, c, d = row.split()
        hdate.append(a)
        htime.append(b)
        solarRad.append(c)
        surfaceHflux.append(d)
solarRad = np.asarray(solarRad); solarRad = solarRad.astype('float64')
surfaceHflux = np.asarray(surfaceHflux); surfaceHflux = surfaceHflux.astype('float64')

## ------------------- heat flux times --------------------- ##
concat_func = lambda x,y: str(x) + " " + str(y)
DT = list(map(concat_func,hdate,htime)) # list the map function
DT = np.asarray(DT)

hdate = np.asarray(hdate); htime = np.asarray(htime)
# print(DT.shape)
yearDaysH = np.empty((2929,0), int)
for i in range(len(DT)):
    yearDaysH= np.append(yearDaysH, datetime.strptime(DT[i],'%Y:%m:%d %H:%M:%S').timetuple().tm_yday)
# print(yearDaysH.shape)
hoursH = np.empty((2929,0), int)
for tim in range(len(htime)):
    timeDays = datetime.strptime(htime[tim], '%H:%M:%S')#.split(':')
    hoursH = np.append(hoursH, timeDays.hour)

hourstodaysH = hoursH/24
# print(hourstodaysH)
convertedDaysH = yearDaysH + hourstodaysH

for data in range(len(convertedDaysH)):
    if convertedDaysH[data] < 267:
        convertedDaysH[data] = convertedDaysH[data] + 366
    else: 
        convertedDaysH[data] == convertedDaysH[data]
convertedDaysH[2928]=633.0

convertedDaysH = np.asarray(convertedDaysH);

############################################################################################
### ------------------- sum of solar radiation and surface heat flux ------------------- ###
############################################################################################

sum_radHeatflux = np.empty((2929,0),float)

for i in range(len(solarRad)):
    sum_radHeatflux = np.append(sum_radHeatflux,solarRad[i]+surfaceHflux[i])
# print('sum of radiation and surface heat flux',sum_radHeatflux)

# ------------ constants ------------ #
rho = np.full(shape=2929, fill_value=1025.0, dtype=np.float) # density of sea water kg/m^3
c_pHeat = np.full(shape=2929, fill_value=3850.0, dtype=np.float) # specific heat capacity of sea water J/(kg C)

## ------------------- full sum of solar radiation and surface heat flux --------------------- ##
fullsum_radHeatflux = np.empty((2929,0),float)
for i in range(len(solarRad)):
    fullsum_radHeatflux = np.append(fullsum_radHeatflux,sum_radHeatflux[i]/(rho[i]*c_pHeat[i]))
# print('full sum of solar radiation and heat flux',fullsum_radHeatflux)

###################################################################
### ------------------- surface heat budget ------------------- ###
###################################################################

# print('OSMOSIS temperature times:', new_time[204:522].values)
# print('OSMOSIS temperature depth:', GldrDepth[0:200].values)
# print('OSMOSIS heat_flux times:', convertedDaysH[2900:])
# print('GOTM temp depth:', GTMdepth[300:499])

## selecting the GOTM temperature for the mixed layer depth
dz = np.diff(GTMdepth)[0];# print(dz)

# kpptemp2d = GTMtemp.isel(lat=0,lon=0); print('shape of GOTM temp', kpptemp2d[0:10,399:499].shape)
kpptemp2dMLD = GTMtemp.isel(lat=0,lon=0,z=range(299,500));# print('GOTM temp MLD', kpptemp2dMLD.values)
print('GOTM depth :', GTMdepth[299:500].values, GTMdepth[299:500].shape, GTMdepth.shape)

GTM_h = (GTMdepth[499]-GTMdepth[299]).values; print('GTM mixed layer depth : ',GTM_h)
sumkppTempMLD = (kpptemp2dMLD.sum(axis=1)*dz); print('sum of temp under GOTM MLD', sumkppTempMLD.values, sumkppTempMLD.shape)
# sumkppTempMLD = (kpptemp2dMLD.sum(axis=1)*dz)/GTM_h; print('sum of temp under GOTM MLD', sumkppTempMLD.values, sumkppTempMLD.shape)

heatBudget = []
midtime = []
    # for i in range(1,Nz):
    #   sumTempMLD[i] = np.sum(GTMtempA.isel(lon=0,lat=0)[i,400:500]) 
dt = np.diff(time2)[0]*86400;# print('dt : ', dt)

for n in range(0, 52595):
    CT = (sumkppTempMLD[n+1].values-sumkppTempMLD[n].values)/dt
    # CT = (GTMtempA.isel(lon=0,lat=0)[n+1,:]-GTMtempA.isel(lon=0,lat=0)[n,:])/dt
    diffT = (time2[n+1].values+time2[n].values)/2
    heatBudget.append(CT)
    midtime.append(diffT)
heatBudget = np.asarray(heatBudget); midtime = np.asarray(midtime)
# print('GOTM heatbudget array : ', heatBudget); print('GOTM HB time array : ', midtime)


## selecting the Glider temperature for mixed layer depth
dz_OS = -np.diff(GldrDepth)[0];# print('OSMOSIS glider dz',dz_OS)
dt_OS = np.diff(new_time)*86400;# print('OSMOSIS glider dt', dt_OS[0:10], len(dt_OS))

## remove nan elements in the Glider OSMOSIS data
# gliderTempMLD = GldrTemp.iloc[:,0:200]; print('list of glider temp  mld', gliderTempMLD.values, gliderTempMLD.shape)

gliderTempMLD = GldrTemp.isel(pressure=range(0,201));# print('list of glider temp mld', gliderTempMLD)
np.savetxt('heatFlux_osmosisAnnual/gliderTempMLDwNaN.txt', gliderTempMLD, fmt = '%1.10f');
# print('Glider depth :', GldrDepth[0:201].values, GldrDepth[0:201].shape)

# gliderTempMLD[np.isnan(gliderTempMLD)] = -9999; gliderTempMLD[gliderTempMLD==np.nan] = -9999; gliderTempMLD = gliderTempMLD.fillna(-9999)
# print('list of glider temp mld', gliderTempMLD[0:10,0:200])
tempIdxNaN = np.argwhere(np.isnan(gliderTempMLD.values));# print('where are the nan values in glider MLD temp', tempIdxNaN)
tempIdxRowNaN = np.unique(tempIdxNaN[:,0]);# print('row indices of Glider temp NaN values', tempIdxRowNaN)
new_gliderTempMLD = np.delete(gliderTempMLD, tempIdxRowNaN, axis=0);
new_glidertime = np.delete(new_time, tempIdxRowNaN, axis=0);
# print('compare shapes of original and new glider temp array', gliderTempMLD.shape, new_gliderTempMLD.shape)
# print('compare shapes of original and new glider time array', new_time.shape, new_glidertime.shape)

new_dt_OS = np.diff(new_glidertime)*86400;# print('OSMOSIS glider dt', dt_OS[0:10], len(dt_OS))
Gldr_h = -(GldrDepth[201]-GldrDepth[0]).values; print('Gldr mixed layer depth : ',Gldr_h)

## Vertical integral of temperature ##
# sumGliderTempMLD = (gliderTempMLD.sum(axis=1)*dz_OS);# print('sum of temp under Glider MLD', sumGliderTempMLD.values, sumGliderTempMLD.shape)
sumGliderTempMLD = (gliderTempMLD.sum(axis=1)*dz_OS)/Gldr_h;# print('sum of temp under Glider MLD', sumGliderTempMLD.values, sumGliderTempMLD.shape)
# new_sumGliderTempMLD = (new_gliderTempMLD.sum(axis=1)*dz_OS); print('sum of NEW temp under Glider MLD', new_sumGliderTempMLD.values, new_sumGliderTempMLD.shape)
new_sumGliderTempMLD = (new_gliderTempMLD.sum(axis=1)*dz_OS)/Gldr_h; print('sum of NEW temp under Glider MLD', new_sumGliderTempMLD.values, new_sumGliderTempMLD.shape)

# print('compare glider temp array at depth =-9 m', gliderTempMLD[4,:], new_gliderTempMLD[0,:])
# print('compare glider time array at depth =-9 m', new_time[0:10], new_glidertime[0:10])

# ## ---------------------- quantitative analysis (temperature) ------------------------- ##

## ------------ correlation ------------ ##



# interpolate

kpptemplist = []
GldrTempIntlist = []
correlation_KPP_Glider = []
mean_abs_error = []
root_mean_square_error = []


for depth in range(len(GldrDepth)):
    GldrTemp = new_gliderTempMLD.isel()[:,depth]
    Gldrdepth = GldrDepth[depth]; Gldrtime = new_glidertime[:]
    GldrTempIntlist = np.append(GldrTempIntlist, GldrTemp.values)
    for depthG in range(len(GTMdepth)):
        kpptempinterp = np.interp(new_glidertime,time2,GTMtemp.isel(lon=0,lat=0)[:,depthG]); 
        kppdepth = GTMdepth[depthG]; kpptime = time2[:]
        kpptemplist = np.append(kpptemplist, kpptempinterp)
    sumkppGldr = sum(kpptemplist*GldrTempIntlist)
    sumkpp = sum((kpptemplist)**(2))
    sumGldr = sum((GldrTempIntlist)**(2))
    correlation_KPP_Glider = np.append(correlation_KPP_Glider, sumkppGldr/(sumkpp*sumGldr)**(0.5) )
    mean_abs_error = np.append(mean_abs_error, 1/np.size(kpptemplist)*sum(np.absolute(kpptemplist-GldrTempIntlist)))
    root_mean_square_error = np.append(root_mean_square_error, (1/np.size(kpptemplist)*sum((kpptemplist-GldrTempIntlist)**(2)))**(0.5))
    print('KPP temp interpolated : ',kpptemplist)
    print('OSMOSIS temp',GldrTempIntlist)
    print('correlation coefficient between KPP and OSMOSIS Temperature :', correlation_KPP_Glider)
    print('mean absolute err between KPP and OSMOSIS Temperature :', mean_abs_error)
    print('root mean sqare err between KPP and OSMOSIS Temperature :', root_mean_square_error)
plt.plot(Gldrtime,GldrTempIntlist, '-k', label='OSMOSIS Temperature at z=%0.2f'%Gldrdepth.values)
plt.plot(Gldrtime,kpptemplist, '-b', label=' interpolated GOTM KPP Temperature at z=%0.2f'%kppdepth.values)
plt.legend(loc='best',fontsize='small')
plt.xlabel('Time -Julian days'); plt.ylabel('Temperature -C')
plt.savefig('osmosis_autumn_variablePlots/Gldr_GOTM_Temperature_interpolated_quantitative_timeseries.png');plt.show()

exit()
## ---------------------- comparison plots ------------------------- ##
# # TEMPERATURE-timeseries
# plt.plot(new_time, gliderTempMLD[:,0], 'k-',label='Original OSMOSIS Temperature at depth=%0.2f m'%GldrDepth[0].values); plt.legend(loc='best', fontsize='small')
# plt.plot(time2, kpptemp2dMLD[:,   200], 'k-',label='GOTM Temperature at depth=%0.2f m'%GTMdepth[200].values); plt.legend(loc='best', fontsize='small')
# plt.plot(new_glidertime, new_gliderTempMLD[:,0], 'b-',label='New OSMOSIS Temperature at depth=%0.2f m'%GldrDepth[0].values); plt.legend(loc='best', fontsize='small')
# plt.xlabel('Time -Julian days'); plt.ylabel('Temperature -C')
# # plt.xlim(270,320)
# # plt.axis([250,275,14,18])
# plt.savefig('heatFlux_osmosisAnnual/compare_GOTM_new_GlidertempTimeseries_removedNAN.png');
# plt.show()
# exit()
# # TEMPERATURE-depth
# plt.plot(gliderTempMLD[4,:], GldrDepth[0:200], 'k-',label='Original Glider Temperature at time=248.767801 days'); plt.legend(loc='lower right', fontsize='small')
# plt.plot(new_gliderTempMLD[0,:], GldrDepth[0:200], 'b-',label='New Glider Temperature at time=248.767801 days'); plt.legend(loc='lower right', fontsize='small')
# plt.xlabel('Glider Temperature -C'); plt.ylabel('Glider depth -m')
# plt.savefig('heatFlux_osmosisAnnual/compare_original_new_GlidertempvsDepth_removedNAN.png');
# plt.show()

# INTEGRAL OF TEMPERATURE - timeseries
# # plt.plot(new_time, sumGliderTempMLD, 'r-',label='Original OSMOSIS Temperature integral'); plt.legend(loc='best', fontsize='small')
# plt.plot(new_glidertime, new_sumGliderTempMLD, 'b-',label='New OSMOSIS Temperature integral'); plt.legend(loc='best', fontsize='small')
# plt.plot(time2, sumkppTempMLD, 'k-',label='GOTM Temperature integral'); plt.legend(loc='best', fontsize='small')
# plt.xlabel('Time -Julian days'); plt.ylabel('Temperature Integral -C')
# # plt.axis([250,275,14,18])
# # plt.xlim(270,320)
# plt.savefig('heatFlux_osmosisAnnual/compare_GOTM_new_Glidertempintegral_Timeseries_removedNAN.png');
# plt.show()

# test1 = np.argwhere(~np.isnan(gliderTempMLD.values)); print('where are the NONnan values in glider MLD temp', test1)
# print(test1.shape)
# gliderTempMLD = gliderTempMLD.dropna('time'); print('list of glider temp mld without NaN', gliderTempMLD.shape)

# NaN_sumGliderTemp = sumGliderTempMLD[sumGliderTempMLD < 4000]; print('testing- OSMOSIS sum(temp) below 4000', NaN_sumGliderTemp[50:300])
# NaN_times = new_time[sumGliderTempMLD < 4000]; print('testing- OSMOSIS times of sum(temp) below 4000', NaN_times[50:300], NaN_times.shape)

# avgsumGliderTempMLD = (gliderTempMLD.sum(axis=1)*dz_OS)/dz_OS;# print('average sum of temp under Glider MLD', avgsumGliderTempMLD[0:200], avgsumGliderTempMLD.shape)
heatBudget_OS = []
midtime_OS = []

for s in range(0, 4095):
    CT_OS = (sumGliderTempMLD[s+1].values - sumGliderTempMLD[s].values)/dt_OS[s]
    diffT_OS = (new_time[s+1].values + new_time[s].values)/2
    heatBudget_OS.append(CT_OS)
    midtime_OS.append(diffT_OS)

heatBudget_OS = np.asarray(heatBudget_OS); midtime_OS = np.asarray(midtime_OS)
# print('OSMOSIS heat budget', heatBudget_OS.shape)

new_heatBudget_OS = []
new_midtime_OS = []

for x in range(0, 2931):
    nCT_OS = (new_sumGliderTempMLD[x+1].values - new_sumGliderTempMLD[x].values)/new_dt_OS[x]
    ndiffT_OS = (new_glidertime[x+1].values + new_glidertime[x].values)/2
    new_heatBudget_OS.append(nCT_OS)
    new_midtime_OS.append(ndiffT_OS)

new_heatBudget_OS = np.asarray(new_heatBudget_OS); new_midtime_OS = np.asarray(new_midtime_OS)
print('New OSMOSIS heat budget', new_heatBudget_OS, new_heatBudget_OS.shape)
Glidermidtime = midtime_OS[212:520]; GliderheatBudget = heatBudget_OS[212:520];
# GliderHB2d = np.c_[Glidermidtime,GliderheatBudget]; print('2D HB :',GliderHB2d[64:71])
# print('ORIGINAL Glider mid time', Glidermidtime);# print('ORIGINAL Glider heat budget', GliderheatBudget)
newGlidermidtime = new_midtime_OS[119:409]; newGliderheatBudget = new_heatBudget_OS[119:409];
# newGliderHB2d = np.c_[newGlidermidtime,newGliderheatBudget]; print('2D new HB :', newGliderHB2d[59:64])
# print('NEW Glider mid time', newGlidermidtime);# print('NEW Glider heat budget', newGliderheatBudget)


# print(Glidermidtime[0:5])
# print(newGlidermidtime[0:5])

### ----------------- PLOTS -compare original and new Glider heat budget ---------------------- ###

# -------------------------------------------- #
### annual comparison of heat energy budgets ###
# -------------------------------------------- #

GOTMHB_interpAnnual = np.interp(convertedDaysH, midtime, heatBudget)
newGliderHB_interpAnnual = np.interp(convertedDaysH, new_midtime_OS, new_heatBudget_OS)

fig, axs = plt.subplots(3,1)
fig.add_subplot(111, frameon=False) # add a big axis, hide frame

axs[0].plot(convertedDaysH, GOTMHB_interpAnnual, 'b-.', label = 'interpolated GOTM heat budget');axs[0].legend(loc='best', fontsize = 'small')
axs[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[0].set_xlim(260,640)
axs[0].set_ylim(-0.00025,0.00025)
axs[0].grid(True)

axs[1].plot(convertedDaysH, newGliderHB_interpAnnual, 'g-.', label = 'interpolated new Glider heat budget'); axs[1].legend(loc='best', fontsize = 'small')
axs[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[1].set_xlim(260,640)
axs[1].set_ylim(-0.00025,0.00025)
axs[1].grid(True)

axs[2].plot(convertedDaysH, fullsum_radHeatflux, 'k-.', label = 'ERA-interim heat budget'); axs[2].legend(loc='best', fontsize = 'small')
axs[2].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[2].set_xlim(260,640)
axs[2].set_ylim(-0.00025,0.00025)
axs[2].grid(True)

plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False) # hide tick@Jand tick label of the big axis
plt.xlabel('Time -Julian days');
plt.ylabel('Estimated ocean surface heat budget - W/m2')
fig.tight_layout()
plt.savefig('heatFlux_osmosisAnnual/compare_interpolatedAnnualHeatBudget.png'); plt.show()
exit()
# plt.plot(midtime, heatBudget, 'r-.', label ='GOTM heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.plot(new_midtime_OS, new_heatBudget_OS, 'k-.', label ='new Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(convertedDaysH, GOTMHB_interpAnnual, 'b-.', label = 'interpolated GOTM heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(convertedDaysH, newGliderHB_interpAnnual, 'g-.', label = 'interpolated new Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(convertedDaysH, fullsum_radHeatflux, 'k-.', label = 'ERA-interim heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.plot(np.linspace(start=265, stop=300, num=160), NaN_times[0], '.r', label = 'Glider HB error due to missing values'); plt.legend(loc='best', fontsize = 'small')
# plt.yscale('log')
plt.xlabel('Time -Julian days'); plt.ylabel('Estimated ocean surface heat budget - W/m2')
# plt.axis([260,650,-0.005,0.005])
# plt.ylim(-0.0005,0.0005)
plt.savefig('heatFlux_osmosisAnnual/test__interpolatedAnnualHeatBudget.png'); plt.show()
# new Glider heat budget timeseries
# plt.plot(new_midtime_OS, new_heatBudget_OS, 'b-')#, label ='New Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.xlabel('Time -Julian days'); plt.ylabel('Glider ocean surface heat budget (New) - W/m2')
# plt.axis([270,275,-0.000075,0.00005])
# plt.savefig('heatFlux_osmosisAnnual/new_GliderheatBudget_removedNAN.png'); plt.show()
# # new Glider heat budget timeseries (Autumn)
# plt.plot(newGlidermidtime, newGliderheatBudget, 'b-')#, label ='New Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.xlabel('Time -Julian days'); plt.ylabel('Glider ocean surface heat budget (New) - W/m2')
# # plt.axis([270,275,-0.000075,0.00005])
# plt.savefig('heatFlux_osmosisAnnual/new_GliderheatBudget_Autumn_removedNAN.png'); plt.show()

# timeseries
plt.plot(Glidermidtime, GliderheatBudget, 'k-', label ='Original Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(newGlidermidtime, newGliderheatBudget, 'b-', label ='New Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.plot(datfileTime, GliderHB_interpolated, 'b.', label = 'interpolated Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.plot(np.linspace(start=265, stop=300, num=160), NaN_times[0], '.r', label = 'Glider HB error due to missing values'); plt.legend(loc='best', fontsize = 'small')
plt.xlabel('Time -Julian days'); plt.ylabel('Glider ocean surface heat budget - W/m2')
plt.axis([265,300,-0.00015,0.00015])
# plt.axis([270,275,-0.000075,0.00005])
plt.savefig('heatFlux_osmosisAnnual/compare_original_new_GliderheatBudget_removedNAN.png'); plt.show()

# # profile
# plt.plot(GliderheatBudget[2,:], GldrDepth[0:200], 'k-', label ='Original Glider heat budget at time= 267.999838 days'); plt.legend(loc='upper right', fontsize = 'small')
# plt.plot(newGliderheatBudget[1,:], GldrDepth[0:200], 'b-', label ='New Glider heat budget at time= 267.999838 days'); plt.legend(loc='upper right', fontsize = 'small')
# # plt.plot(datfileTime, GliderHB_interpolated, 'b.', label = 'interpolated Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
# # plt.plot(np.linspace(start=265, stop=300, num=160), NaN_times[0], '.r', label = 'Glider HB error due to missing values'); plt.legend(loc='best', fontsize = 'small')
# plt.xlabel('Glider ocean surface heat budget - W/m2'); plt.ylabel('Depth -m')
# # plt.axis([270,275,-0.000075,0.00005])
# plt.savefig('heatFlux_osmosisAnnual/compare_original_new_GliderheatBudgetvsdepth_removedNAN.png'); plt.show()



# print('GOTM time:', convertedDaysH[0:240])
# print(midtime[0:4463]) #,midtime[52594], midtime.shape)
# print(heatBudget[0],heatBudget1[52594], midtime.shape)
# print('Glider midtime', midtime_OS[204:522])
# avgsumGliderTempMLD2d = np.vstack((new_time, avgsumGliderTempMLD)).T; 
# print('2d average sum of temp under Glider MLD', avgsumGliderTempMLD2d[277:295])
lowHBv1 = GldrTemp[277,:]; #print(lowHBv1.values)
lowHBv2 = GldrTemp[295,:]; #print(lowHBv2.values)

# # Glider averaged heat budget within MLD
# plt.plot(new_time[204:522], avgsumGliderTempMLD[204:522], 'b-', label = 'int_mld^0(Glider_theta)'); plt.legend(loc='lower right') 
# plt.xlabel('Time -Julian days'); plt.ylabel('Glider averaged temperature for MLD')
# # plt.axis([265,300,-0.02,0.03])
# # plt.axis([250,650,-0.0002,0.0003])
# plt.savefig('averaged_temperature_Glider_mixed_layer_depth.png'); plt.show()

### ----------------------------- calculate high and lowla_t ----------------------------- ###

Aut1HighLat = np.full(shape=1465, fill_value=277.0, dtype=np.float)
Aut2HighLat = np.full(shape=1465, fill_value=283.625, dtype=np.float)
Aut1LowLat = np.full(shape=1465, fill_value=301.875, dtype=np.float)
Aut2LowLat = np.full(shape=1465, fill_value=310.5, dtype=np.float)

Spr1HighLat = np.full(shape=1465, fill_value=524.0, dtype=np.float)
Spr2HighLat = np.full(shape=1465, fill_value=531.75, dtype=np.float)
Spr1LowLat = np.full(shape=1465, fill_value=493.0, dtype=np.float)
Spr2LowLat = np.full(shape=1465, fill_value=501.75, dtype=np.float)

############################################################
### -------- reduce GOTM heat budget for Glider -------- ###
# This is used to find the correlation between heat budget #
############################################################
print('OSMOSIS midtime', midtime_OS.shape)
print('OSMOSIS heat budget', heatBudget_OS.shape)
timeGliderHB = np.c_[midtime_OS, heatBudget_OS]
print('2d Glider time and heatbudget',timeGliderHB[0:10])
timeGTMHB = np.c_[midtime, heatBudget]
print('2d GOTM time and heatbudget', timeGTMHB[0:10])

fmt = '%1.10f', '%1.10f'
np.savetxt('heatFlux_osmosisAnnual/2dGlider_HB.txt', timeGliderHB, fmt=fmt)

def closestGlidHBArray(timeGliderHB, timeGTMHB, elem = 0):
    def find_nearest(value):
        return min(timeGTMHB, key=lambda x:abs(x[elem]-value))
    return [find_nearest(timeGliderHB[i][elem]) for i in range(len(timeGliderHB))]

newGOTMHBArray = closestGlidHBArray(timeGliderHB,timeGTMHB); newGOTMHBArray = np.asarray(newGOTMHBArray)
# print('reduced GOTM heat budget for Glider :', newGOTMHBArray[0:10])
np.savetxt('heatFlux_osmosisAnnual/reducedGOTMHB.txt', newGOTMHBArray,fmt=fmt)
###################################################################
### -------- reduce GOTM heat budget for heatflux data -------- ###
# This allows to perform the correlation between the heat budgets #
###################################################################

datfle_timeOSHB = np.vstack((convertedDaysH,fullsum_radHeatflux)).T
np.savetxt('heatFlux_osmosisAnnual/2dDAT_HB.txt', datfle_timeOSHB, fmt=fmt)
# print('2d heatflux.dat time and heatbudget', datfle_timeOSHB[0:10])

def closestheatdatArray(datfle_timeOSHB,timeGTMHB, elem = 0):
    def find_nearestt(valuee):
        return min(timeGTMHB, key=lambda x:abs(x[elem]-valuee))
    return [find_nearestt(datfle_timeOSHB[i][elem]) for i in range(len(datfle_timeOSHB))]

new2GOTMHBArray = closestheatdatArray(datfle_timeOSHB,timeGTMHB); new2GOTMHBArray = np.asarray(new2GOTMHBArray)
# print('reduced GOTM heat budget for heatflux.dat :', new2GOTMHBArray[0:10])
np.savetxt('heatFlux_osmosisAnnual/reduced2GOTMHB.txt', new2GOTMHBArray,fmt=fmt)


### ----------------- INTERPOLATION ---------------------- ###
# print('Interpolation information :')
# print('heat flux time', convertedDaysH[7:247])
# print('GOTM time', midtime[0:4319])
# print('Glider time', midtime_OS[212:520])

datfileTime = convertedDaysH[7:247]; datfileHBudget = fullsum_radHeatflux[7:247]

#GOTM#
GOTMmidtime = midtime[0:4319]; GOTMheatBudget = heatBudget[0:4319]; 
GOTMHB_interpolated = np.interp(datfileTime, GOTMmidtime, GOTMheatBudget)

plt.plot(GOTMmidtime, GOTMheatBudget, 'g-.', label = 'GOTM heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.plot(datfileTime, GOTMHB_interpolated, 'b.', label = 'interpolated GOTM heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(newGlidermidtime, newGliderheatBudget, 'b-.', label ='New Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(datfileTime, datfileHBudget, 'k-.', label = 'ERA-interim heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.axis([265,300,-0.00015,0.00015])
plt.xlabel('Time -Julian days'); plt.ylabel('Ocean surface heat budget - W/m2')
plt.savefig('heatFlux_osmosisAnnual/test3__interpolateGOTMheatbudget.png'); plt.show()

#Glider#
# Glidermidtime = midtime_OS[212:520]; GliderheatBudget = heatBudget_OS[212:520];
GliderHB_interpolated = np.interp(datfileTime, Glidermidtime, GliderheatBudget)
newGliderHB_interpolated = np.interp(datfileTime, newGlidermidtime, newGliderheatBudget)


# plt.plot(Glidermidtime, GliderheatBudget, 'k.-', label ='Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(datfileTime, GliderHB_interpolated, 'b-.', label = 'interpolated Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
plt.plot(datfileTime, newGliderHB_interpolated, 'g-.', label = 'interpolated new Glider heat budget'); plt.legend(loc='best', fontsize = 'small')
# plt.plot(np.linspace(start=265, stop=300, num=160), NaN_times[0], '.r', label = 'Glider HB error due to missing values'); plt.legend(loc='best', fontsize = 'small')
plt.xlabel('Time -Julian days'); plt.ylabel('Ocean surface heat budget - W/m2')
# plt.axis([260,650,-0.005,0.005])
plt.ylim(-0.0005,0.0005)
plt.savefig('heatFlux_osmosisAnnual/test__interpolateGliderheatbudget.png'); plt.show()


exit()

# ### --------------------- Correlation -------------------------- ###

# print('length of interpolated GOTM heat budget', len(GOTMHB_interpolated))
# print('length of interpolated Glider heat budget', len(GliderHB_interpolated))

# correlation_GOTM_fluxFile = sum(np.multiply(datfileHBudget, GOTMHB_interpolated))/np.sqrt(np.sum(np.square(datfileHBudget))*np.sum(np.square(GOTMHB_interpolated)))
# print('(Autumn 2012) correlation coefficient between GOTM and heatflux.dat heat budget :', correlation_GOTM_fluxFile)

# correlation_Glider_fluxfile = sum(np.multiply(datfileHBudget, GliderHB_interpolated))/np.sqrt(np.sum(np.square(datfileHBudget))*np.sum(np.square(GliderHB_interpolated)))
# print('(Autumn 2012) correlation coefficient between OSMOSIS and heatflux.dat heat budget :', correlation_Glider_fluxfile)

# correlation_GOTM_Glider = sum(np.multiply(GOTMHB_interpolated,GliderHB_interpolated))/np.sqrt(np.sum(np.square(GOTMHB_interpolated))*np.sum(np.square(GliderHB_interpolated)))
# print('correlation coefficient between GOTM and OSMOSIS heat budget :', correlation_GOTM_Glider)
# exit()


### --------------------- PLOTS -------------------------- ###

## [H(0) + lambda*E(0) + R_n(0)]/rho*c_p

# plt.plot(convertedDaysH[0:247],fullsum_radHeatflux[0:247], 'g-',label = '[H + R_n(0)]/rho*c_p'); plt.legend(loc='best')
plt.plot(convertedDaysH,fullsum_radHeatflux, '--k',label = '[H(0) + R_n(0)]/rho*c_p'); plt.legend(loc=(0.25,0.7))
plt.xlabel('Time -Julian days'); plt.ylabel('Sum {solar radiation & surface heat flux} -W/m2');# plt.title('Surface Heat Flux - OSMOSIS Autumn observations')
# plot the ranges for high and low ranges of La_t in Spring and Autumn
plt.plot(Aut1HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'salmon',linestyle='dotted',label = 'Autumn high La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Aut2HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'salmon',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75),
plt.plot(Aut1LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'slategrey',linestyle='dotted',label = 'Autumn low La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Aut2LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'slategrey',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75))
plt.plot(Spr1HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'peru',linestyle='dotted',label = 'Spring high La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Spr2HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'peru',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75))
plt.plot(Spr1LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'cadetblue',linestyle='dotted',label = 'Spring low La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Spr2LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'cadetblue',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75))
# plt.axis([265,300,-0.00015,0.0002])
# plt.axis([445,539,-400,100])
plt.savefig('heatFlux_osmosisAnnual/GOTM_fullsumSurface_heat_fluxAnnual.png'); plt.show()

# GOTM heat budget within MLD
# plt.plot(midtime[0:4463], heatBudget[0:4463], 'k-', label = 'd/dt int_mld^0(GOTM_theta)'); plt.legend(loc='best') 
plt.plot(midtime, heatBudget, '--k', label = 'd/dt int_mld^0(GOTM_theta)'); plt.legend(loc=(0.25,0.7))
# plot the ranges for high and low ranges of La_t in Spring and Autumn
plt.plot(Aut1HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'salmon',linestyle='dotted',label = 'Autumn high La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Aut2HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'salmon',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75),
plt.plot(Aut1LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'slategrey',linestyle='dotted',label = 'Autumn low La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Aut2LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'slategrey',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75))
plt.plot(Spr1HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'peru',linestyle='dotted',label = 'Spring high La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Spr2HighLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'peru',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75))
plt.plot(Spr1LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'cadetblue',linestyle='dotted',label = 'Spring low La_t range'); plt.legend(loc=(0.25,0.7), fontsize = 'small')
plt.plot(Spr2LowLat, np.linspace(start=-0.00015, stop=0.0002, num=1465), color = 'cadetblue',linestyle='dotted')#,label = 'La_t = 0.3'); plt.legend(loc=(0.05,0.75))
plt.xlabel('Time -Julian days'); plt.ylabel('GOTM surface heat budget for mixed layer depth')
# plt.axis([265,300,-0.00015,0.0002])
plt.savefig('heatFlux_osmosisAnnual/derived_heat_budget_mixed_layer_depthAnnual.png'); plt.show()


# Glider heat budget within MLD
plt.plot(midtime_OS[204:522], heatBudget_OS[204:522], 'b-', label = 'd/dt int_mld^0(Glider_theta)'); plt.legend(loc='lower right') 
plt.xlabel('Time -Julian days'); plt.ylabel('Glider heat budget for mixed layer depth')
# plt.axis([265,300,-0.02,0.03])
plt.axis([265,300,-0.01,0.02])
plt.savefig('heatFlux_osmosisAnnual/derived_heat_budget_Glider_mixed_layer_depth.png'); plt.show()

## compare surface heat budgets
# plt.plot(midtime, heatBudget, 'k-', label = 'd/dt int_mld^0(GOTM_theta)'); plt.legend(loc='lower right')
plt.plot(midtime[0:4463], heatBudget[0:4463], 'k-', label = 'd/dt int_mld^0(GOTM_theta)'); plt.legend(loc='lower right')
# plt.plot(convertedDaysH,fullsum_radHeatflux, 'g-',label = '[H + R_n(0)]/rho*c_p'); plt.legend(loc='lower right')
plt.plot(convertedDaysH[0:247],fullsum_radHeatflux[0:247], 'g-',label = '[H + R_n(0)]/rho*c_p'); plt.legend(loc='lower right')
# plt.plot(midtime_OS[204:522], heatBudget_OS[204:522], 'b-.', label = 'd/dt int_mld^0(Glider_theta)'); plt.legend(loc='lower right') 
plt.xlabel('Time -Julian days'); plt.ylabel('Surface heat budget -W/m2');# plt.title('Surface Heat Flux - OSMOSIS Autumn observations')
# plt.axis([265,300,-0.02,0.03])
plt.savefig('heatFlux_osmosisAnnual/heatBudgetMLD_GOTM_Glider_solarRadheatflux_comparison30days.png'); plt.show()

exit()


#############################################################
#### --------- rolling correlation coefficient --------- ####
#############################################################
## manual calculation of correlation coefficient
x_bar = np.mean(fullsum_radHeatflux); print(x_bar)
y_bar = np.mean(heatBudget); print(y_bar)

## correlation between surface heat budgets

df = pd.DataFrame(dict(x=heatBudget))

CORR_VALS = np.array(fullsum_radHeatflux)
def get_correlation(vals):
    return pearsonr(vals, CORR_VALS)[0]

df['correlation'] = df.rolling(len(CORR_VALS)).apply(get_correlation)
# f = open("Pearsoncorrelation_GOTM_OSMOSIS_heatbudget.txt", "a")
# print('Pearson correlation applied to GOTM and OSMOSIS heat budgets', df['correlation'].values, file=f)
# f.close()
print(df['correlation'].values)
exit()
print(df)
plt.scatter(np.arange(len(heatBudget)),df['correlation'].values,c='k',marker='.',s=2);
plt.xlabel('length of GOTM heat budget'); plt.ylabel('Correlation coefficient between GOTM & OSMOSIS heat budget');# plt.title('Surface Heat Flux - OSMOSIS Autumn observations')
plt.savefig('CorrCoeff_heatBudget_GOTM_OSMOSIS.png'); plt.show()

exit()



## H(0) + lambda*E(0) + R_n(0)

plt.plot(convertedDaysH,sum_radHeatflux,'k-.',label = 'H + R_n(0)'); plt.legend(loc=(0.65,0.05), fontsize = 'small')
plt.xlabel('Time -Julian days'); plt.ylabel('Sum of solar radiation and surface heat flux -W/m2');# plt.title('Surface Heat Flux - OSMOSIS Autumn observations')
# plt.axis([445,539,-400,100])
plt.axis([260,650,-0.0002,0.0003])
plt.savefig('GOTM_sumSurface_heat_flux.png'); plt.show()

exit()

### --------------------- Average values for low and high La_t  -------------------------- ###

AnnualHeatBudavg = np.c_[convertedDaysH, fullsum_radHeatflux];# print(AnnualHeatBudavg)
AnnualHeatBud = np.c_[convertedDaysH, sum_radHeatflux];# print(AnnualHeatBud)
DateAHB = np.c_[hdate,AnnualHeatBud];# print(DateAHB)


## Autumn - low Lat
# DaysAutlowLat = np.where((convertedDaysH <= 310.5) & (convertedDaysH >= 301.875)); 
# RowAutlowlat = convertedDaysH[DaysAutlowLat]
# print(DaysAutlowLat)
# print(RowAutlowlat)
# AnnHBAutlowLatavg = AnnualHeatBudavg[278:348,:]; print('heat budget (avg) values for Aut low Lat',AnnHBAutlowLatavg)
# AnnHBAutlowLatavgMean = np.mean(AnnHBAutlowLatavg[:,1]); print('mean HB(avg) for Aut low Lat',AnnHBAutlowLatavgMean)

# AnnHBAutlowLat = AnnualHeatBud[278:348,:]; print('heat budget values for Aut low Lat',AnnHBAutlowLat)
# AnnHBAutlowLatMean = np.mean(AnnHBAutlowLat[:,1]); print('mean HB for Aut low Lat',AnnHBAutlowLatMean)
# print(DateAHB[278:348])


## Autumn - high Lat
# DaysAutlowLat = np.where((convertedDaysH <= 283.625) & (convertedDaysH >= 277.0))
# RowAutlowlat = convertedDaysH[DaysAutlowLat]
# print('days with high Lat in Autumn', DaysAutlowLat)
# print('row# with high Lat in Autumn', RowAutlowlat)

# AnnHBAutlowLatavg = AnnualHeatBudavg[80:133,:]; print('heat budget (avg) values for Aut high Lat',AnnHBAutlowLatavg)
# AnnHBAutlowLatavgMean = np.mean(AnnHBAutlowLatavg[:,1]); print('mean HB(avg) for Aut high Lat',AnnHBAutlowLatavgMean)

# AnnHBAutlowLat = AnnualHeatBud[80:133,:]; print('heat budget values for Aut high Lat',AnnHBAutlowLat)
# AnnHBAutlowLatMean = np.mean(AnnHBAutlowLat[:,1]); print('mean HB for Aut high Lat',AnnHBAutlowLatMean)
# print(DateAHB[80:133])

## Spring - low Lat
# DaysAutlowLat = np.where((convertedDaysH <= 503.75) & (convertedDaysH >= 493.0))
# RowAutlowlat = convertedDaysH[DaysAutlowLat]
# print('days with low Lat in Spring', DaysAutlowLat)
# print('row# with low Lat in Spring', RowAutlowlat)

# AnnHBAutlowLatavg = AnnualHeatBudavg[1808:1894,:]; print('heat budget (avg) values for Spr low Lat',AnnHBAutlowLatavg)
# AnnHBAutlowLatavgMean = np.mean(AnnHBAutlowLatavg[:,1]); print('mean HB(avg) for Spr low Lat',AnnHBAutlowLatavgMean)

# AnnHBAutlowLat = AnnualHeatBud[1808:1894,:]; print('heat budget values for Spr low Lat',AnnHBAutlowLat)
# AnnHBAutlowLatMean = np.mean(AnnHBAutlowLat[:,1]); print('mean HB for Spr low Lat',AnnHBAutlowLatMean)
# print(DateAHB[1808:1894])

## Spring - high Lat
DaysAutlowLat = np.where((convertedDaysH <= 531.75) & (convertedDaysH >= 524.0))
RowAutlowlat = convertedDaysH[DaysAutlowLat]
print('days with high Lat in Spring', DaysAutlowLat)
print('row# with high Lat in Spring', RowAutlowlat)

AnnHBAutlowLatavg = AnnualHeatBudavg[2056:2118,:]; print('heat budget (avg) values for Spr high Lat',AnnHBAutlowLatavg)
AnnHBAutlowLatavgMean = np.mean(AnnHBAutlowLatavg[:,1]); print('mean HB(avg) for Spr high Lat',AnnHBAutlowLatavgMean)

AnnHBAutlowLat = AnnualHeatBud[2056:2118,:]; print('heat budget values for Spr high Lat',AnnHBAutlowLat)
AnnHBAutlowLatMean = np.mean(AnnHBAutlowLat[:,1]); print('mean HB for Spr high Lat',AnnHBAutlowLatMean)
print(DateAHB[2056:2118])
