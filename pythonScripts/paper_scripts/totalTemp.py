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
plt.switch_backend('agg')
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

gotmkpp = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_placebo1.nc", decode_times=False)
# gotmkeps = xr.open_dataset("osmosis_annual_surface_forcingkeps/OSMOSIS_glider_comparison.nc", decode_times=False)
# changed the definition of eddy viscosity (num(k)) in kpp.F90
gotmkppLa = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_nu_Lat.nc",decode_times=False)
# gotmkppLaRib = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib.nc",decode_times=False)
# gotmkppLaRib_g2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_g2.nc",decode_times=False)
# gotmkppLaRib_g025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_g025.nc",decode_times=False)
# gotmkppLaRib_g05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_g05.nc",decode_times=False)
# gotmkppLaRib_g075 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_g075.nc",decode_times=False)
# gotmkppLaRib_g125 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_g125.nc",decode_times=False)
# gotmkppLaRib_a2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_a2.nc",decode_times=False)
# gotmkppLaRib_a025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_a025.nc",decode_times=False)
# gotmkppLaRib_a05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_a05.nc",decode_times=False)
# gotmkppLaRib_a075 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_a075.nc",decode_times=False)
# gotmkppLaRib_a125 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nu_Lat_Rib_a125.nc",decode_times=False)
# changed the definition of heat diffusitivity (nuh(k)) in kpp.F90
# gotmkpp_heat = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_KPPHD_PAP.nc",decode_times=False)
# gotmkpp_heat125dlta = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_KPPHD_8*dlt_PAP.nc",decode_times=False)
# changed the definition of heat diffusitivity (nuh(k)) and eddy viscosity (num(k)) in kpp.F90 ## gamma=1
# gotmkpp_m_h = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/osmosis_annual_surface_forcingKPP/OSMOSIS_new_Visc_Diff.nc",decode_times=False)
# multiply factors to originial num and nuh
# gotmkpp_025_025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_025_025_PAP.nc", decode_times=False)
# gotmkpp_025_05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_025_05_PAP.nc", decode_times=False)
# gotmkpp_025_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_025_1_PAP.nc", decode_times=False)
# gotmkpp_025_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_025_2_PAP.nc", decode_times=False)
# gotmkpp_025_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_025_4_PAP.nc", decode_times=False)
# gotmkpp_05_025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_05_025_PAP.nc", decode_times=False)
# gotmkpp_05_05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_05_05_PAP.nc", decode_times=False)
# gotmkpp_05_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_05_1_PAP.nc", decode_times=False)
# gotmkpp_05_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_05_2_PAP.nc", decode_times=False)
# gotmkpp_05_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_05_4_PAP.nc", decode_times=False)
# gotmkpp_1_025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_1_025_PAP.nc", decode_times=False)
# gotmkpp_1_05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_1_05_PAP.nc", decode_times=False)
# gotmkpp_1_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_placebo2.nc", decode_times=False)
gotmkpp_1_2 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_2Kh.nc", decode_times=False)
# gotmkpp_R_1_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nuh2.nc", decode_times=False)
# gotmkpp_RO_1_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nuh2_Ric2.nc", decode_times=False)
# gotmkpp_RR_1_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_nuh2_Ric05.nc", decode_times=False)
# gotmkpp_1_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_1_4_PAP.nc", decode_times=False)
# gotmkpp_HD_1_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_KPPHD_1_4_PAP.nc", decode_times=False)
# gotmkpp_2_025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_2_025_PAP.nc", decode_times=False)
# gotmkpp_2_05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_2_05_PAP.nc", decode_times=False)
gotmkpp_2_1 = xr.open_dataset("~/Documents/gotm-4.0.0-edit/simulations/annual_edited/OSMOSIS_2Km.nc", decode_times=False)
# gotmkpp_R_2_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_num2.nc", decode_times=False)
# gotmkpp_RO_2_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_num2_Ric05.nc", decode_times=False)
# gotmkpp_RR_2_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_num2_Ric2.nc", decode_times=False)
# gotmkpp_RR_2_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_05multi_RiG.nc", decode_times=False)
# gotmkpp_2_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_2_2_PAP.nc", decode_times=False)
# gotmkpp_2_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_2_4_PAP.nc", decode_times=False)
# gotmkpp_4_025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_4_025_PAP.nc", decode_times=False)
# gotmkpp_4_05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_4_05_PAP.nc", decode_times=False)
# gotmkpp_4_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_4_1_PAP.nc", decode_times=False)
# gotmkpp_R_4_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_num4.nc", decode_times=False)
# gotmkpp_4_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_4_2_PAP.nc", decode_times=False)
# gotmkpp_4_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_4_4_PAP.nc", decode_times=False)
# gotmkpp_58_025 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_58_025_PAP.nc", decode_times=False)
# gotmkpp_58_05 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_58_05_PAP.nc", decode_times=False)
# gotmkpp_58_1 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_58_1_PAP.nc", decode_times=False)
# gotmkpp_58_2 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_58_2_PAP.nc", decode_times=False)
# gotmkpp_58_4 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/spring/OSMOSIS_numnuh_58_4_PAP.nc", decode_times=False)

# changed the definition of heat diffusitivity (nuh(k)) in kpp.F90 ## gamma=0.5
# gotmkpp_heatg5 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/osmosis_annual_surface_forcingKPP/OSMOSIS_newHeatDiff_gam0.5.nc",decode_times=False)
# changed the definition of heat diffusitivity (nuh(k)) in kpp.F90 ## gamma=1.5
# gotmkpp_heatg15 = xr.open_dataset("../../../gotm-4.0.0-edit/simulations/osmosis_annual_surface_forcingKPP/OSMOSIS_newHeatDiff_gam1.5.nc",decode_times=False)
# print(gotmkpp_heat125dlta.variables.keys())

GTMtemp = gotmkpp.temp; GTMdepth = gotmkpp.z; GTMnum = gotmkpp.num; GTMnuh = gotmkpp.nuh; GTMsal = gotmkpp.salt; GTMo2 = gotmkpp.o2_obs; GTMu = gotmkpp.u; GTMv = gotmkpp.v; GTMxflx = gotmkpp.gamu; GTMyflx = gotmkpp.gamv; time = gotmkpp.time
# GTMketemp = gotmkeps.temp; GTMkedepth = gotmkeps.z; GTMkenum = gotmkeps.num; GTMkenuh = gotmkeps.nuh; GTMkeu = gotmkeps.u; GTMkev = gotmkeps.v; ketime = gotmkeps.time
GTMtempLa = gotmkppLa.temp; GTMdepthLa = gotmkppLa.z; GTMnumLa = gotmkppLa.num; GTMnuhLa = gotmkppLa.nuh; GTMuLa = gotmkppLa.u; GTMvLa = gotmkppLa.v; timeLa = gotmkppLa.time
# GTMtempLaRib = gotmkppLaRib.temp; GTMdepthLaRib = gotmkppLaRib.z; GTMnumLaRib = gotmkppLaRib.num; GTMnuhLaRib = gotmkppLaRib.nuh; GTMuLaRib = gotmkppLaRib.u; GTMvLaRib = gotmkppLaRib.v; timeLaRib = gotmkppLaRib.time
# GTMtempLaRiba2 = gotmkppLaRib_a2.temp; GTMdepthLaRiba2 = gotmkppLaRib_a2.z; GTMnumLaRiba2 = gotmkppLaRib_a2.num; GTMnuhLaRiba2 = gotmkppLaRib_a2.nuh; GTMuLaRiba2 = gotmkppLaRib_a2.u; GTMvLaRiba2 = gotmkppLaRib_a2.v; timeLaRiba2 = gotmkppLaRib_a2.time
# GTMtempLaRiba025 = gotmkppLaRib_a025.temp; GTMdepthLaRiba025 = gotmkppLaRib_a025.z; GTMnumLaRiba025 = gotmkppLaRib_a025.num; GTMnuhLaRiba025 = gotmkppLaRib_a025.nuh; GTMuLaRiba025 = gotmkppLaRib_a025.u; GTMvLaRiba025 = gotmkppLaRib_a025.v; timeLaRiba025 = gotmkppLaRib_a025.time
# GTMtempLaRiba05 = gotmkppLaRib_a05.temp; GTMdepthLaRiba05 = gotmkppLaRib_a05.z; GTMnumLaRiba05 = gotmkppLaRib_a05.num; GTMnuhLaRiba05 = gotmkppLaRib_a05.nuh; GTMuLaRiba05 = gotmkppLaRib_a05.u; GTMvLaRiba05 = gotmkppLaRib_a05.v; timeLaRiba05 = gotmkppLaRib_a05.time
# GTMtempLaRiba075 = gotmkppLaRib_a075.temp; GTMdepthLaRiba075 = gotmkppLaRib_a075.z; GTMnumLaRiba075 = gotmkppLaRib_a075.num; GTMnuhLaRiba075 = gotmkppLaRib_a075.nuh; GTMuLaRiba075 = gotmkppLaRib_a075.u; GTMvLaRiba075 = gotmkppLaRib_a075.v; timeLaRiba075 = gotmkppLaRib_a075.time
# GTMtempLaRiba125 = gotmkppLaRib_a125.temp; GTMdepthLaRiba125 = gotmkppLaRib_a125.z; GTMnumLaRiba125 = gotmkppLaRib_a125.num; GTMnuhLaRiba125 = gotmkppLaRib_a125.nuh; GTMuLaRiba125 = gotmkppLaRib_a125.u; GTMvLaRiba125 = gotmkppLaRib_a125.v; timeLaRiba125 = gotmkppLaRib_a125.time
# GTMtempLaRibg2 = gotmkppLaRib_g2.temp; GTMdepthLaRibg2 = gotmkppLaRib_g2.z; GTMnumLaRibg2 = gotmkppLaRib_g2.num; GTMnuhLaRibg2 = gotmkppLaRib_g2.nuh; GTMuLaRibg2 = gotmkppLaRib_g2.u; GTMvLaRibg2 = gotmkppLaRib_g2.v; timeLaRibg2 = gotmkppLaRib_g2.time
# GTMtempLaRibg025 = gotmkppLaRib_g025.temp; GTMdepthLaRibg025 = gotmkppLaRib_g025.z; GTMnumLaRibg025 = gotmkppLaRib_g025.num; GTMnuhLaRibg025 = gotmkppLaRib_g025.nuh; GTMuLaRibg025 = gotmkppLaRib_g025.u; GTMvLaRibg025 = gotmkppLaRib_g025.v; timeLaRibg025 = gotmkppLaRib_g025.time
# GTMtempLaRibg05 = gotmkppLaRib_g05.temp; GTMdepthLaRibg05 = gotmkppLaRib_g05.z; GTMnumLaRibg05 = gotmkppLaRib_g05.num; GTMnuhLaRibg05 = gotmkppLaRib_g05.nuh; GTMuLaRibg05 = gotmkppLaRib_g05.u; GTMvLaRibg05 = gotmkppLaRib_g05.v; timeLaRibg05 = gotmkppLaRib_g05.time
# GTMtempLaRibg075 = gotmkppLaRib_g075.temp; GTMdepthLaRibg075 = gotmkppLaRib_g075.z; GTMnumLaRibg075 = gotmkppLaRib_g075.num; GTMnuhLaRibg075 = gotmkppLaRib_g075.nuh; GTMuLaRibg075 = gotmkppLaRib_g075.u; GTMvLaRibg075 = gotmkppLaRib_g075.v; timeLaRibg075 = gotmkppLaRib_g075.time
# GTMtempLaRibg125 = gotmkppLaRib_g125.temp; GTMdepthLaRibg125 = gotmkppLaRib_g125.z; GTMnumLaRibg125 = gotmkppLaRib_g125.num; GTMnuhLaRibg125 = gotmkppLaRib_g125.nuh; GTMuLaRibg125 = gotmkppLaRib_g125.u; GTMvLaRibg125 = gotmkppLaRib_g125.v; timeLaRibg125 = gotmkppLaRib_g125.time
# GTMtemp_heat = gotmkpp_heat.temp; GTMdepth_heat = gotmkpp_heat.z; GTMnum_heat = gotmkpp_heat.num; GTMnuh_heat = gotmkpp_heat.nuh; GTMu_heat = gotmkpp_heat.u; GTMv_heat = gotmkpp_heat.v; time_heat = gotmkpp_heat.time
# GTMtemp_heat125dlta = gotmkpp_heat125dlta.temp; GTMdepth_heat125dlta = gotmkpp_heat125dlta.z; GTMnum_heat125dlta = gotmkpp_heat125dlta.num; GTMnuh_heat125dlta = gotmkpp_heat125dlta.nuh; GTMu_heat125dlta = gotmkpp_heat125dlta.u; GTMv_heat125dlta = gotmkpp_heat125dlta.v; time_heat125dlta = gotmkpp_heat125dlta.time
# GTMtemp_mh = gotmkpp_m_h.temp; GTMdepth_mh = gotmkpp_m_h.z; GTMnum_mh = gotmkpp_m_h.num; GTMnuh_mh = gotmkpp_m_h.nuh; GTMu_mh = gotmkpp_m_h.u; GTMv_mh = gotmkpp_m_h.v; time_mh = gotmkpp_m_h.time
# GTMtemp_heatg5 = gotmkpp_heatg5.temp; GTMdepth_heatg5 = gotmkpp_heatg5.z; GTMnum_heatg5 = gotmkpp_heatg5.num; GTMnuh_heatg5 = gotmkpp_heatg5.nuh; time_heatg5 = gotmkpp_heatg5.time
# GTMtemp_heatg15 = gotmkpp_heatg15.temp; GTMdepth_heatg15 = gotmkpp_heatg15.z; GTMnum_heatg15 = gotmkpp_heatg15.num; GTMnuh_heatg15 = gotmkpp_heatg15.nuh; time_heatg15 = gotmkpp_heatg15.time
# GTMtemp025025 = gotmkpp_025_025.temp; GTMdepth025025 = gotmkpp_025_025.z; GTMnum025025 = gotmkpp_025_025.num; GTMnuh025025 = gotmkpp_025_025.nuh; GTMu025025 = gotmkpp_025_025.u; GTMv025025 = gotmkpp_025_025.v; time025025 = gotmkpp_025_025.time
# GTMtemp02505 = gotmkpp_025_05.temp; GTMdepth02505 = gotmkpp_025_05.z; GTMnum02505 = gotmkpp_025_05.num; GTMnuh02505 = gotmkpp_025_05.nuh; GTMu02505 = gotmkpp_025_05.u; GTMv02505 = gotmkpp_025_05.v; time02505 = gotmkpp_025_05.time
# GTMtemp0251 = gotmkpp_025_1.temp; GTMdepth0251 = gotmkpp_025_1.z; GTMnum0251 = gotmkpp_025_1.num; GTMnuh0251 = gotmkpp_025_1.nuh; GTMu0251 = gotmkpp_025_1.u; GTMv0251 = gotmkpp_025_1.v; time0251 = gotmkpp_025_1.time
# GTMtemp0252 = gotmkpp_025_2.temp; GTMdepth0252 = gotmkpp_025_2.z; GTMnum0252 = gotmkpp_025_2.num; GTMnuh0252 = gotmkpp_025_2.nuh; GTMu0252 = gotmkpp_025_2.u; GTMv0252 = gotmkpp_025_2.v; time0252 = gotmkpp_025_2.time
# GTMtemp0254 = gotmkpp_025_4.temp; GTMdepth0254 = gotmkpp_025_4.z; GTMnum0254 = gotmkpp_025_4.num; GTMnuh0254 = gotmkpp_025_4.nuh; GTMu0254 = gotmkpp_025_4.u; GTMv0254 = gotmkpp_025_4.v; time0254 = gotmkpp_025_4.time
# GTMtemp05025 = gotmkpp_05_025.temp; GTMdepth05025 = gotmkpp_05_025.z; GTMnum05025 = gotmkpp_05_025.num; GTMnuh05025 = gotmkpp_05_025.nuh; GTMu05025 = gotmkpp_05_025.u; GTMv05025 = gotmkpp_05_025.v; time05025 = gotmkpp_05_025.time
# GTMtemp0505 = gotmkpp_05_05.temp; GTMdepth0505 = gotmkpp_05_05.z; GTMnum0505 = gotmkpp_05_05.num; GTMnuh0505 = gotmkpp_05_05.nuh; GTMu0505 = gotmkpp_05_05.u; GTMv0505 = gotmkpp_05_05.v; time0505 = gotmkpp_05_05.time
# GTMtemp051 = gotmkpp_05_1.temp; GTMdepth051 = gotmkpp_05_1.z; GTMnum051 = gotmkpp_05_1.num; GTMnuh051 = gotmkpp_05_1.nuh; GTMu051 = gotmkpp_05_1.u; GTMv051 = gotmkpp_05_1.v; time051 = gotmkpp_05_1.time
# GTMtemp052 = gotmkpp_05_2.temp; GTMdepth052 = gotmkpp_05_2.z; GTMnum052 = gotmkpp_05_2.num; GTMnuh052 = gotmkpp_05_2.nuh; GTMu052 = gotmkpp_05_2.u; GTMv052 = gotmkpp_05_2.v; time052 = gotmkpp_05_2.time
# GTMtemp054 = gotmkpp_05_4.temp; GTMdepth054 = gotmkpp_05_4.z; GTMnum054 = gotmkpp_05_4.num; GTMnuh054 = gotmkpp_05_4.nuh; GTMu054 = gotmkpp_05_4.u; GTMv054 = gotmkpp_05_4.v; time054 = gotmkpp_05_4.time
# GTMtemp1025 = gotmkpp_1_025.temp; GTMdepth1025 = gotmkpp_1_025.z; GTMnum1025 = gotmkpp_1_025.num; GTMnuh1025 = gotmkpp_1_025.nuh; GTMu1025 = gotmkpp_1_025.u; GTMv1025 = gotmkpp_1_025.v; time1025 = gotmkpp_1_025.time
# GTMtemp105 = gotmkpp_1_05.temp; GTMdepth105 = gotmkpp_1_05.z; GTMnum105 = gotmkpp_1_05.num; GTMnuh105 = gotmkpp_1_05.nuh; GTMu105 = gotmkpp_1_05.u; GTMv105 = gotmkpp_1_05.v; time105 = gotmkpp_1_05.time
# GTMtemp11 = gotmkpp_1_1.temp; GTMdepth11 = gotmkpp_1_1.z; GTMnum11 = gotmkpp_1_1.num; GTMnuh11 = gotmkpp_1_1.nuh; GTMu11 = gotmkpp_1_1.u; GTMv11 = gotmkpp_1_1.v; time11 = gotmkpp_1_1.time
GTMtemp12 = gotmkpp_1_2.temp; GTMdepth12 = gotmkpp_1_2.z; GTMnum12 = gotmkpp_1_2.num; GTMnuh12 = gotmkpp_1_2.nuh; GTMu12 = gotmkpp_1_2.u; GTMv12 = gotmkpp_1_2.v; time12 = gotmkpp_1_2.time
# GTMtempR12 = gotmkpp_R_1_2.temp; GTMdepthR12 = gotmkpp_R_1_2.z; GTMnumR12 = gotmkpp_R_1_2.num; GTMnuhR12 = gotmkpp_R_1_2.nuh; GTMuR12 = gotmkpp_R_1_2.u; GTMvR12 = gotmkpp_R_1_2.v; timeR12 = gotmkpp_R_1_2.time
# GTMtempRO12 = gotmkpp_RO_1_2.temp; GTMdepthRO12 = gotmkpp_RO_1_2.z; GTMnumRO12 = gotmkpp_RO_1_2.num; GTMnuhRO12 = gotmkpp_RO_1_2.nuh; GTMuRO12 = gotmkpp_RO_1_2.u; GTMvRO12 = gotmkpp_RO_1_2.v; timeRO12 = gotmkpp_RO_1_2.time
# GTMtempRR12 = gotmkpp_RR_1_2.temp; GTMdepthRR12 = gotmkpp_RR_1_2.z; GTMnumRR12 = gotmkpp_RR_1_2.num; GTMnuhRR12 = gotmkpp_RR_1_2.nuh; GTMuRR12 = gotmkpp_RR_1_2.u; GTMvRR12 = gotmkpp_RR_1_2.v; timeRR12 = gotmkpp_RR_1_2.time
# GTMtempHD13 = gotmkpp_HD_1_3.temp; GTMdepthHD13 = gotmkpp_HD_1_3.z; GTMnumHD13 = gotmkpp_HD_1_3.num; GTMnuhHD13 = gotmkpp_HD_1_3.nuh; GTMuHD13 = gotmkpp_HD_1_3.u; GTMvHD13 = gotmkpp_HD_1_3.v; timeHD13 = gotmkpp_HD_1_3.time
# GTMtemp14 = gotmkpp_1_4.temp; GTMdepth14 = gotmkpp_1_4.z; GTMnum14 = gotmkpp_1_4.num; GTMnuh14 = gotmkpp_1_4.nuh; GTMu14 = gotmkpp_1_4.u; GTMv14 = gotmkpp_1_4.v; time14 = gotmkpp_1_4.time
# GTMtempHD14 = gotmkpp_HD_1_4.temp; GTMdepthHD14 = gotmkpp_HD_1_4.z; GTMnumHD14 = gotmkpp_HD_1_4.num; GTMnuhHD14 = gotmkpp_HD_1_4.nuh; GTMuHD14 = gotmkpp_HD_1_4.u; GTMvHD14 = gotmkpp_HD_1_4.v; timeHD14 = gotmkpp_HD_1_4.time
# GTMtemp2025 = gotmkpp_2_025.temp; GTMdepth2025 = gotmkpp_2_025.z; GTMnum2025 = gotmkpp_2_025.num; GTMnuh2025 = gotmkpp_2_025.nuh; GTMu2025 = gotmkpp_2_025.u; GTMv2025 = gotmkpp_2_025.v; time2025 = gotmkpp_2_025.time
# GTMtemp205 = gotmkpp_2_05.temp; GTMdepth205 = gotmkpp_2_05.z; GTMnum205 = gotmkpp_2_05.num; GTMnuh205 = gotmkpp_2_05.nuh; GTMu205 = gotmkpp_2_05.u; GTMv205 = gotmkpp_2_05.v; time205 = gotmkpp_2_05.time
GTMtemp21 = gotmkpp_2_1.temp; GTMdepth21 = gotmkpp_2_1.z; GTMnum21 = gotmkpp_2_1.num; GTMnuh21 = gotmkpp_2_1.nuh; GTMu21 = gotmkpp_2_1.u; GTMv21 = gotmkpp_2_1.v; time21 = gotmkpp_2_1.time
# GTMtemp22 = gotmkpp_2_2.temp; GTMdepth22 = gotmkpp_2_2.z; GTMnum22 = gotmkpp_2_2.num; GTMnuh22 = gotmkpp_2_2.nuh; GTMu22 = gotmkpp_2_2.u; GTMv22 = gotmkpp_2_2.v; time22 = gotmkpp_2_2.time
# GTMtemp24 = gotmkpp_2_4.temp; GTMdepth24 = gotmkpp_2_4.z; GTMnum24 = gotmkpp_2_4.num; GTMnuh24 = gotmkpp_2_4.nuh; GTMu24 = gotmkpp_2_4.u; GTMv24 = gotmkpp_2_4.v; time24 = gotmkpp_2_4.time
# GTMtempR21 = gotmkpp_R_2_1.temp; GTMdepthR21 = gotmkpp_R_2_1.z; GTMnumR21 = gotmkpp_R_2_1.num; GTMnuhR21 = gotmkpp_R_2_1.nuh; GTMuR21 = gotmkpp_R_2_1.u; GTMvR21 = gotmkpp_R_2_1.v; timeR21 = gotmkpp_R_2_1.time
# GTMtempRO21 = gotmkpp_RO_2_1.temp; GTMdepthRO21 = gotmkpp_RO_2_1.z; GTMnumRO21 = gotmkpp_RO_2_1.num; GTMnuhRO21 = gotmkpp_RO_2_1.nuh; GTMuRO21 = gotmkpp_RO_2_1.u; GTMvRO21 = gotmkpp_RO_2_1.v; timeRO21 = gotmkpp_RO_2_1.time
# GTMtempRR21 = gotmkpp_RR_2_1.temp; GTMdepthRR21 = gotmkpp_RR_2_1.z; GTMnumRR21 = gotmkpp_RR_2_1.num; GTMnuhRR21 = gotmkpp_RR_2_1.nuh; GTMuRR21 = gotmkpp_RR_2_1.u; GTMvRR21 = gotmkpp_RR_2_1.v; timeRR21 = gotmkpp_RR_2_1.time
# GTMtemp4025 = gotmkpp_4_025.temp; GTMdepth4025 = gotmkpp_4_025.z; GTMnum4025 = gotmkpp_4_025.num; GTMnuh4025 = gotmkpp_4_025.nuh; GTMu4025 = gotmkpp_4_025.u; GTMv4025 = gotmkpp_4_025.v; time4025 = gotmkpp_4_025.time
# GTMtemp405 = gotmkpp_4_05.temp; GTMdepth405 = gotmkpp_4_05.z; GTMnum405 = gotmkpp_4_05.num; GTMnuh405 = gotmkpp_4_05.nuh; GTMu405 = gotmkpp_4_05.u; GTMv405 = gotmkpp_4_05.v; time405 = gotmkpp_4_05.time
# GTMtemp41 = gotmkpp_4_1.temp; GTMdepth41 = gotmkpp_4_1.z; GTMnum41 = gotmkpp_4_1.num; GTMnuh41 = gotmkpp_4_1.nuh; GTMu41 = gotmkpp_4_1.u; GTMv41 = gotmkpp_4_1.v; time41 = gotmkpp_4_1.time
# GTMtempR41 = gotmkpp_R_4_1.temp; GTMdepthR41 = gotmkpp_R_4_1.z; GTMnumR41 = gotmkpp_R_4_1.num; GTMnuhR41 = gotmkpp_R_4_1.nuh; GTMuR41 = gotmkpp_R_4_1.u; GTMvR41 = gotmkpp_R_4_1.v; timeR41 = gotmkpp_R_4_1.time
# GTMtemp42 = gotmkpp_4_2.temp; GTMdepth42 = gotmkpp_4_2.z; GTMnum42 = gotmkpp_4_2.num; GTMnuh42 = gotmkpp_4_2.nuh; GTMu42 = gotmkpp_4_2.u; GTMv42 = gotmkpp_4_2.v; time42 = gotmkpp_4_2.time
# GTMtemp44 = gotmkpp_4_4.temp; GTMdepth44 = gotmkpp_4_4.z; GTMnum44 = gotmkpp_4_4.num; GTMnuh44 = gotmkpp_4_4.nuh; GTMu44 = gotmkpp_4_4.u; GTMv44 = gotmkpp_4_4.v; time44 = gotmkpp_4_4.time
# GTMtemp58025 = gotmkpp_58_025.temp; GTMdepth58025 = gotmkpp_58_025.z; GTMnum58025 = gotmkpp_58_025.num; GTMnuh58025 = gotmkpp_58_025.nuh; GTMu58025 = gotmkpp_58_025.u; GTMv58025 = gotmkpp_58_025.v; time58025 = gotmkpp_58_025.time
# GTMtemp5805 = gotmkpp_58_05.temp; GTMdepth5805 = gotmkpp_58_05.z; GTMnum5805 = gotmkpp_58_05.num; GTMnuh5805 = gotmkpp_58_05.nuh; GTMu5805 = gotmkpp_58_05.u; GTMv5805 = gotmkpp_58_05.v; time5805 = gotmkpp_58_05.time
# GTMtemp581 = gotmkpp_58_1.temp; GTMdepth581 = gotmkpp_58_1.z; GTMnum581 = gotmkpp_58_1.num; GTMnuh581 = gotmkpp_58_1.nuh; GTMu581 = gotmkpp_58_1.u; GTMv581 = gotmkpp_58_1.v; time581 = gotmkpp_58_1.time
# GTMtemp582 = gotmkpp_58_2.temp; GTMdepth582 = gotmkpp_58_2.z; GTMnum582 = gotmkpp_58_2.num; GTMnuh582 = gotmkpp_58_2.nuh; GTMu582 = gotmkpp_58_2.u; GTMv582 = gotmkpp_58_2.v; time582 = gotmkpp_58_2.time
# GTMtemp584 = gotmkpp_58_4.temp; GTMdepth584 = gotmkpp_58_4.z; GTMnum584 = gotmkpp_58_4.num; GTMnuh584 = gotmkpp_58_4.nuh; GTMu584 = gotmkpp_58_4.u; GTMv584 = gotmkpp_58_4.v; time584 = gotmkpp_58_4.time

## Convert time
time2 = 267.75 + (time[:]/86400);# ketime2 = 267.75 + (ketime[:]/86400);
time2La = 267.75 + (timeLa[:]/86400); 
# time2LaRib = 445.0 + (timeLaRib[:]/86400); 
# time2LaRiba2 = 445.0 + (timeLaRiba2[:]/86400); 
# time2LaRiba025 = 445.0 + (timeLaRiba025[:]/86400); 
# time2LaRiba05 = 445.0 + (timeLaRiba05[:]/86400); 
# time2LaRiba075 = 445.0 + (timeLaRiba075[:]/86400); 
# time2LaRiba125 = 445.0 + (timeLaRiba125[:]/86400); 
# time2LaRibg2 = 445.0 + (timeLaRibg2[:]/86400); 
# time2LaRibg025 = 445.0 + (timeLaRibg025[:]/86400); 
# time2LaRibg05 = 445.0 + (timeLaRibg05[:]/86400); 
# time2LaRibg075 = 445.0 + (timeLaRibg075[:]/86400); 
# time2LaRibg125 = 445.0 + (timeLaRibg125[:]/86400); 
# time2_heat = 445.0 + (time_heat[:]/86400); time2_heat125dlta = 445.0 + (time_heat125dlta[:]/86400);
# time2_mh = 445.0 + (time_mh[:]/86400); 
# time2_heatg5 = 445.0 + (time_heatg5[:]/86400); time2_heatg15 = 445.0 + (time_heatg15[:]/86400); 
# time_025025 = 445.0 + (time02505[:]/86400);
# time_02505 = 445.0 + (time02505[:]/86400);
# time_0251 = 445.0 + (time0251[:]/86400);
# time_0252 = 445.0 + (time0252[:]/86400);
# time_0254 = 445.0 + (time0254[:]/86400);

# time_05025 = 445.0 + (time05025[:]/86400);
# time_0505 = 445.0 + (time0505[:]/86400);
# time_051 = 445.0 + (time051[:]/86400);
# time_052 = 445.0 + (time052[:]/86400);
# time_054 = 445.0 + (time054[:]/86400);

# time_1025 = 445.0 + (time1025[:]/86400);
# time_105 = 445.0 + (time105[:]/86400);
# time_11 = 445.0 + (time11[:]/86400);
time_12 = 267.75 + (time12[:]/86400);
# # time_HD13 = 445.0 + (timeHD13[:]/86400);
# time_14 = 445.0 + (time14[:]/86400);
# # time_HD14 = 445.0 + (timeHD14[:]/86400);
# time_R12 = 445.0 + (timeR12[:]/86400);
# time_RO12 = 445.0 + (timeRO12[:]/86400);
# time_RR12 = 445.0 + (timeRR12[:]/86400);

# time_2025 = 445.0 + (time2025[:]/86400);
# time_205 = 445.0 + (time205[:]/86400);
time_21 = 267.75 + (time21[:]/86400);
# time_22 = 445.0 + (time22[:]/86400);
# time_24 = 445.0 + (time24[:]/86400);
# time_R21 = 445.0 + (timeR21[:]/86400);
# time_RO21 = 445.0 + (timeRO21[:]/86400);
# time_RR21 = 445.0 + (timeRR21[:]/86400);

# time_4025 = 445.0 + (time4025[:]/86400);
# time_405 = 445.0 + (time405[:]/86400);
# time_41 = 445.0 + (time41[:]/86400);
# time_R41 = 445.0 + (timeR41[:]/86400);
# time_42 = 445.0 + (time42[:]/86400);
# time_44 = 445.0 + (time44[:]/86400);

# time_58025 = 445.0 + (time58025[:]/86400);
# time_5805 = 445.0 + (time5805[:]/86400);
# time_581 = 445.0 + (time581[:]/86400);
# time_582 = 445.0 + (time582[:]/86400);
# time_584 = 445.0 + (time584[:]/86400);

# times = np.c_[time,time2];
# print('compare times', times)
# print('GOTM times :', time2.values);

gotmkpp = gotmkpp.assign(time2=("time", time2)); gotmkpp = gotmkpp.swap_dims({"time" : "time2"})
gotmkppLa = gotmkppLa.assign(time2La=("time", time2La)); gotmkppLa = gotmkppLa.swap_dims({"time" : "time2La"})
# gotmkppLaRib = gotmkppLaRib.assign(time2LaRib=("time", time2LaRib)); gotmkppLaRib = gotmkppLaRib.swap_dims({"time" : "time2LaRib"})
# gotmkppLaRib_a2 = gotmkppLaRib_a2.assign(time2LaRiba2=("time", time2LaRiba2)); gotmkppLaRib_a2 = gotmkppLaRib_a2.swap_dims({"time" : "time2LaRiba2"})
# gotmkppLaRib_a025 = gotmkppLaRib_a025.assign(time2LaRiba025=("time", time2LaRiba025)); gotmkppLaRib_a025 = gotmkppLaRib_a025.swap_dims({"time" : "time2LaRiba025"})
# gotmkppLaRib_a05 = gotmkppLaRib_a05.assign(time2LaRiba05=("time", time2LaRiba05)); gotmkppLaRib_a05 = gotmkppLaRib_a05.swap_dims({"time" : "time2LaRiba05"})
# gotmkppLaRib_a075 = gotmkppLaRib_a075.assign(time2LaRiba075=("time", time2LaRiba075)); gotmkppLaRib_a075 = gotmkppLaRib_a075.swap_dims({"time" : "time2LaRiba075"})
# gotmkppLaRib_a125 = gotmkppLaRib_a125.assign(time2LaRiba125=("time", time2LaRiba125)); gotmkppLaRib_a15 = gotmkppLaRib_a125.swap_dims({"time" : "time2LaRiba125"})
# gotmkppLaRib_g2 = gotmkppLaRib_g2.assign(time2LaRibg2=("time", time2LaRibg2)); gotmkppLaRib_g2 = gotmkppLaRib_g2.swap_dims({"time" : "time2LaRibg2"})
# gotmkppLaRib_g025 = gotmkppLaRib_g025.assign(time2LaRibg025=("time", time2LaRibg025)); gotmkppLaRib_g025 = gotmkppLaRib_g025.swap_dims({"time" : "time2LaRibg025"})
# gotmkppLaRib_g05 = gotmkppLaRib_g05.assign(time2LaRibg05=("time", time2LaRibg05)); gotmkppLaRib_g05 = gotmkppLaRib_g05.swap_dims({"time" : "time2LaRibg05"})
# gotmkppLaRib_g075 = gotmkppLaRib_g075.assign(time2LaRibg075=("time", time2LaRibg075)); gotmkppLaRib_g075 = gotmkppLaRib_g075.swap_dims({"time" : "time2LaRibg075"})
# gotmkppLaRib_g125 = gotmkppLaRib_g125.assign(time2LaRibg125=("time", time2LaRibg125)); gotmkppLaRib_g125 = gotmkppLaRib_g125.swap_dims({"time" : "time2LaRibg125"})
# gotmkpp_heat = gotmkpp_heat.assign(time2_heat=("time", time2_heat)); gotmkpp_heat = gotmkpp_heat.swap_dims({"time" : "time2_heat"})
# gotmkpp_heat125dlta = gotmkpp_heat125dlta.assign(time2_heat125dlta=("time", time2_heat125dlta)); gotmkpp_heat125dlta = gotmkpp_heat125dlta.swap_dims({"time" : "time2_heat125dlta"})
# gotmkpp_m_h = gotmkpp_m_h.assign(time2_mh=("time", time2_mh)); gotmkpp_m_h = gotmkpp_m_h.swap_dims({"time" : "time2_mh"})
# gotmkpp_heatg5 = gotmkpp_heatg5.assign(time2_heatg5=("time", time2_heatg5)); gotmkpp_heatg5 = gotmkpp_heatg5.swap_dims({"time" : "time2_heatg5"})
# gotmkpp_heatg15 = gotmkpp_heatg15.assign(time2_heatg15=("time", time2_heatg15)); gotmkpp_heatg15 = gotmkpp_heatg15.swap_dims({"time" : "time2_heatg15"})
# gotmkeps = gotmkeps.assign(ketime2=("time", ketime2)); gotmkeLaps = gotmkeps.swap_dims({"time" : "ketime2"})
# gotmkpp_025_025 = gotmkpp_025_025.assign(time_025025=("time", time_025025)); gotmkpp_025_025 = gotmkpp_025_025.swap_dims({"time" : "time_025025"})
# gotmkpp_025_05 = gotmkpp_025_05.assign(time_02505=("time", time_02505)); gotmkpp_025_05 = gotmkpp_025_05.swap_dims({"time" : "time_02505"})
# gotmkpp_025_1 = gotmkpp_025_1.assign(time_0251=("time", time_0251)); gotmkpp_025_1 = gotmkpp_025_1.swap_dims({"time" : "time_0251"})
# gotmkpp_025_2 = gotmkpp_025_2.assign(time_0252=("time", time_0252)); gotmkpp_025_2 = gotmkpp_025_2.swap_dims({"time" : "time_0252"})
# gotmkpp_025_4 = gotmkpp_025_4.assign(time_0254=("time", time_0254)); gotmkpp_025_4 = gotmkpp_025_4.swap_dims({"time" : "time_0254"})
# gotmkpp_05_025 = gotmkpp_05_025.assign(time_05025=("time", time_05025)); gotmkpp_05_025 = gotmkpp_05_025.swap_dims({"time" : "time_05025"})
# gotmkpp_05_05 = gotmkpp_05_05.assign(time_0505=("time", time_0505)); gotmkpp_05_05 = gotmkpp_05_05.swap_dims({"time" : "time_0505"})
# gotmkpp_05_1 = gotmkpp_05_1.assign(time_051=("time", time_051)); gotmkpp_05_1 = gotmkpp_05_1.swap_dims({"time" : "time_051"})
# gotmkpp_05_2 = gotmkpp_05_2.assign(time_052=("time", time_052)); gotmkpp_05_2 = gotmkpp_05_2.swap_dims({"time" : "time_052"})
# gotmkpp_05_4 = gotmkpp_05_4.assign(time_054=("time", time_054)); gotmkpp_05_4 = gotmkpp_05_4.swap_dims({"time" : "time_054"})
# gotmkpp_1_025 = gotmkpp_1_025.assign(time_1025=("time", time_1025)); gotmkpp_1_025 = gotmkpp_1_025.swap_dims({"time" : "time_1025"})
# gotmkpp_1_05 = gotmkpp_1_05.assign(time_105=("time", time_105)); gotmkpp_1_05 = gotmkpp_1_05.swap_dims({"time" : "time_105"})
# gotmkpp_1_1 = gotmkpp_1_1.assign(time_11=("time", time_11)); gotmkpp_1_1 = gotmkpp_1_1.swap_dims({"time" : "time_11"})
gotmkpp_1_2 = gotmkpp_1_2.assign(time_12=("time", time_12)); gotmkpp_1_2 = gotmkpp_1_2.swap_dims({"time" : "time_12"})
# gotmkpp_R_1_2 = gotmkpp_R_1_2.assign(time_R12=("time", time_R12)); gotmkpp_R_1_2 = gotmkpp_R_1_2.swap_dims({"time" : "time_R12"})
# gotmkpp_RO_1_2 = gotmkpp_RO_1_2.assign(time_RO12=("time", time_RO12)); gotmkpp_RO_1_2 = gotmkpp_RO_1_2.swap_dims({"time" : "time_RO12"})
# gotmkpp_RR_1_2 = gotmkpp_RR_1_2.assign(time_RR12=("time", time_RR12)); gotmkpp_RR_1_2 = gotmkpp_RR_1_2.swap_dims({"time" : "time_RR12"})
# gotmkpp_HD_1_3 = gotmkpp_HD_1_3.assign(time_HD13=("time", time_HD13)); gotmkpp_HD_1_3 = gotmkpp_HD_1_3.swap_dims({"time" : "time_HD13"})
# gotmkpp_1_4 = gotmkpp_1_4.assign(time_14=("time", time_14)); gotmkpp_1_4 = gotmkpp_1_4.swap_dims({"time" : "time_14"})
# gotmkpp_HD_1_4 = gotmkpp_HD_1_4.assign(time_HD14=("time", time_HD14)); gotmkpp_HD_1_4 = gotmkpp_HD_1_4.swap_dims({"time" : "time_HD14"})
# gotmkpp_2_025 = gotmkpp_2_025.assign(time_2025=("time", time_2025)); gotmkpp_2_025 = gotmkpp_2_025.swap_dims({"time" : "time_2025"})
# gotmkpp_2_05 = gotmkpp_2_05.assign(time_205=("time", time_205)); gotmkpp_2_05 = gotmkpp_2_05.swap_dims({"time" : "time_205"})
gotmkpp_2_1 = gotmkpp_2_1.assign(time_21=("time", time_21)); gotmkpp_2_1 = gotmkpp_2_1.swap_dims({"time" : "time_21"})
# gotmkpp_2_2 = gotmkpp_2_2.assign(time_22=("time", time_22)); gotmkpp_2_2 = gotmkpp_2_2.swap_dims({"time" : "time_22"})
# gotmkpp_2_4 = gotmkpp_2_4.assign(time_24=("time", time_24)); gotmkpp_2_4 = gotmkpp_2_4.swap_dims({"time" : "time_24"})
# gotmkpp_RO_2_1 = gotmkpp_RO_2_1.assign(time_RO21=("time", time_RO21)); gotmkpp_RO_2_1 = gotmkpp_RO_2_1.swap_dims({"time" : "time_RO21"})
# gotmkpp_RR_2_1 = gotmkpp_RR_2_1.assign(time_RR21=("time", time_RR21)); gotmkpp_RR_2_1 = gotmkpp_RR_2_1.swap_dims({"time" : "time_RR21"})
# gotmkpp_4_025 = gotmkpp_4_025.assign(time_4025=("time", time_4025)); gotmkpp_4_025 = gotmkpp_4_025.swap_dims({"time" : "time_4025"})
# gotmkpp_4_05 = gotmkpp_4_05.assign(time_405=("time", time_405)); gotmkpp_4_05 = gotmkpp_4_05.swap_dims({"time" : "time_405"})
# gotmkpp_4_1 = gotmkpp_4_1.assign(time_41=("time", time_41)); gotmkpp_4_1 = gotmkpp_4_1.swap_dims({"time" : "time_41"})
# gotmkpp_R_4_1 = gotmkpp_R_4_1.assign(time_R41=("time", time_R41)); gotmkpp_R_4_1 = gotmkpp_R_4_1.swap_dims({"time" : "time_R41"})
# gotmkpp_4_2 = gotmkpp_4_2.assign(time_42=("time", time_42)); gotmkpp_4_2 = gotmkpp_4_2.swap_dims({"time" : "time_42"})
# gotmkpp_4_4 = gotmkpp_4_4.assign(time_44=("time", time_44)); gotmkpp_4_4 = gotmkpp_4_4.swap_dims({"time" : "time_44"})
# gotmkpp_58_025 = gotmkpp_58_025.assign(time_58025=("time", time_58025)); gotmkpp_58_025 = gotmkpp_58_025.swap_dims({"time" : "time_58025"})
# gotmkpp_58_05 = gotmkpp_58_05.assign(time_5805=("time", time_5805)); gotmkpp_58_05 = gotmkpp_58_05.swap_dims({"time" : "time_5805"})
# gotmkpp_58_1 = gotmkpp_58_1.assign(time_581=("time", time_581)); gotmkpp_58_1 = gotmkpp_58_1.swap_dims({"time" : "time_581"})
# gotmkpp_58_2 = gotmkpp_58_2.assign(time_582=("time", time_582)); gotmkpp_58_2 = gotmkpp_58_2.swap_dims({"time" : "time_582"})
# gotmkpp_58_4 = gotmkpp_58_4.assign(time_584=("time", time_584)); gotmkpp_58_4 = gotmkpp_58_4.swap_dims({"time" : "time_584"})

print(gotmkpp.variables.items())

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

# plt.plot(time2,mldKPP); plt.xlabel('time -julian day'); plt.ylabel('mixed layer depth -m')
# plt.axis([260,341,-140,0]); plt.title('GOTM Mixed Layer Depth')
# plt.savefig('GOTM_Mixed_Layer_Depth_test6.png'); plt.show()

## using KPP EV model
mld_tempKPPLa = GTMtempLa.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtempLa-mld_tempKPPLa).argmin('z'); # print(z_indexes)
mldKPPLa = GTMtempLa.z[z_indexes]; 
mldKPPLa = mldKPPLa.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib original
# mld_tempKPPLaRib = GTMtempLaRib.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRib-mld_tempKPPLaRib).argmin('z'); # print(z_indexes)
# mldKPPLaRib = GTMtempLaRib.z[z_indexes]; 
# mldKPPLaRib = mldKPPLaRib.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=1;Ch=2
# mld_tempKPPLaRiba2 = GTMtempLaRiba2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRiba2-mld_tempKPPLaRiba2).argmin('z'); # print(z_indexes)
# mldKPPLaRiba2 = GTMtempLaRiba2.z[z_indexes]; 
# mldKPPLaRiba2 = mldKPPLaRiba2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=1;Ch=0.5
# mld_tempKPPLaRiba05 = GTMtempLaRiba05.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRiba05-mld_tempKPPLaRiba05).argmin('z'); # print(z_indexes)
# mldKPPLaRiba05 = GTMtempLaRiba05.z[z_indexes]; 
# mldKPPLaRiba05 = mldKPPLaRiba05.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=1;Ch=0.25
# mld_tempKPPLaRiba025 = GTMtempLaRiba025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRiba025-mld_tempKPPLaRiba025).argmin('z'); # print(z_indexes)
# mldKPPLaRiba025 = GTMtempLaRiba025.z[z_indexes]; 
# mldKPPLaRiba025 = mldKPPLaRiba025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=1;Ch=0.75
# mld_tempKPPLaRiba075 = GTMtempLaRiba075.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRiba075-mld_tempKPPLaRiba075).argmin('z'); # print(z_indexes)
# mldKPPLaRiba075 = GTMtempLaRiba075.z[z_indexes]; 
# mldKPPLaRiba075 = mldKPPLaRiba075.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=1;Ch=1.25
# mld_tempKPPLaRiba125 = GTMtempLaRiba125.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRiba125-mld_tempKPPLaRiba125).argmin('z'); # print(z_indexes)
# mldKPPLaRiba125 = GTMtempLaRiba125.z[z_indexes]; 
# mldKPPLaRiba125 = mldKPPLaRiba125.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=2;Ch=1
# mld_tempKPPLaRibg2 = GTMtempLaRibg2.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg2-mld_tempKPPLaRibg2).argmin('z'); # print(z_indexes)
# mldKPPLaRibg2 = GTMtempLaRibg2.z[z_indexes]; 
# mldKPPLaRibg2 = mldKPPLaRibg2.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=0.25;Ch=1
# mld_tempKPPLaRibg025 = GTMtempLaRibg025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg025-mld_tempKPPLaRibg025).argmin('z'); # print(z_indexes)
# mldKPPLaRibg025 = GTMtempLaRibg025.z[z_indexes]; 
# mldKPPLaRibg025 = mldKPPLaRibg025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=0.5;Ch=1
# mld_tempKPPLaRibg05 = GTMtempLaRibg05.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg05-mld_tempKPPLaRibg05).argmin('z'); # print(z_indexes)
# mldKPPLaRibg05 = GTMtempLaRibg05.z[z_indexes]; 
# mldKPPLaRibg05 = mldKPPLaRibg05.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=0.75;Ch=1
# mld_tempKPPLaRibg075 = GTMtempLaRibg075.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg075-mld_tempKPPLaRibg075).argmin('z'); # print(z_indexes)
# mldKPPLaRibg075 = GTMtempLaRibg075.z[z_indexes]; 
# mldKPPLaRibg075 = mldKPPLaRibg075.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model - use Lat and Rib Cm=1.25;Ch=1
# mld_tempKPPLaRibg125 = GTMtempLaRibg125.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempLaRibg125-mld_tempKPPLaRibg125).argmin('z'); # print(z_indexes)
# mldKPPLaRibg125 = GTMtempLaRibg125.z[z_indexes]; 
# mldKPPLaRibg125 = mldKPPLaRibg125.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP HD model
# mld_tempKPPheat = GTMtemp_heat.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp_heat-mld_tempKPPheat).argmin('z'); # print(z_indexes)
# mldKPPheat = GTMtemp_heat.z[z_indexes]; 
# mldKPPheat = mldKPPheat.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP HD 8*dlta model
# mld_tempKPPheat125dlta = GTMtemp_heat125dlta.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp_heat125dlta-mld_tempKPPheat125dlta).argmin('z'); # print(z_indexes)
# mldKPPheat125dlta = GTMtemp_heat125dlta.z[z_indexes]; 
# mldKPPheat125dlta = mldKPPheat125dlta.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using k-epsilon model
# mld_tempKeps = GTMketemp.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMketemp-mld_tempKeps).argmin('z'); # print(z_indexes)
# mldKeps = GTMketemp.z[z_indexes]; 
# mldKeps = mldKeps.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ###  num == 0.25
# ## using KPP model (num=0.25; nuh=0.25)
# mld_tempKPP025025 = GTMtemp025025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp025025-mld_tempKPP025025).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP025025 = GTMtemp025025.z[z_indexes]; 
# mldKPP025025 = mldKPP025025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.25; nuh=0.5)
# mld_tempKPP02505 = GTMtemp02505.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp02505-mld_tempKPP02505).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP02505 = GTMtemp02505.z[z_indexes]; 
# mldKPP02505 = mldKPP02505.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.25; nuh=1)
# mld_tempKPP0251 = GTMtemp0251.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp0251-mld_tempKPP0251).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP0251 = GTMtemp0251.z[z_indexes]; 
# mldKPP0251 = mldKPP0251.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.25; nuh=2)
# mld_tempKPP0252 = GTMtemp0252.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp0252-mld_tempKPP0252).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP0252 = GTMtemp0252.z[z_indexes]; 
# mldKPP0252 = mldKPP0252.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.25; nuh=4)
# mld_tempKPP0254 = GTMtemp0254.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp0254-mld_tempKPP0254).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP0254 = GTMtemp0254.z[z_indexes]; 
# mldKPP0254 = mldKPP0254.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ###  num == 0.5
# ## using KPP model (num=0.5; nuh=0.25)
# mld_tempKPP05025 = GTMtemp05025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp05025-mld_tempKPP05025).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP05025 = GTMtemp05025.z[z_indexes]; 
# mldKPP05025 = mldKPP05025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.5; nuh=0.5)
# mld_tempKPP0505 = GTMtemp0505.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp0505-mld_tempKPP0505).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP0505 = GTMtemp0505.z[z_indexes]; 
# mldKPP0505 = mldKPP0505.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.5; nuh=1)
# mld_tempKPP051 = GTMtemp051.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp051-mld_tempKPP051).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP051 = GTMtemp051.z[z_indexes]; 
# mldKPP051 = mldKPP051.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.5; nuh=2)
# mld_tempKPP052 = GTMtemp052.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp052-mld_tempKPP052).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP052 = GTMtemp0252.z[z_indexes]; 
# mldKPP052 = mldKPP052.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=0.5; nuh=4)
# mld_tempKPP054 = GTMtemp054.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp054-mld_tempKPP054).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP054 = GTMtemp054.z[z_indexes]; 
# mldKPP054 = mldKPP054.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# ###  num == 1
# ## using KPP model (num=1; nuh=0.25)
# mld_tempKPP1025 = GTMtemp1025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp1025-mld_tempKPP1025).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP1025 = GTMtemp1025.z[z_indexes]; 
# mldKPP1025 = mldKPP1025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=1; nuh=0.5)
# mld_tempKPP105 = GTMtemp105.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp105-mld_tempKPP105).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP105 = GTMtemp105.z[z_indexes]; 
# mldKPP105 = mldKPP105.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=1; nuh=1)
# mld_tempKPP11 = GTMtemp11.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp11-mld_tempKPP11).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP11 = GTMtemp11.z[z_indexes]; 
# mldKPP11 = mldKPP11.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model (num=1; nuh=2)
mld_tempKPP12 = GTMtemp12.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtemp12-mld_tempKPP12).argmin('z'); # print(z_indexes)
print('GOTM temp MLD indexes', z_indexes)
mldKPP12 = GTMtemp12.z[z_indexes]; 
mldKPP12 = mldKPP12.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP HD model (num=1; nuh=3)
# mld_tempKPPHD13 = GTMtempHD13.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()j
# z_indexes = abs(GTMtempHD13-mld_tempKPPHD13).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPHD13 = GTMtempHD13.z[z_indexes]; 
# mldKPPHD13 = mldKPPHD13.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=1; nuh=4)
# mld_tempKPP14 = GTMtemp14.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()j
# z_indexes = abs(GTMtemp14-mld_tempKPP14).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP14 = GTMtemp14.z[z_indexes]; 
# mldKPP14 = mldKPP14.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP HD model (num=1; nuh=4)
# mld_tempKPPHD14 = GTMtempHD14.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()j
# z_indexes = abs(GTMtempHD14-mld_tempKPPHD14).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPHD14 = GTMtempHD14.z[z_indexes]; 
# mldKPPHD14 = mldKPPHD14.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=1; nuh=2) modified RiG #
# mld_tempKPPR12 = GTMtempR12.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempR12-mld_tempKPPR12).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPR12 = GTMtempR12.z[z_indexes]; 
# mldKPPR12 = mldKPPR12.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=1; nuh=2) modified oppFac RiG #
# mld_tempKPPRO12 = GTMtempRO12.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempRO12-mld_tempKPPRO12).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPRO12 = GTMtempRO12.z[z_indexes]; 
# mldKPPRO12 = mldKPPRO12.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=1; nuh=2) modified RiG with factor#
# mld_tempKPPRR12 = GTMtempRR12.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempRR12-mld_tempKPPRR12).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPRR12 = GTMtempRR12.z[z_indexes]; 
# mldKPPRR12 = mldKPPRR12.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# ###  num == 2
# ## using KPP model (num=2; nuh=0.25)
# mld_tempKPP2025 = GTMtemp2025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp2025-mld_tempKPP2025).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP2025 = GTMtemp2025.z[z_indexes]; 
# mldKPP2025 = mldKPP2025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=2; nuh=0.5)
# mld_tempKPP205 = GTMtemp205.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp205-mld_tempKPP205).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP205 = GTMtemp205.z[z_indexes]; 
# mldKPP205 = mldKPP205.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

## using KPP model (num=2; nuh=1)
mld_tempKPP21 = GTMtemp21.sel(z=-11.0, method='nearest') - 0.2; 
# z_indexes = abs(GTMtemp - mld_temp).z#.load()
z_indexes = abs(GTMtemp21-mld_tempKPP21).argmin('z'); # print(z_indexes)
print('GOTM temp MLD indexes', z_indexes)
mldKPP21 = GTMtemp21.z[z_indexes]; 
mldKPP21 = mldKPP21.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=2; nuh=2)
# mld_tempKPP22 = GTMtemp22.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp22-mld_tempKPP22).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP22 = GTMtemp22.z[z_indexes]; 
# mldKPP22 = mldKPP22.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=2; nuh=4)
# mld_tempKPP24 = GTMtemp24.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp24-mld_tempKPP24).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP24 = GTMtemp24.z[z_indexes]; 
# mldKPP24 = mldKPP24.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=2; nuh=1) modified RiG
# mld_tempKPPR21 = GTMtempR21.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempR21-mld_tempKPPR21).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPR21 = GTMtempR21.z[z_indexes]; 
# mldKPPR21 = mldKPPR21.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=2; nuh=1) modified RiG
# mld_tempKPPRO21 = GTMtempRO21.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempRO21-mld_tempKPPRO21).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPRO21 = GTMtempRO21.z[z_indexes]; 
# mldKPPRO21 = mldKPPRO21.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=2; nuh=1) modified RiG with factor
# mld_tempKPPRR21 = GTMtempRR21.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempRR21-mld_tempKPPRR21).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPRR21 = GTMtempRR21.z[z_indexes]; 
# mldKPPRR21 = mldKPPRR21.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# ###  num == 4
# ## using KPP model (num=4; nuh=0.25)
# mld_tempKPP4025 = GTMtemp4025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp4025-mld_tempKPP4025).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP4025 = GTMtemp4025.z[z_indexes]; 
# mldKPP4025 = mldKPP4025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=4; nuh=0.5)
# mld_tempKPP405 = GTMtemp405.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp405-mld_tempKPP405).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP405 = GTMtemp405.z[z_indexes]; 
# mldKPP405 = mldKPP405.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=4; nuh=1)
# mld_tempKPP41 = GTMtemp41.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp41-mld_tempKPP41).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP41 = GTMtemp41.z[z_indexes]; 
# mldKPP41 = mldKPP41.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=4; nuh=1) modified Rig
# mld_tempKPPR41 = GTMtempR41.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtempR41-mld_tempKPPR41).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPPR41 = GTMtempR41.z[z_indexes]; 
# mldKPPR41 = mldKPPR41.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=4; nuh=2)
# mld_tempKPP42 = GTMtemp42.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp42-mld_tempKPP42).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP42 = GTMtemp42.z[z_indexes]; 
# mldKPP42 = mldKPP42.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=4; nuh=4)
# mld_tempKPP44 = GTMtemp44.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp44-mld_tempKPP44).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP44 = GTMtemp44.z[z_indexes]; 
# mldKPP44 = mldKPP44.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


# ###  num == 5.8
# ## using KPP model (num=5.8; nuh=0.25)
# mld_tempKPP58025 = GTMtemp58025.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp58025-mld_tempKPP58025).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP58025 = GTMtemp58025.z[z_indexes]; 
# mldKPP58025 = mldKPP58025.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=5.8; nuh=0.5)
# mld_tempKPP5805 = GTMtemp5805.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp5805-mld_tempKPP5805).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP5805 = GTMtemp5805.z[z_indexes]; 
# mldKPP5805 = mldKPP5805.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=5.8; nuh=1)
# mld_tempKPP581 = GTMtemp581.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp581-mld_tempKPP581).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP581 = GTMtemp581.z[z_indexes]; 
# mldKPP581 = mldKPP581.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=5.8; nuh=2)
# mld_tempKPP582 = GTMtemp582.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp582-mld_tempKPP582).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP582 = GTMtemp582.z[z_indexes]; 
# mldKPP582 = mldKPP582.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)

# ## using KPP model (num=5.8; nuh=4)
# mld_tempKPP584 = GTMtemp584.sel(z=-11.0, method='nearest') - 0.2; 
# # z_indexes = abs(GTMtemp - mld_temp).z#.load()
# z_indexes = abs(GTMtemp584-mld_tempKPP584).argmin('z'); # print(z_indexes)
# print('GOTM temp MLD indexes', z_indexes)
# mldKPP584 = GTMtemp584.z[z_indexes]; 
# mldKPP584 = mldKPP584.isel(lat=0,lon=0); # print('mixed layer depth :',mld, mld.shape)


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
# to get a list, instead of a generator, use
# x_y = list(read__lines())
### Plot of MLD from OSMOSIS Observations
# plt.plot(j,np.negative(i),color='olivedrab',linestyle='dotted'); plt.xlabel('Time -Julian day'); plt.ylabel('Mixed Layer Depth -m')
# plt.grid()
# # plt.axis([260,341,-400,0]); 
# # plt.title('OSMOSIS Mixed Layer Depth')
# plt.savefig('OSMOSIS_Mixed_Layer_Depth_.png'); plt.show()

## ------------------------------- moving average MLDs--------------------------------- ##

def moving_average(i, n=3) :
    ret = np.cumsum(i, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# i = np.arange(10); print('set :',i); example = moving_average(i); print('moving average example :', example)

df_mld = moving_average(i)
df_mld5 = moving_average(i,n=5)
df_mld7 = moving_average(i,n=7)

df_dats = moving_average(j)
df_dats5 = moving_average(j,n=5)
df_dats7 = moving_average(j,n=7)


print('manualy moving average :', len(df_mld))
# df_mld = pd.read_csv('mld.csv')
# print(df_mld.info())
# print(df_mld)
# df_mld = pd.DataFrame(i)
# df_mld = i.rolling(3,min_periods=1).mean()
# print('moving average MLD', df_mld)

### Plots-compare all MLDs ###
plt.rcParams['figure.figsize'] = [15, 10]
# plt.plot(j,np.negative(i), color='gold',linestyle='solid',label = 'OSMOSIS');
# plt.plot(df_dats,np.negative(df_mld), color='green',linestyle='solid',label = 'OSMOSIS moving average n=3'); 
plt.plot(df_dats5,np.negative(df_mld5), color='red',linestyle='solid',label = 'OSMOSIS moving average n=5'); 
# plt.plot(df_dats7,np.negative(df_mld7), color='red',linestyle='solid',label = 'OSMOSIS moving average n=7'); 
# plt.plot(new_time,mldOS, color='salmon',linestyle='dotted',label = 'OSMOSIS-calculated'); plt.legend(loc=(0.10,0.15)); 
plt.plot(time2,mldKPP,color='olivedrab',linestyle='solid',label = r'GOTM KPP $K_m$; $K_h$');
# plt.plot(time2La,mldKPPLa,color='gray',linestyle='dotted',label = 'KPP EV');
# plt.plot(time2LaRib,mldKPPLaRib,color='gray',linestyle='dashed',label = 'KPP Lat Rib Cm=Ch=1');
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
# plt.plot(time_025025,mldKPP025025,color='steelblue',linestyle='solid',label = 'KPP num=0.25; nuh=0.25');
# plt.plot(time_02505,mldKPP02505,color='olivedrab',linestyle='solid',label = 'KPP num=0.25; nuh=0.5');
# plt.plot(time_0251,mldKPP0251,color='salmon',linestyle='solid',label = 'KPP num=0.25; nuh=1');
# plt.plot(time_0252,mldKPP0252,color='darkorange',linestyle='solid',label = 'KPP num=0.25; nuh=2');
# plt.plot(time_0254,mldKPP0254,color='rosybrown',linestyle='solid',label = 'KPP num=0.25; nuh=4');
# plt.plot(time_05025,mldKPP05025,color='steelblue',linestyle='solid',label = 'KPP num=0.5; nuh=0.25');
# plt.plot(time_0505,mldKPP0505,color='olivedrab',linestyle='solid',label = 'KPP num=0.5; nuh=0.5');
# plt.plot(time_051,mldKPP051,color='salmon',linestyle='solid',label = 'KPP num=0.5; nuh=1');
# plt.plot(time_052,mldKPP052,color='darkorange',linestyle='solid',label = 'KPP num=0.5; nuh=2');
# plt.plot(time_054,mldKPP054,color='rosybrown',linestyle='solid',label = 'KPP num=0.5; nuh=4');
# plt.plot(time_1025,mldKPP1025,color='steelblue',linestyle='solid',label = 'KPP num=1; nuh=0.25');
# plt.plot(time_105,mldKPP105,color='olivedrab',linestyle='solid',label = 'KPP num=1; nuh=0.5');
# plt.plot(time_11,mldKPP11,color='salmon',linestyle='solid',label = 'KPP num=1; nuh=1');
plt.plot(time_12,mldKPP12,color='steelblue',linestyle='solid',label = r'GOTM KPP $K_m$; $2*K_h$');
# plt.plot(time_R12,mldKPPR12,color='slateblue',linestyle='dashed',label = 'KPP num=1; nuh=2 mod RiG');
# plt.plot(time_RO12,mldKPPRO12,color='violet',linestyle='dotted',label = 'KPP num=1; nuh=2 mod Ric=0.6');
# plt.plot(time_RR12,mldKPPRR12,color='lawngreen',linestyle='dotted',label = 'KPP num=1; nuh=2 mod Ric=0.15');
# plt.plot(time_HD13,mldKPPHD13,color='dimgrey',linestyle='solid',label = 'KPP HD num=1; nuh=3');
# plt.plot(time_14,mldKPP14,color='rosybrown',linestyle='solid',label = 'KPP num=1; nuh=4');
# plt.plot(time_HD14,mldKPPHD14,color='chocolate',linestyle='solid',label = 'KPP HD num=1; nuh=4');
# plt.plot(time_2025,mldKPP2025,color='steelblue',linestyle='solid',label = 'KPP num=2; nuh=0.25');
# plt.plot(time_205,mldKPP205,color='olivedrab',linestyle='solid',label = 'KPP num=2; nuh=0.5');
plt.plot(time_21,mldKPP21,color='magenta',linestyle='solid',label = r'KPP $2*K_m$; $K_h$');
# plt.plot(time_R21,mldKPPR21,color='gray',linestyle='dashed',label = 'KPP num=2; nuh=1 mod RiG');
# plt.plot(time_RO21,mldKPPRO21,color='chocolate',linestyle='dotted',label = 'KPP num=2; nuh=1 mod Ric=0.15');
# plt.plot(time_RR21,mldKPPRR21,color='lawngreen',linestyle='dotted',label = 'KPP num=2; nuh=1 mod Ric=0.6');
# plt.plot(time_22,mldKPP22,color='darkorange',linestyle='solid',label = 'KPP num=2; nuh=2');
# plt.plot(time_24,mldKPP24,color='gray',linestyle='solid',label = 'KPP num=2; nuh=4');
# plt.plot(time_4025,mldKPP4025,color='steelblue',linestyle='solid',label = 'KPP num=4; nuh=0.25');
# plt.plot(time_405,mldKPP405,color='olivedrab',linestyle='solid',label = 'KPP num=4; nuh=0.5');
# plt.plot(time_41,mldKPP41,color='olivedrab',linestyle='solid',label = 'KPP num=4; nuh=1');
# plt.plot(time_R41,mldKPPR41,color='steelblue',linestyle='dashed',label = 'KPP num=4; nuh=1 mod RiG');
# plt.plot(time_42,mldKPP42,color='darkorange',linestyle='solid',label = 'KPP num=4; nuh=2');
# plt.plot(time_44,mldKPP44,color='rosybrown',linestyle='solid',label = 'KPP num=4; nuh=4');
# plt.plot(time_58025,mldKPP58025,color='steelblue',linestyle='solid',label = 'KPP num=5.8; nuh=0.25');
# plt.plot(time_5805,mldKPP5805,color='olivedrab',linestyle='solid',label = 'KPP num=5.8; nuh=0.5');
# plt.plot(time_581,mldKPP581,color='salmon',linestyle='solid',label = 'KPP num=5.8; nuh=1');
# plt.plot(time_582,mldKPP582,color='darkorange',linestyle='solid',label = 'KPP num=5.8; nuh=2');
# plt.plot(time_584,mldKPP584,color='rosybrown',linestyle='solid',label = 'KPP num=5.8; nuh=4');
# plt.plot(time2,mldKeps,color='goldenrod',linestyle='solid',label = 'GOTM k-eps');
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize='small');
plt.legend(loc='best',fontsize=10)
plt.grid()
plt.xlabel('Time -Julian days',fontsize=15); plt.ylabel('Mixed Layer Depth -m',fontsize=15)
# plt.xlim(530.9,539.1)
# plt.ylim(-200,0)
plt.xlim(480,540)
# plt.xlim(480,540)
# plt.axis([470,633,-150,0]);
# plt.axis([440,540,400,0]);
# plt.xlim(445.0,537.0)
# plt.axis([445.0,537.0,-150,0]);# plt.title('Compare OSMOSIS & GOTM-KPP,-k-epsilon MLDs')
# plt.tight_layout(rect=[0,0,0.75,1])
# plt.title(' z_ref=-11.0; DeltaT=0.2 ')
plt.savefig("GOTM_OSMOSIS__Mixed_Layer_Depth_Km_Kh.png"); 
plt.show(); #,bbox_inches='tight'

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
