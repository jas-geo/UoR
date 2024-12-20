import numpy as np
import math
import sys
import csv
from scipy.stats.stats import pearsonr
from scipy import interpolate
from sklearn.metrics import pairwise_distances_argmin
import os
import xarray as xr
import pandas as pd
import tarfile
import netCDF4
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
import matplotlib.dates as mdates
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, ScalarFormatter,LogFormatter)

# plt.switch_backend('agg')
np.set_printoptions(suppress=True, threshold=sys.maxsize, formatter={'float_kind':'{:f}'.format}) #  

#############################################################################################
### -------------------------------------- OSMOSIS -------------------------------------- ###
#############################################################################################
glider = xr.open_dataset('../gotm-4.0.0-kpp/simulations/glider_timeseries.nc', decode_times=False);

GldrDepth = np.negative(glider.pres_grid); GldrTemp = glider.potemp; GldrDens = glider.pot_den; 
Gldrdays = glider.dats; GldrSalt = glider.prac_sal; Gldro2 = glider.oxygen

## Convert data time into year days and assign it as a dimension ##
init_time = datetime(2012,1,1,0,0,0);# print('Jan 1 2012 : ',init_time)
initime_jul = init_time.toordinal();# print('Jan 1 2012 Julian : ',initime_jul)
initime_RD = initime_jul+365 ;# print('Jan 1 2012 RATA DIE : ',initime_RD)
# 
new_time = Gldrdays - np.full_like(Gldrdays, initime_RD) # print(new_time[194:938]) year days
glider = glider.assign(new_time=("time", new_time))
glider.new_time.attrs['Unit'] = 'days since 01/01/2012 00:00:00 UTC'
glider = glider.set_coords('new_time', inplace=True); glider = glider.swap_dims({"time" : "new_time"})
# print('time', new_time[0:100].values)

print(glider['potemp'][:,:].shape)

import seaborn as sns
 
# plt.rcParams["figure.figsize"]=(7,5)

# # Function to show the heat map --------------------
# ax = sns.heatmap( GldrTemp.T , cmap = 'inferno' )
  
# # Adding details to the plot
# plt.title( "2-D OSMOSIS Heat Map" )
# plt.xlabel('time elements')
# plt.ylabel('depth elements')
# plt.show()

# # Function to show the saline map ----------------------
# ax = sns.heatmap( GldrSalt.T , cmap = 'inferno' )
  
# # Adding details to the plot
# plt.title( "2-D OSMOSIS Saline Map" )
# plt.xlabel('time elements')
# plt.ylabel('depth elements')
# plt.show()

# # Function to show the oxygen map ----------------------
# ax = sns.heatmap( Gldro2.T , cmap = 'inferno' )
  
# # Adding details to the plot
# plt.title( "2-D OSMOSIS Oxygen Map" )
# plt.xlabel('time elements')
# plt.ylabel('depth elements')
# plt.show()

# # Function to show the density map ----------------------
# ax = sns.heatmap( GldrDens.T , cmap = 'inferno' )
  
# # Adding details to the plot
# plt.title( "2-D OSMOSIS Density Map" )
# plt.xlabel('time elements')
# plt.ylabel('depth elements')
# plt.show()




#############################################################################################
### ----------------------------- Ocean Reanalysis System 5 ----------------------------- ###
#############################################################################################

ORA5_temp = xr.open_dataset('votemper_combined.nc', decode_times=False);
ORA5_salt = xr.open_dataset('vosaline_combined.nc', decode_times=False);
ORA5_wtrflx = xr.open_dataset('sowaflup_combined.nc', decode_times=False)

# # combine netcdf4 files with same predixes
# ORA5_wtrflx = xr.open_mfdataset('sowaflup_control_monthly_highres_2D_*.nc')
# ORA5_wtrflx = ORA5_wtrflx.to_netcdf('sowaflup_combined.nc')

ORA_t_Depth = ORA5_temp.deptht; ORA_t_Temp2 = ORA5_temp.votemper; ORA_t_time2 = ORA5_temp.time_counter
ORA_t_lat = ORA5_temp.nav_lat; ORA_t_lon = ORA5_temp.nav_lon;

# print('ORA5 time', ORAtime2.values)
# print('ORA5 dataset', ORA5_temp.variables.items())
# print('ORA5 dataset', ORA5_wtrflx.variables.items())
# print('ORAtemp3D items', ORA5_temp.variables.items())
# print('ORAtemp3D latitude', ORAlat[:,0:10], ORAlat.shape)
# print('ORAtemp3D longitude', ORAlon[:,0:10].values, ORAlon.shape)
# print('latitude', np.min(ORA_t_lat), np.max(ORA_t_lat))
# print('longitude', np.min(ORA_t_lon), np.max(ORA_t_lon)
#------------------------------------------------#
##### extract data for selected co-ordinates #####
#------------------------------------------------#

#TEMPERATURE
def extract(ORA5_temp, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_temp['nav_lon'].values.ravel(), ORA5_temp['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_temp['nav_lon'].shape))
    return ORA5_temp.isel(x=j0, y=i0).squeeze()
selectORATempO = extract(ORA5_temp, -16.25, 48.75);# print('extracts',selectORATemp)
selectORATemp1 = extract(ORA5_temp, -16.25, 48.50);#North
selectORATemp2 = extract(ORA5_temp, -16.00, 48.75);#East
selectORATemp3 = extract(ORA5_temp, -16.00, 48.50);#West
selectORATemp4 = extract(ORA5_temp, -16.25, 46.75);#South
print('select ORAtemp', selectORATempO.deptht.values)

diff = np.diff(selectORATempO.deptht.values)
print('difference',diff)
# selectORADepthO = selectORATempO.deptht; selectORAtimeO = selectORATempO.time_counter; 
# timeDepthO =np.meshgrid(selectORAtimeO,selectORADepthO)

# ### heat map ###--------------------------------------------------------------------------------------

# plt.rcParams["figure.figsize"]=(7,5)

# # Function to show the heat map
# ax = sns.heatmap( selectORATempO.votemper.T , cmap = 'inferno' )
  
# # Adding details to the plot
# plt.title( "2-D ORA5 Heat Map" )
# plt.xlabel('time elements')
# plt.ylabel('depth elements')
# plt.show()

#SALINITY
def extract(ORA5_salt, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_salt['nav_lon'].values.ravel(), ORA5_salt['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_salt['nav_lon'].shape))
    return ORA5_salt.isel(x=j0, y=i0).squeeze()
selectORASaltO = extract(ORA5_salt, -16.25, 48.75);# print('extracts',selectORASalt)
selectORASalt1 = extract(ORA5_salt, -16.25, 48.50);#North
selectORASalt2 = extract(ORA5_salt, -16.00, 48.75);#East
selectORASalt3 = extract(ORA5_salt, -16.00, 48.50);#West
selectORASalt4 = extract(ORA5_salt, -16.25, 46.75);#South

### saline map ###-----------------------------------------------------------------------------

plt.rcParams["figure.figsize"]=(7,5)

# Function to show the saline map
ax = sns.heatmap( selectORASaltO.vosaline.T , cmap = 'inferno' )
  
# Adding details to the plot
plt.title( "2-D ORA5 Saline Map" )
plt.xlabel('time elements')
plt.ylabel('depth elements')
plt.show()

#WATER FLUX
def extract(ORA5_wtrflx, nav_lon, nav_lat):
    index = pairwise_distances_argmin(X=[[nav_lon, nav_lat]], Y=np.c_[ORA5_wtrflx['nav_lon'].values.ravel(), ORA5_wtrflx['nav_lat'].values.ravel()])
    i0, j0 = np.unravel_index(index, (ORA5_wtrflx['nav_lon'].shape))
    return ORA5_wtrflx.isel(x=j0, y=i0).squeeze()
selectORAwFluxO = extract(ORA5_wtrflx, -16.25, 48.75);# print('extracts',selectORASalt)
selectORAwFlux1 = extract(ORA5_wtrflx, -16.25, 48.5);#North
selectORAwFlux2 = extract(ORA5_wtrflx, -16.0, 48.75);#East
selectORAwFlux3 = extract(ORA5_wtrflx, -16.0, 48.5);#West
# selectORAwFlux4 = extract(ORA5_wtrflx, -16.25, 46.75);#South

# print(selectORATempO.info())


## extract depth information for specific time counters
selectORATempO = selectORATempO.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORATemp1 = selectORATemp1.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORATemp2 = selectORATemp2.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORATemp3 = selectORATemp3.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORATemp4 = selectORATemp4.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)


selectORASaltO = selectORASaltO.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORASalt1 = selectORASalt1.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORASalt2 = selectORASalt2.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORASalt3 = selectORASalt3.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)
selectORASalt4 = selectORASalt4.sel(time_counter=212)# (S=0,O=30,N=61,D=91,J=122,F=152,M=181,A=212,M=242,J=273,J=303,A=334,S=365)

# --------------------------#
# #### plotting figures #####
# --------------------------#
#TEMPERATURE
plt.plot(selectORATempO.votemper,np.negative(selectORATempO.deptht), label = r'ORA5 (48.7 $\degree$N, 16.2 $\degree$W) 2013-04')
plt.plot(selectORATemp1.votemper,np.negative(selectORATemp1.deptht), linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W) 2013-04')
plt.plot(selectORATemp2.votemper,np.negative(selectORATemp2.deptht), linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W) 2013-04')
plt.plot(selectORATemp3.votemper,np.negative(selectORATemp3.deptht), linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W) 2013-04')
plt.plot(selectORATemp4.votemper,np.negative(selectORATemp4.deptht), linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W) 2013-04')
plt.plot(Apr13_temp,np.negative(glider.pres_grid), label = 'O5MOSIS 2013-04')
plt.legend(loc='lower right')
plt.ylabel('Depth -m'); plt.xlabel(r'Temperature - $\degree$ C')
plt.ylim(-1000,0); plt.xlim(7,19)
plt.grid()
plt.show()

#SALINITY
plt.plot(selectORASaltO.vosaline,np.negative(selectORASaltO.deptht), label = r'ORA5 (48.7 $\degree$N, 16.2 $\degree$W) 2013-04')
plt.plot(selectORASalt1.vosaline,np.negative(selectORASalt1.deptht), linestyle = 'dashed', label = r'ORA5 (50.7 $\degree$N, 16.2 $\degree$W) 2013-04')
plt.plot(selectORASalt2.vosaline,np.negative(selectORASalt2.deptht), linestyle = 'dotted', label = r'ORA5 (48.7 $\degree$N, 14.2 $\degree$W) 2013-04')
plt.plot(selectORASalt3.vosaline,np.negative(selectORASalt3.deptht), linestyle = 'dashed', label = r'ORA5 (48.7 $\degree$N, 18.2 $\degree$W) 2013-04')
plt.plot(selectORASalt4.vosaline,np.negative(selectORASalt4.deptht), linestyle = 'dotted', label = r'ORA5 (46.7 $\degree$N, 16.2 $\degree$W) 2013-04')
plt.plot(Apr13_sal,np.negative(glider.pres_grid), label = 'O5MOSIS 2013-04')
plt.legend(loc='lower right')
plt.ylabel('Depth -m'); plt.xlabel(r'Salinity - PSU')
plt.ylim(-1000,0); plt.xlim(35.1,35.8)
plt.grid()
plt.show()
month = np.array([datetime(2012,9,1),datetime(2012,10,1),datetime(2012,11,1),datetime(2012,12,1),
    datetime(2013,1,1),datetime(2013,2,1),datetime(2013,3,1),datetime(2013,4,1),datetime(2013,5,1),
    datetime(2013,6,1),datetime(2013,7,1),datetime(2013,8,1)])