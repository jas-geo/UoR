
import numpy as np
import math
import sys
import os.path
import xarray as xr
from definitionsOSMOSIS import compareKPP, localKPPvisualize
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime, timedelta, date
from time import strftime
from netcdftime import utime#, #datetime
import matplotlib.pyplot as plt

# OSMOSIS time
Gldr_sel = [3910] 
# 850  -> 330.573599 (Autumn)
# 1588 -> 399.913923 (Winter)
# 2700 -> 500.523923 (Spring)
# 3910 -> 600.658726 (Summer)

GTM_sel = [48900]
# 10000 -> 330.451388 (Autumn)
# 20000 -> 399.895833 (Winter)
# 34500 -> 500.590277 (Spring)
# 48900 -> 600.590277 (Summer)

kppdata, kppdataL = compareKPP(Gldr_sel, GTM_sel)

# print('annual OSMOSIS observations',new_time.values)
# print('annual KPP simulations with non-local fluxes',kpptime2.values)
# print('annual KPP simulations with only local fluxes',localtime2.values)

localKPPvisualize(plt,kppdata,kppdataL, 'Glider day %s')