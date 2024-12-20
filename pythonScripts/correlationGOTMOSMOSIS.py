import numpy as np
import os
import sys
from scipy import stats
import matplotlib.pyplot as plt
# sys.path.append('/home/users/mc837749/Documents/gotm-4.0.0/simulations/annualOSMOSISforcing_2502/')
# from heatFlux_osmosisAnnual import *

#################################################
### -------- correlation coefficient -------- ###
#################################################

## r = sum(GOTM*Glider)/sqrt(sum(GOTM^2)*sum(Glider^2))

##### Glider heat Budget #####

Gliderdays = []; GliderHeatBudget = [];

with open('heatFlux_osmosisAnnual/2dGlider_HB.txt', 'r') as f:
	for row in f:
		a, b = row.split()
		Gliderdays.append(a)
		GliderHeatBudget.append(b)
Gliderdays = np.array(Gliderdays, dtype=float); GliderHeatBudget = np.array(GliderHeatBudget, dtype=float)

# print('Glider heat budget :\n', GliderHeatBudget)

##### Glider heat Budget #####

datFiledays = []; datFileHeatBudget = [];

with open('heatFlux_osmosisAnnual/2dDAT_HB.txt', 'r') as k:
	for row in k:
		i, j = row.split()
		datFiledays.append(i)
		datFileHeatBudget.append(j)
datFiledays = np.array(datFiledays, dtype=float); datFileHeatBudget = np.array(datFileHeatBudget, dtype=float)

# print('Glider heat budget :\n', GliderHeatBudget)


##### GOTM heat Budget #####

GOTMdays = []; GOTMheatBudget = []

with open('heatFlux_osmosisAnnual/reducedGOTMHB.txt', 'r') as g:
	for row in g:
		c, d = row.split()
		GOTMdays.append(c)
		GOTMheatBudget.append(d)
GOTMdays = np.array(GOTMdays, dtype=float); GOTMheatBudget = np.array(GOTMheatBudget, dtype=float)
print('GOTM heat budget : \n', GOTMheatBudget)

GOTM2days = []; GOTM2heatBudget = []

with open('heatFlux_osmosisAnnual/reduced2GOTMHB.txt', 'r') as h:
	for row in h:
		e, f = row.split()
		GOTM2days.append(e)
		GOTM2heatBudget.append(f)
GOTM2days = np.array(GOTM2days, dtype=float); GOTM2heatBudget = np.array(GOTM2heatBudget, dtype=float)
print('GOTM heat budget : \n', GOTM2heatBudget)
print('print shapes (1) datfile (2) GOTM : ', datFileHeatBudget.shape, GOTM2heatBudget.shape)

corrCoeff1 = sum(np.multiply(GOTMheatBudget,GliderHeatBudget))/np.sqrt(np.sum(np.square(GOTMheatBudget))*np.sum(np.square(GliderHeatBudget)))
print('correlation coefficient between GOTM and OSMOSIS heat budget :', corrCoeff1)	

# Autumn 2012
corrCoeffAut12 = sum(np.multiply(datFileHeatBudget[8:274],GOTM2heatBudget[8:274]))/np.sqrt(np.sum(np.square(datFileHeatBudget[8:274]))*np.sum(np.square(GOTM2heatBudget[8:274])))
print('(Autumn 2012) correlation coefficient between GOTM and heatflux.dat heat budget :', corrCoeffAut12)

print('gotm heat budget',GOTM2heatBudget[8:274])
print('.dat file heat budget', datFileHeatBudget[8:274])


plt.plot(GOTM2days[8:274], GOTM2heatBudget[8:274], 'k-')#, label = 'd/dt int_mld^0(GOTM_theta)'); plt.legend(loc='lower right')
plt.plot(datFiledays[8:274], datFileHeatBudget[8:274], 'g-')#,label = '[H + R_n(0)]/rho*c_p'); plt.legend(loc='lower right')
plt.xlabel('Time -Julian days'); plt.ylabel('Surface heat budget (reduced) -W/m2');# plt.title('Surface Heat Flux - OSMOSIS Autumn observations')
# plt.axis([265,300,-0.02,0.03])
plt.savefig('heatBudgetMLD_GOTM_Glider_solarRadheatflux_comparison30days.png'); plt.show()

exit()

# Winter 2012
corrCoeffWin12 = sum(np.multiply(datFileHeatBudget[275:1424],GOTM2heatBudget[275:1424]))/np.sqrt(np.sum(np.square(datFileHeatBudget[275:1424]))*np.sum(np.square(GOTM2heatBudget[275:1424])))
print('(Winter 2012) correlation coefficient between GOTM and heatflux.dat heat budget :', corrCoeffWin12)	
# Spring 2013
corrCoeffSpr13 = sum(np.multiply(datFileHeatBudget[1425:2176],GOTM2heatBudget[1425:2176]))/np.sqrt(np.sum(np.square(datFileHeatBudget[1425:2176]))*np.sum(np.square(GOTM2heatBudget[1425:2176])))
print('(Spring 2013) correlation coefficient between GOTM and heatflux.dat heat budget :', corrCoeffSpr13)	
# Summer 2013
corrCoeffSum13 = sum(np.multiply(datFileHeatBudget[2177:2928],GOTM2heatBudget[2177:2928]))/np.sqrt(np.sum(np.square(datFileHeatBudget[2177:2928]))*np.sum(np.square(GOTM2heatBudget[2177:2928])))
print('(Summer 2013) correlation coefficient between GOTM and heatflux.dat heat budget :', corrCoeffSum13)	
exit()
corrCoeff = []

def correlationHB():
	for i in arange(1, 4095):
		corrCoeff = np.append(corrCoeff, sum(GOTMheatBudget[i]*GliderHeatBudget[i])/sqrt(sum(GOTMheatBudget[i]^2)*sum(GliderHeatBudget[i]^2)))



exit()
PearsoncorrCoeff = stats.pearsonr(GOTMheatBudget, GliderHeatBudget)
print('correlation coefficient between GOTM and OSMOSIS heat budget :\n', PearsoncorrCoeff)

exit()

