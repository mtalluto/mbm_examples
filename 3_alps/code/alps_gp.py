#!/usr/bin/env python

#### if working in sublimetext
# import os
# os.chdir("..")


import sys
import warnings
# gpy prints lots of warnings during optimization; normally it is safe to ignore these
if '--warn' in sys.argv:
    warn = True
    warnings.simplefilter('default')
else:
    warn = False
    warnings.filterwarnings("ignore")


import re
# find where to look got mbmtools
try:
    p = re.compile('--dir=(.+)')
    mbmdir = filter(None, [p.match(x) for x in sys.argv])[0].group(1)
    sys.path.append(mbmdir)
except IndexError:
    mbmdir = None



import numpy as np
import mbmtools as mbm
import GPy

# xVars = ['dist_pToPET']
xVars = ['distance', 'precip_PET_ratio', 'Temp_Isothermallity', 'Temp_ColdestPeriod_min', 'Rad_Jan_total', 'geogDistance100km']
# xVars = ['dist_pToPET', 'pToPET', 'dist_isothermality', 'isothermality', 'dist_tempColdest', 'tempColdest', 'dist_radJan', 'radJan', 'geogDistance100km']

file = 'dat/tasFit_1.csv'
validDat = np.genfromtxt('dat/tasValid.csv', delimiter=',', skip_header=0, names=True, dtype=float)
dat = np.genfromtxt(file, delimiter=',', skip_header=0, names=True, dtype=float)

modTax = mbm.MBM(GPy.likelihoods.link_functions.Probit())
modTax.add_data(dat, xVars, 'taxSorensen')
modTax.add_data(validDat, datType='validation')
modTax.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
modTax.fit_model()
modTax.predict_all('res/taxo')

modFun = mbm.MBM(GPy.likelihoods.link_functions.Probit())
modFun.add_data(dat, xVars, 'functionalDistance')
modFun.add_data(validDat, datType='validation')
modFun.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
modFun.fit_model()
modFun.predict_all('res/func')
    
    

