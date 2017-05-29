#!/usr/bin/env python
import os
import re
import numpy as np
import sys
import warnings

#### if working in sublimetext
# os.chdir("..")

# gpy prints lots of warnings during optimization; normally it is safe to ignore these
if '--warn' in sys.argv:
    warn = True
    warnings.simplefilter('default')
else:
    warn = False
    warnings.filterwarnings("ignore")


# find where to look got mbmtools
try:
    p = re.compile('--dir=(.+)')
    mbmdir = filter(None, [p.match(x) for x in sys.argv])[0].group(1)
    sys.path.append(mbmdir)
except IndexError:
    mbmdir = None



import GPy
import mbmtools as mbm

xVars = ['distance', 'bio_4', 'bio_6', 'bio_7', 'bio_15']
mods = ['LDMC_sor', 'LDMC_mpd', 'PL_VEG_H_sor', 'PL_VEG_H_mpd', 'SEEDM_sor', 'SEEDM_mpd', 'SLA_sor',
    'SLA_mpd', 'sor', 'f_sor', 'f_mpd', 'p_sor', 'p_mpd']
file = 'dat/betaDiv.csv'
validFile = 'dat/betaDiv_valid.csv'
# predictFile = 'dat/betaDiv_predict.csv'
dat = np.genfromtxt(file, delimiter=',', skip_header=0, names=True, dtype=float)
validDat = np.genfromtxt(validFile, delimiter=',', skip_header=0, names=True, dtype=float)
# predictDat = np.genfromtxt(predictFile, delimiter=',', skip_header=0, names=True, dtype=float)

for m in mods:
    print "starting model " + m
    # create result directories
    rdir = 'res/' + m
    if not os.path.isdir(rdir):
        os.makedirs(rdir)
    link = GPy.likelihoods.link_functions.Log() if re.match('mpd', m) else GPy.likelihoods.link_functions.Probit() 
    mod = mbm.MBM(link)
    mod.add_data(dat, xVars, m)
    mod.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
    mod.fit_model()
    mod.save_params(rdir)
    print "   predicting to " + str(np.shape(mod.X)[0]) + " datapoints"
    datPredict = mod.predict(header=True)
    np.savetxt(rdir + "/dat_predict.csv", datPredict[1], header=datPredict[0] )
    print "   predicting to " + str(np.shape(validDat)[0]) + " validation points"
    validPredict = mod.predict(newX=validDat, header=True)
    np.savetxt(rdir + "/valid_predict.csv", validPredict[1], header=validPredict[0] )
    # print "   predicting to " + str(np.shape(predictDat)[0]) + " other points"
    # predict = mod.predict(newX = predictDat, header=True)
    # np.savetxt(rdir + "/predict.csv", predict[1], header=predict[0] )

    

