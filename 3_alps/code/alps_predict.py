#!/usr/bin/env python
import numpy as np
import re
import sys
import os

# find where to look got mbmtools
try:
    p = re.compile('--dir=(.+)')
    mbmdir = filter(None, [p.match(x) for x in sys.argv])[0].group(1)
    sys.path.append(mbmdir)
except IndexError:
    mbmdir = None

# get prediction file(s)
try:
    p = re.compile('--predict=(.+)')
    predictFiles = [y.group(1) for y in filter(None, [p.match(x) for x in sys.argv])]
except IndexError:
    print("Usage: " + sys.argv[0] + " --predict=predictFile1.csv --predict=predictFile2.csv ...")
    sys.exit(1)

import GPy
import mbmtools as mbm

xVars = ['distance', 'bio_4', 'bio_6', 'bio_7', 'bio_15']
# mods = ['sor', 'f_mpd', 'p_mpd']
mods = ['sor']
file = 'dat/betaDiv.csv'
dat = np.genfromtxt(file, delimiter=',', skip_header=0, names=True, dtype=float)
predictDatSets = [np.genfromtxt(f, delimiter=',', skip_header=0, names=True, dtype=float) for f in predictFiles]

for m in mods:
    print "starting model " + m
    # create result directories
    rdir = 'res/' + m
    if not os.path.isdir(rdir):
        os.makedirs(rdir)
    if m == 'sor':
        link = GPy.likelihoods.link_functions.Probit()
        mf = True
    else:
        link = GPy.likelihoods.link_functions.Log()
        mf = False
    mod = mbm.MBM(link, linearMeanFunction = mf, meanFunctionSlope = 0.25)
    mod.add_data(dat, xVars, m)
    # mod.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
    par_file = rdir + '/params.csv'
    mod.load_params(par_file)
    prCounter = 0
    for p_dat in predictDatSets:
        print "   predicting to " + str(np.shape(p_dat)[0]) + " other points"
        predict = mod.predict(newX = p_dat, header=True)
        np.savetxt(rdir + "/predict" + str(prCounter) + ".csv", predict[1], header=predict[0], delimiter=',', comments='')
        prCounter += 1
