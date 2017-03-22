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


import numpy as np
import mbmtools as mbm
import GPy


# take_samples = ('--sample' in sys.argv)

# decide which models to run
doRich = '--alpha' in sys.argv or '--richness' in sys.argv
doSimp = '--alpha' in sys.argv or '--simp' in sys.argv
doFA = '--alpha' in sys.argv or '--fa' in sys.argv
doPA = '--alpha' in sys.argv or '--pa' in sys.argv


fname = "dat/alpha_fit.csv"
validFname = "dat/alpha_valid.csv"
dat = np.genfromtxt(fname, delimiter=',', skip_header=0, names=True, dtype=float)
validDat = np.genfromtxt(validFname, delimiter=',', skip_header=0, names=True, dtype=float)

if doRich:
    print "Starting richness"
    rich = mbm.MBM(GPy.likelihoods.link_functions.Log(), GPy.likelihoods.Poisson)
    rich.add_data(dat, ['bio6', 'bio15'], 'richness')
    rich.add_data(validDat, datType='validation')
    # fit optimal model with vague priors
    rich.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
    rich.fit_model()
    print "  predicting"
    rich.predict_all("res/richness")

    # turn off optimization and try tweaking lengthscale
    rich.set_lengthscale(1, which = 0)
    rich.predict_all("res/richness_1")

    rich.set_lengthscale(0.4, which = 0)
    rich.predict_all("res/richness_0.4")



if doSimp:
    print "Starting simpson"
    simp = mbm.MBM(linkFunction = GPy.likelihoods.link_functions.Probit())
    simp.add_data(dat, ['bio3', 'bio13'], 'simpson')
    simp.add_data(validDat, datType='validation')
    simp.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
    simp.fit_model()
    print "  predicting"
    simp.predict_all("res/simpson")

if doFA:
    print "Starting FA"
    funA = mbm.MBM()
    funA.add_data(dat, ['bio6', 'bio15'], 'mpd_af')
    funA.add_data(validDat, datType='validation')
    funA.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
    funA.fit_model()
    print "  predicting"
    funA.predict_all("res/funcAlpha")

if doPA:
    print "Starting PA"
    phyA = mbm.MBM()
    phyA.add_data(dat, ['bio5', 'bio15'], 'mpd_ap')
    phyA.add_data(validDat, datType='validation')
    phyA.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,3.))
    phyA.fit_model()
    print "  predicting"

    # phyA.set_lengthscale(0.4, which = 0)
    # phyA.set_lengthscale(0.4, which = 1)
    phyA.predict_all("res/phyloAlpha")




## some old code for plotting here, might incorporate into class later
# import matplotlib;matplotlib.rcParams['figure.figsize'] = (8,5)
# from matplotlib import pyplot as plt
# import pylab
# pylab.ion()

# model.kern.plot()
# model.plot()
# slices = [-1, 0, 1.5]
# figure = GPy.plotting.plotting_library().figure(3, 1)
# for i, j in zip(range(3), slices):
#     canvas = model.plot(figure=figure, fixed_inputs=[(0,j)], row=(i+1), plot_data=True)

# figure = GPy.plotting.plotting_library().figure(3, 1)
# for i, j in zip(range(3), slices):
#     canvas = model.plot(figure=figure, fixed_inputs=[(1,j)], row=(i+1), plot_data=True)

    
    
    
    

