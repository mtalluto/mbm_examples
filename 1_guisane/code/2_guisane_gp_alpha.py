#!/usr/bin/env python

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


# import matplotlib;matplotlib.rcParams['figure.figsize'] = (8,5)
# from matplotlib import pyplot as plt
# import pylab
# pylab.ion()


#### if working in sublimetext
# import os
# os.chdir("..")

    
# take_samples = ('--sample' in sys.argv)

fname = "dat/alpha_fit.csv"
validFname = "dat/alpha_holdout.csv"
dat = np.genfromtxt(fname, delimiter=',', skip_header=0, names=True, dtype=float)
validDat = np.genfromtxt(validFname, delimiter=',', skip_header=0, names=True, dtype=float)

rich = mbm.MBM(linkFunction = GPy.likelihoods.link_functions.Log())
rich.add_data(dat, ['bio6', 'bio15'], 'richness')
rich.add_data(validDat, datType='validation')
rich.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,10.))
rich.fit_model()
rich.predict_all("res/richness")

# rich.sample = True
# rich.predict_all("res/richnessSamples")

richPois = mbm.MBM(GPy.likelihoods.link_functions.Log(), GPy.likelihoods.Poisson)
richPois.add_data(dat, ['bio6', 'bio15'], 'richness')
richPois.add_data(validDat, datType='validation')
richPois.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,10.))
richPois.fit_model()
richPois.predict_all("res/richnessPois")

richPois.sample = True
richPois.predict_all("res/richnessPoisSamples")


if "--simpson" in sys.argv:
    simp = mbm.MBM(linkFunction = GPy.likelihoods.link_functions.Probit())
    simp.add_data(dat, ['bio3', 'bio15'], 'simpson')
    simp.add_data(validDat, datType='validation')
    simp.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,10.))
    simp.fit_model()
    simp.predict_all("res/simpson")

if "--fa" in sys.argv:
    funA = mbm.MBM()
    funA.add_data(dat, ['bio6', 'bio15'], 'mpd_af')
    funA.add_data(validDat, datType='validation')
    funA.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,10.))
    funA.fit_model()
    funA.predict_all("res/funcAlpha")

if "--pa" in sys.argv:
    phyA = mbm.MBM()
    phyA.add_data(dat, ['bio5', 'bio15'], 'mpd_ap')
    phyA.add_data(validDat, datType='validation')
    phyA.add_prior(allPriors = GPy.priors.Gamma.from_EV(1.,10.))
    phyA.fit_model()
    phyA.predict_all("res/phyloAlpha")




## some old code for plotting here, might incorporate into class later
    # model.kern.plot()
    # model.plot()
    # slices = [-1, 0, 1.5]
    # figure = GPy.plotting.plotting_library().figure(3, 1)
    # for i, j in zip(range(3), slices):
    #     canvas = model.plot(figure=figure, fixed_inputs=[(0,j)], row=(i+1), plot_data=False)

 #    figure = GPy.plotting.plotting_library().figure(3, 1)
 #    for i, j in zip(range(3), slices):
 #        canvas = model.plot(figure=figure, fixed_inputs=[(1,j)], row=(i+1), plot_data=False)

    
    
    
    

