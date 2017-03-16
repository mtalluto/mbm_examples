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


take_samples = ('--sample' in sys.argv)

fname = "dat/beta_fit.csv"
holdout_fname = "dat/beta_holdout.csv"
xVars = ['distance', 'bio3', 'bio6', 'bio12', 'bio15']
dat = np.genfromtxt(fname, delimiter=',', skip_header=0, names=True, dtype=float)
holdout_dat = np.genfromtxt(holdout_fname, delimiter=',', skip_header=0, names=True, dtype=float)
xData = mbm.prep_x(dat, xVars)
hldDat = mbm.prep_x(holdout_dat, xVars)

y = {'t': dat['sorensen_t'],
    'p': dat['mpd_bp'],
    'f': dat['mpd_bf']}
links = {'t': GPy.likelihoods.link_functions.Probit(),
        'p': GPy.likelihoods.link_functions.Identity(),
        'f': GPy.likelihoods.link_functions.Identity()}
inference = {'t': GPy.inference.latent_function_inference.Laplace(),
        'p': GPy.inference.latent_function_inference.ExactGaussianInference(),
        'f': GPy.inference.latent_function_inference.ExactGaussianInference()}
output_dir = {'t': 'res/betaTaxo', 'p': 'res/betaPhylo', 'f': 'res/betaFunc'}
likelihood = { k: GPy.likelihoods.Gaussian(gp_link=links[k]) for k in links}


for mName in y:
    print("Starting " + mName)
    xx = xData
    yy = y[mName]
    yy = np.reshape(yy, (np.shape(yy)[0], 1))
    kern = GPy.kern.RBF(input_dim=np.shape(xx)[1], ARD=True)

    # set some priors
    kern.variance.set_prior(GPy.priors.Gamma.from_EV(1.,10.))
    kern.lengthscale.set_prior(GPy.priors.Gamma.from_EV(1.,10.))

    # fit model and set up noise variance prior
    model = GPy.core.GP(xx, Y=yy, likelihood=likelihood[mName], inference_method=inference[mName], kernel=kern)
    model.Gaussian_noise.variance.set_prior(GPy.priors.Gamma.from_EV(1.,10.))
    print("  optimizing")
    model.optimize()

    # model.kern.plot()
    # model.plot()
    # slice = 0
    # figure = GPy.plotting.plotting_library().figure(3, 1)
    # canvas = model.plot(figure=figure, fixed_inputs=[(1,slice), (2, slice)], row=(1), plot_data=True)
    # canvas = model.plot(figure=figure, fixed_inputs=[(0,3.), (2, slice)], row=(2), plot_data=True)
    # canvas = model.plot(figure=figure, fixed_inputs=[(0,3.), (1, slice)], row=(3), plot_data=True)

    # save parameters
    mbm.make_dir(output_dir[mName])
    np.savetxt(output_dir[mName] + '/params.csv', model.param_array, delimiter=',')

    #### PREDICTION
    header = ','.join(xVars) + ',' + ','.join(['mean', 'sd', 'lower', 'upper'])
    
    ## predict to the data points and holdout data
    print("  predicting to data")
    samp = take_samples or not isinstance(model.likelihood, GPy.likelihoods.Gaussian)
    datPredict = mbm.mod_predict(model, xData, samp)
    hoPredict = mbm.mod_predict(model, hldDat, samp)
    np.savetxt(output_dir[mName] + '/datPredict.csv', datPredict, delimiter=',', header=header, comments='')
    np.savetxt(output_dir[mName] + '/holdoutPredict.csv', hoPredict, delimiter=',', header=header, comments='')
    
    
    ## predict 1-D response curves
    print("  response curve")
    sliceLen = 250
    predDat = [np.empty((sliceLen, xx.shape[1])) for i in range(xx.shape[1])]
    for i in range(xx.shape[1]):
        for j in range(xx.shape[1]):
            if i == j:
                predDat[i][:,j] = np.linspace(min(xx[:,i]), max(xx[:,i]), num=sliceLen)
            else:
                predDat[i][:,j] = np.full(sliceLen, np.round(np.mean(xx[:,j]), 3))
    predDat = np.concatenate(predDat, axis=0)
    respCurve = mbm.mod_predict(model, predDat, samp)
    np.savetxt(output_dir[mName] + '/respCurve.csv', respCurve, delimiter=',', header=header, comments='')
    
    ## predict 2-D grids -- note, this will need some work
    # gridDim = 100
    # grids = []
    # grMargins = [np.linspace(min(xx[:,i]), max(xx[:,i]), num=gridDim) for i in range(xx.shape[1])]
    # grids.append(mbm.expand_grid(grMargins[0], grMargins[1], np.round(np.mean(xx[:,2]), 3), ind=(0,1,2)))
    # grids.append(mbm.expand_grid(grMargins[0], grMargins[2], np.round(np.mean(xx[:,1]), 3), ind=(0,2,1)))
    # gridPredict = [mbm.mod_predict(model, gr) for gr in grids]
    # null = [np.savetxt(output_dir[mName] + '/grid' + str(i) + '.csv', gridPredict[i], delimiter=',', header=header, comments='') for i in range(len(gridPredict))]

