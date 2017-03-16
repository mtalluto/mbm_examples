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

fname = "dat/alpha_fit.csv"
holdout_fname = "dat/alpha_holdout.csv"
xVars = {'rich': ['bio6', 'bio15'], 'simp': ['bio3', 'bio15'], 
        'f':['bio6', 'bio15'], 'p':['bio5', 'bio15']}

dat = np.genfromtxt(fname, delimiter=',', skip_header=0, names=True, dtype=float)
holdout_dat = np.genfromtxt(holdout_fname, delimiter=',', skip_header=0, names=True, dtype=float)
xData = {k: mbm.prep_x(dat, xVars[k]) for k in xVars}
hldDat = {k: mbm.prep_x(holdout_dat, xVars[k]) for k in xVars}

y = {'rich': dat['richness'],
    'simp': dat['simpson'],
    'p': dat['mpd_ap'],
    'f': dat['mpd_af']}
links = {'rich': GPy.likelihoods.link_functions.Log(),
        'simp': GPy.likelihoods.link_functions.Probit(), 
        'p': GPy.likelihoods.link_functions.Identity(),
        'f': GPy.likelihoods.link_functions.Identity()}
inference = {'rich': GPy.inference.latent_function_inference.Laplace(),
        'simp': GPy.inference.latent_function_inference.Laplace(),
        'p': GPy.inference.latent_function_inference.ExactGaussianInference(),
        'f': GPy.inference.latent_function_inference.ExactGaussianInference()}
output_dir = {
    'rich': 'res/alphaRichness',
    'simp': 'res/alphaSimpson',
    'p': 'res/alphaPhylo',
    'f': 'res/alphaFunc'}
likelihood = { k: GPy.likelihoods.Gaussian(gp_link=links[k]) for k in links}
print("Note: using log-link and Gaussian likelihood for species richness to avoid slow samples")
# likelihood['rich'] = GPy.likelihoods.Poisson(gp_link=links["rich"])

for mName in y:
    xx = xData[mName]
    yy = y[mName]
    yy = np.reshape(yy, (np.shape(yy)[0], 1))
    kern = GPy.kern.RBF(input_dim=np.shape(xx)[1], ARD=True)

    # set some priors
    kern.variance.set_prior(GPy.priors.Gamma.from_EV(1.,10.))
    kern.lengthscale.set_prior(GPy.priors.Gamma.from_EV(1.,10.))


    # fit model and set up noise variance prior
    model = GPy.core.GP(xx, Y=yy, likelihood=likelihood[mName], inference_method=inference[mName], kernel=kern)
    if isinstance(model.likelihood, GPy.likelihoods.Gaussian):
        model.Gaussian_noise.variance.set_prior(GPy.priors.Gamma.from_EV(1.,10.))
    model.optimize()

    # model.kern.plot()
    # model.plot()
    # slices = [-1, 0, 1.5]
    # figure = GPy.plotting.plotting_library().figure(3, 1)
    # for i, j in zip(range(3), slices):
    #     canvas = model.plot(figure=figure, fixed_inputs=[(0,j)], row=(i+1), plot_data=False)

 #    figure = GPy.plotting.plotting_library().figure(3, 1)
 #    for i, j in zip(range(3), slices):
 #        canvas = model.plot(figure=figure, fixed_inputs=[(1,j)], row=(i+1), plot_data=False)

    #### PREDICTION
    # save parameters
    mbm.make_dir(output_dir[mName])
    np.savetxt(output_dir[mName] + '/params.csv', model.param_array, delimiter=',')

    #### PREDICTION
    header = ','.join(xVars[mName]) + ',' + ','.join(['mean', 'sd', 'lower', 'upper'])
    
    ## predict to the data points and holdout data
    samp = take_samples or not isinstance(model.likelihood, GPy.likelihoods.Gaussian)
    datPredict = mbm.mod_predict(model, xData[mName], samp=samp)
    hoPredict = mbm.mod_predict(model, hldDat[mName], samp=samp)
    np.savetxt(output_dir[mName] + '/datPredict.csv', datPredict, delimiter=',', header=header, comments='')
    np.savetxt(output_dir[mName] + '/holdoutPredict.csv', hoPredict, delimiter=',', header=header, comments='')
    
    ## predict 1-D response curves
    sliceLen = 250
    predDat = [np.empty((sliceLen, xx.shape[1])) for i in range(xx.shape[1])]
    for i in range(xx.shape[1]):
        for j in range(xx.shape[1]):
            if i == j:
                predDat[i][:,j] = np.linspace(min(xx[:,i]), max(xx[:,i]), num=sliceLen)
            else:
                predDat[i][:,j] = np.full(sliceLen, np.round(np.mean(xx[:,j]), 3))
    predDat = np.concatenate(predDat, axis=0)
    respCurve = mbm.mod_predict(model, predDat, samp=samp)
    np.savetxt(output_dir[mName] + '/respCurve.csv', respCurve, delimiter=',', header=header, comments='')
    
    
    ## predict 2-D grids
    gridDim = 100
    grMargins = [np.linspace(min(xx[:,0]), max(xx[:,0]), num=gridDim), np.linspace(min(xx[:,1]), max(xx[:,1]), num=gridDim)]
    grid = mbm.expand_grid(grMargins[0], grMargins[1], ind=(0,1))
    gridPredict = mbm.mod_predict(model, grid, samp=samp)
    np.savetxt(output_dir[mName] + '/grid.csv', gridPredict, delimiter=',', header=header, comments='')
    

