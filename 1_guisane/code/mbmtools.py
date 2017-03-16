import os
import numpy as np
import errno


def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def mod_predict(mod, newX, samp=False, nsamp=1000, pcts=(2.5, 97.5)):
    if samp:
        samples = mod.posterior_samples_f(newX, nsamp)
        mean = np.reshape(np.mean(samples, axis=1), (-1,1))
        variance = np.reshape(np.var(samples, axis=1), (-1,1))
        quants = np.transpose(np.percentile(samples, pcts, axis=1))
    else:
        mean, variance = mod.predict_noiseless(newX)
        quants = np.transpose(np.squeeze(mod.predict_quantiles(newX, pcts)))
    sd = np.sqrt(variance)
    
    return np.concatenate((newX, mean, sd, quants), axis=1)

def expand_grid(d1, d2, c = None, ind = (0,1,2)):
    """
    Returns a list of all permutations of the arrays contained in values with an optional constant c
    
    d1, d2: 1D numpy arrays
    c: optional constant
    ind: indices in the output array for d1, d2, and c
    """
    rDim = 2 if c is None else 3
    res = np.empty((np.shape(d1)[0] * np.shape(d2)[0], rDim))
    k = 0
    for i in range(d1.shape[0]):
        for j in range(d2.shape[0]):
            res[k, ind[0]] = d1[i]
            res[k, ind[1]] = d2[j]
            if c is not None:
                res[k, ind[2]] = c
            k += 1
    return res


def prep_x(arr, xvars):
    """
    sets up x array for GPy by pulling out xvars from arr and reshaping appropriately
    """
    result = arr[xvars]
    return result.view(np.float).reshape(result.shape + (-1,))

