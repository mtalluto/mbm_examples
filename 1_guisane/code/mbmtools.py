import os
import numpy as np
import errno
import GPy

class MBM(object):
    def __init__(self, linkFunction = None, likelihoodClass = None):
        self.link = GPy.likelihoods.link_functions.Identity() if linkFunction is None else linkFunction
        self.set_likelihood(likelihoodClass)
        self.priors = dict()
        self.model = None
        self.sample = False
        self.validX = None

    def add_data(self, dat, xvars=None, yvar=None, datType='fit'):
        if xvars is not None:
            self.xvars = xvars
        if yvar is not None:
            self.yvar = yvar

        if datType == 'fit':
            self.X = prep_x(dat, self.xvars)
            self.set_y(dat[self.yvar])
            self.kern = GPy.kern.RBF(input_dim=np.shape(self.X)[1], ARD=True)
        elif datType == 'validation':
            self.validX = prep_x(dat, self.xvars)
        else:
            raise ValueError("datType must be one of 'fit' or 'validation'")

    def set_likelihood(self, likelihoodClass, inference = None):
        if likelihoodClass is None:
            likelihoodClass = GPy.likelihoods.Gaussian
        self.likelihood = likelihoodClass(gp_link = self.link)

        if inference is None:
            if isinstance(self.link, GPy.likelihoods.link_functions.Identity) and \
                    isinstance(self.likelihood, GPy.likelihoods.Gaussian):
                self.inference = GPy.inference.latent_function_inference.ExactGaussianInference()
            else:
                self.inference = GPy.inference.latent_function_inference.Laplace()
        else:
            self.inference = inference

    def set_y(self, yy):
        if len(np.shape(yy)) != 2:
            yy = np.reshape(yy, (-1, 1))
        self.Y = yy

    def add_prior(self, variance = None, lengthscale = None, noiseVariance = None, allPriors = None):
        if allPriors is not None:
            variance = lengthscale = noiseVariance = allPriors
        if variance is not None:
            self.priors["variance"] = variance
        if lengthscale is not None:
            self.priors["lengthscale"] = lengthscale
        if noiseVariance is not None:
            self.priors["Gaussian_noise.variance"] = noiseVariance

    def set_priors(self):
        if 'variance' in self.priors.keys():
            self.kern.variance.set_prior(self.priors['variance'])
        if 'lengthscale' in self.priors.keys():
            self.kern.lengthscale.set_prior(self.priors['lengthscale'])
        if 'Gaussian_noise.variance' in self.priors.keys() and self.model is not None and \
                isinstance(self.model.likelihood, GPy.likelihoods.Gaussian):
            self.model.Gaussian_noise.variance.set_prior(self.priors['Gaussian_noise.variance'])

    def fit_model(self, optimize=True):
        self.model = GPy.core.GP(X=self.X, Y=self.Y, likelihood=self.likelihood, \
                inference_method=self.inference, kernel=self.kern)
        self.set_priors()
        if optimize:
            self.model.optimize()

    def predict_all(self, outputDir):
        make_dir(outputDir)
        self.save_params(outputDir)
        self.predict_to_data(outputDir)
        self.resp_curve_1d(outputDir)
        self.resp_curve_2d(outputDir)

    def save_params(self, outputDir):
        """
        Save model parameter array to the specified directory
        """
        np.savetxt(outputDir + '/params.csv', self.model.param_array, delimiter=',')

    def predict_to_data(self, outputDir):
        self.predict_and_save(self.X, outputDir + '/datPredict.csv')
        if self.validX is not None:
            self.predict_and_save(self.validX, outputDir + 'validPredict.csv')

    def predict_and_save(self, newX, outputFile, nsamp=1000, pcts=(2.5, 97.5)):
        """
        Predict the model to new data and save to a file
        GPy.predict_quantiles includes likelihood variance, so instead we use the
        analytical 95% confidence interval if we are predicting without samples

        newX: new data to predict to, must have same dims as original calibration data
        outputFile: destination for saving
        nsamp: if self.sample is True, the number of samples that will be taken
        pcts: if self.sample is True, the percentiles to return

        value: n by 4 numpy array; columns are mean, standard deviation, and quantiles
        """
        header = ','.join(self.xvars) + ',' + ','.join(['mean', 'sd', 'lower', 'upper'])
        if self.sample:
            samples = self.model.posterior_samples_f(newX, nsamp)
            mean = np.reshape(np.mean(samples, axis=1), (-1,1))
            sd = np.sqrt(np.reshape(np.var(samples, axis=1), (-1,1)))
            quants = np.transpose(np.percentile(samples, pcts, axis=1))
        else:
            mean, variance = self.model.predict_noiseless(newX)
            sd = np.sqrt(variance)
            lower = mean - 1.96 * sd
            upper = mean + 1.96 * sd
            quants = np.concatenate([lower, upper], axis=1)
        predictions = np.concatenate((newX, mean, sd, quants), axis=1)
        np.savetxt(outputFile, predictions, delimiter = ',', header=header, comments='')

    def resp_curve_1d(self, outputDir, sliceLen = 250):
        nv = self.X.shape[1]
        predDat = [np.empty((sliceLen, nv)) for i in range(nv)]
        for i in range(nv):
            for j in range(nv):
                if i == j:
                    predDat[i][:,j] = np.linspace(min(self.X[:,i]), max(self.X[:,i]), num=sliceLen)
                else:
                    predDat[i][:,j] = np.full(sliceLen, np.round(np.mean(self.X[:,j]), 3))
        predDat = np.concatenate(predDat, axis=0)
        self.predict_and_save(predDat, outputDir + '/respCurve1D.csv')

    def resp_curve_2d(self, outputDir, gridDim = 100):
        if(self.X.shape[1] > 2):
            print("2d response curves not yet supported with >2 x variables")
            return()
        gradient1 = np.linspace(min(self.X[:,0]), max(self.X[:,0]), num=gridDim)
        gradient2 = np.linspace(min(self.X[:,1]), max(self.X[:,1]), num=gridDim)
        grid = expand_grid(gradient1, gradient2, ind=(0,1))
        self.predict_and_save(grid, outputDir + "/responseCurve2D.csv")





def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise



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

