##############################################################################
# Some General use functions
# Authored by Ammar Mian, 17/06/2018
# e-mail: ammar.mian@centralesupelec.fr
##############################################################################
# Copyright 2018 @CentraleSupelec
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################
import numpy as np
import scipy as sp
from tqdm import tqdm
from multiprocessing import Process, Queue
import scipy.special


def multivariate_complex_normal_samples(mean, covariance, N, pseudo_covariance=0):
    """ A function to generate multivariate complex normal vectos as described in:
        Picinbono, B. (1996). Second-order complex random vectors and normal
        distributions. IEEE Transactions on Signal Processing, 44(10), 2637–2640.
        Inputs:
            * mean = vector of size p, mean of the distribution
            * covariance = the covariance matrix of size p*p(Gamma in the paper)
            * pseudo_covariance = the pseudo-covariance of size p*p (C in the paper)
                for a circular distribution omit the parameter
            * N = number of Samples
        Outputs:
            * Z = Samples from the complex Normal multivariate distribution, size p*N"""

    (p, p) = covariance.shape
    Gamma = covariance
    C = pseudo_covariance

    # Computing elements of matrix Gamma_2r
    Gamma_x = 0.5 * np.real(Gamma + C)
    Gamma_xy = 0.5 * np.imag(-Gamma + C)
    Gamma_yx = 0.5 * np.imag(Gamma + C)
    Gamma_y = 0.5 * np.real(Gamma - C)

    # Matrix Gamma_2r as a block matrix
    Gamma_2r = np.block([[Gamma_x, Gamma_xy], [Gamma_yx, Gamma_y]])

    # Generating the real part and imaginary part
    mu = np.hstack((mean.real, mean.imag))
    v = np.random.multivariate_normal(mu, Gamma_2r, N).T
    X = v[0:p, :]
    Y = v[p:, :]
    return X + 1j * Y


def multivariate_complex_t_samples(mean, covariance, N, df, pseudo_covariance=0):
    """ A function to generate multivariate complex t distributed vectors using the
    definition with a product of a multivaraite normal with an inverse chi2 distributed samples. 
    Inputs:
        * mean = vector of size p, mean of the distribution
        * covariance = the covariance matrix of size p*p
        * pseudo_covariance = the pseudo-covariance of size p*p
            for a circular distribution omit the parameter
        * df = degrees of freedom of the chi-squared distribution
        * N = number of Samples
    Outputs:
        * Z = Samples from the complex multivariate t distribution, size p*N"""


    if df == np.inf:
        x = 1
    else:
        x = np.random.chisquare(df, N)/df
    z = multivariate_complex_normal_samples(np.zeros(mean.shape), covariance, N, pseudo_covariance)
    return np.tile(mean.reshape((len(mean),1)),(1,N)) + z/np.sqrt(x)[None,:] 


def multivariate_complex_K_samples(mean, covariance, N, mu, b, pseudo_covariance=0):
    """ A function to generate multivariate complex K distributed vectors using the
    definition provided at page 27 of the Pd.d thesis:
    "Detection en environement non Gaussien", Emanuelle Jay. 
    Inputs:
        * mean = vector of size p, mean of the distribution
        * covariance = the covariance matrix of size p*p
        * pseudo_covariance = the pseudo-covariance of size p*p
            for a circular distribution omit the parameter
        * mu = Shape parameter
        * b = Scale parameter
        * N = number of Samples
    Outputs:
        * Z = Samples from the complex multivariate t distribution, size p*N"""

    x = np.random.gamma(mu, 2/(b**2), N)
    z = multivariate_complex_normal_samples(np.zeros(mean.shape), covariance, N, pseudo_covariance)
    return np.tile(mean.reshape((len(mean),1)),(1,N)) + z*np.sqrt(x)[None,:]   


def multivariate_complex_Cauchy_samples(mean, covariance, N, mu, b, pseudo_covariance=0):
    """ A function to generate multivariate complex Cauchy distributed vectors using the
    definition provided at page 26 of the Pd.d thesis:
    "Detection en environement non Gaussien", Emanuelle Jay. 
    Inputs:
        * mean = vector of size p, mean of the distribution
        * covariance = the covariance matrix of size p*p
        * pseudo_covariance = the pseudo-covariance of size p*p
            for a circular distribution omit the parameter
        * mu = Shape parameter
        * b = Scale parameter
        * N = number of Samples
    Outputs:
        * Z = Samples from the complex multivariate t distribution, size p*N"""

    x = np.random.gamma(mu, 2/(b**2), N)
    z = multivariate_complex_normal_samples(np.zeros(mean.shape), covariance, N, pseudo_covariance)
    return np.tile(mean.reshape((len(mean),1)),(1,N)) + z/np.sqrt(x)[None,:]    


def multivariate_complex_Laplace_samples(mean, covariance, N, beta, pseudo_covariance=0):
    """ A function to generate multivariate complex Cauchy distributed vectors using the
    definition provided at page 27 of the Pd.d thesis:
    "Detection en environement non Gaussien", Emanuelle Jay. 
    Inputs:
        * mean = vector of size p, mean of the distribution
        * covariance = the covariance matrix of size p*p
        * pseudo_covariance = the pseudo-covariance of size p*p
            for a circular distribution omit the parameter
        * beta = Scale parameter
        * N = number of Samples
    Outputs:
        * Z = Samples from the complex multivariate t distribution, size p*N"""

    x = np.random.exponential(beta, N)
    z = multivariate_complex_normal_samples(np.zeros(mean.shape), covariance, N, pseudo_covariance)
    return np.tile(mean.reshape((len(mean),1)),(1,N)) + z*np.sqrt(x)[None,:]    


def SCM(x, *args):
    """ A function that computes the SCM for covariance matrix estimation
            Inputs:
                * x = a matrix of size p*N with each observation along column dimension
            Outputs:
                * Sigma = the estimate"""

    (p, N) = x.shape
    return (x @ x.conj().T) / N


def TylerFixedPoint(x, tol=0.0001, iterMax=10, *args):
    """ A function that computes the Tyler Fixed Point Estimator for covariance matrix estimation
            Inputs:
                * x = a matrix of size p*N with each observation along column dimension
                * tol = tolerance for convergence of estimator
                * iterMax = number of maximum iterations
            Outputs:
                * Sigma = the estimate
                * error = the final error between two iterations
                * iteration = number of iterations til convergence """

    (p, N) = x.shape
    error = np.inf
    converged = False
    iteration = 0
    Sigma = np.eye(p)
    while not converged and iteration < iterMax:
        iteration = iteration + 1

        # Compute the quadratic form efficiently
        vx = np.dot(np.linalg.inv(np.linalg.cholesky(Sigma).conj().T), x)
        ax = np.mean(vx * np.conj(vx))
        xbis = x / np.sqrt(np.tile(ax, (p, 1)))
        Sigma_new = np.dot(xbis, xbis.conj().T) / N
        Sigma_new = p * Sigma_new / np.trace(Sigma_new)

        # Compute error
        error = np.linalg.norm(Sigma_new - Sigma, 'fro') / np.linalg.norm(Sigma, 'fro')
        converged = error < tol
        Sigma = Sigma_new

    if iteration == iterMax:
        print("Waring: Tyler Fixed point algorithm has not converged")

    return (Sigma, error, iteration)


def ToeplitzMatrix(rho, p):
    """ A function that computes a Hermitian semi-positive matrix.
            Inputs:
                * rho = a scalar
                * p = size of matrix
            Outputs:
                * the matrix """

    return sp.linalg.toeplitz(np.power(rho, np.arange(0, p)))

def vec(M):
    """ Vectorize the matrix M in input """
    return M.reshape(np.prod(M.shape))