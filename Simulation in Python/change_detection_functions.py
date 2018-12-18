##############################################################################
# Functions used to detect a change in the parameters of a SIRV distribution
# Authored by Ammar Mian, 28/09/2018
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
import warnings
from generic_functions import *

##############################################################################
# Gaussian Statistics
##############################################################################
def covariance_equality_glrt_gaussian_statistic(𝐗, args=None):
    """ GLRT statistic for detecting a change of covariance matrix in a multivariate Gaussian Time Series.
        At each time, Ni.i.d samples are available. A description of the statistic can be found in:
        D. Ciuonzo, V. Carotenuto and A. De Maio, 
        "On Multiple Covariance Equality Testing with Application to SAR Change Detection," 
        in IEEE Transactions on Signal Processing, vol. 65, no. 19, pp. 5078-5091, 1 Oct.1, 2017.
        doi: 10.1109/TSP.2017.2712124
        Inputs:
            * 𝐗 = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * args = None
        Outputs:
            * the GLRT statistic given the observations in input"""

    (p, N, T) = 𝐗.shape
    S = 0
    logDenominator = 0
    for t in range(0, T):
        St = SCM(𝐗[:, :, t])
        logDenominator = logDenominator + N * np.log(np.abs(np.linalg.det(St)))
        S = S + St / T
    logNumerator = N * T * np.log(np.abs(np.linalg.det(S)))
    if args is not None:
        if args=='log':
            return np.real(logNumerator - logDenominator)
    return np.exp(np.real(logNumerator - logDenominator))


def covariance_equality_t1_gaussian_statistic(𝐗, args=None):
    """ t1 statistic for detecting a change of covariance matrix in a multivariate Gaussian Time Series.
        At each time, Ni.i.d samples are available. A description of the statistic can be found in:
        D. Ciuonzo, V. Carotenuto and A. De Maio, 
        "On Multiple Covariance Equality Testing with Application to SAR Change Detection," 
        in IEEE Transactions on Signal Processing, vol. 65, no. 19, pp. 5078-5091, 1 Oct.1, 2017.
        doi: 10.1109/TSP.2017.2712124
        Inputs:
            * 𝐗 = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
        Outputs:
            * the t1 statistic given the observations in input"""

    (N, K, M) = 𝐗.shape

    Sigma_10 = SCM(𝐗.reshape((N, K*M)))
    iSigma_10 = np.linalg.inv(Sigma_10)
    t1 = 0
    for t in range(0, M):
        Sigma_m1 = SCM(𝐗[:, :, t])
        S = (iSigma_10 @ Sigma_m1)
        t1 = t1 + np.trace( S @ S )/M;

    if args is not None:
        if args=='log':
            return np.log(np.real(t1))
    return np.real(t1)


def covariance_equality_Wald_gaussian_statistic(𝐗, args=None):
    """ Wald statistic for detecting a change of covariance matrix in a multivariate Gaussian Time Series.
        At each time, Ni.i.d samples are available. A description of the statistic can be found in:
        D. Ciuonzo, V. Carotenuto and A. De Maio, 
        "On Multiple Covariance Equality Testing with Application to SAR Change Detection," 
        in IEEE Transactions on Signal Processing, vol. 65, no. 19, pp. 5078-5091, 1 Oct.1, 2017.
        doi: 10.1109/TSP.2017.2712124
        Inputs:
            * 𝐗 = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
        Outputs:
            * the Wald statistic given the observations in input"""

    (N, K, M) = 𝐗.shape
    L = 0;
    O = 0;
    Q = 0;
    Sigma_11 = SCM(𝐗[:, :, 0])
    for m in range(0,M):
        Sigma_m1 = SCM(𝐗[:, :, m])
        iSigma_m1 = np.linalg.inv(Sigma_m1)
        if m != 0:
            S = np.eye(N) - Sigma_11@iSigma_m1
            L = L + K*np.trace(S@S)
            Q = Q + K*(iSigma_m1 - iSigma_m1@Sigma_11@iSigma_m1)
        O = O + K*np.kron(iSigma_m1.T, iSigma_m1)
    
    if args is not None:
        if args=='log':
            return np.real(np.real(L - vec(Q).conj().T @ (np.linalg.inv(O)@vec(Q)))[0,0])
    return np.real(L - vec(Q).conj().T @ (np.linalg.inv(O)@vec(Q)))[0,0]


##############################################################################
# Robust Statistics
##############################################################################
def student_t_shape_statistic_d_known(X, Args):
    """ GLRT test for testing a change in the shape of a multivariate
        Student-t distribution when the degree of freefom is known.
        Inputs:
            * X = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * Args = d=degree of freedom and tol, iterMax for Tyler, scale
        Outputs:
            * the statistic given the observations in input"""

    d, tol, iter_max, scale = Args
    (p, N, T) = 𝐗.shape

    # Estimating 𝚺_0 using all the observations
    (𝚺_0, δ, niter) = student_t_estimator_covariance_mle(𝐗.reshape((p,T*N)), d, tol, iter_max)
    i𝚺_0 = np.linalg.inv(𝚺_0)

    # Some initialisation
    log_numerator_determinant_terms = T*N*np.log(np.abs(np.linalg.det(𝚺_0)))
    log_denominator_determinant_terms = 0
    log𝛕_0 = 0
    log𝛕_t = 0
    # Iterating on each date to compute the needed terms
    for t in range(0,T):
        # Estimating 𝚺_t
        (𝚺_t, δ, iteration) = student_t_estimator_covariance_mle(𝐗[:,:,t], d, tol, iter_max)

        # Computing determinant add adding it to log_denominator_determinant_terms
        log_denominator_determinant_terms = log_denominator_determinant_terms + \
                                            N*np.log(np.abs(np.linalg.det(𝚺_t)))

        # Computing quadratic terms
        log𝛕_0 =  log𝛕_0 + np.log(d + np.diagonal(𝐗[:,:,t].conj().T@i𝚺_0@𝐗[:,:,t]))
        log𝛕_t = log𝛕_t + np.log(d + np.diagonal(𝐗[:,:,t].conj().T@np.linalg.inv(𝚺_t)@𝐗[:,:,t]))

    # Computing quadratic terms
    log_numerator_quadtratic_terms = (d+p)*np.sum(log𝛕_0)
    log_denominator_quadtratic_terms = (d+p)*np.sum(log𝛕_t)

    # Final expression of the statistic
    if scale=='linear':
        λ = np.exp(np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms))
    else:
        λ = np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms)

    return λ


def scale_and_shape_equality_robust_statistic(𝐗, args):
    """ GLRT test for testing a change in the scale or/and shape of 
        a deterministic SIRV model.
        Inputs:
            * 𝐗 = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * args = tol, iter_max for Tyler, scale
        Outputs:
            * the statistic given the observations in input"""

    tol, iter_max, scale = args
    (p, N, T) = 𝐗.shape

    # Estimating 𝚺_0 using all the observations
    (𝚺_0, δ, niter) = tyler_estimator_covariance_matandtext(𝐗, tol, iter_max)
    i𝚺_0 = np.linalg.inv(𝚺_0)

    # Some initialisation
    log_numerator_determinant_terms = T*N*np.log(np.abs(np.linalg.det(𝚺_0)))
    log_denominator_determinant_terms = 0
    𝛕_0 = 0
    log𝛕_t = 0
    # Iterating on each date to compute the needed terms
    for t in range(0,T):
        # Estimating 𝚺_t
        (𝚺_t, δ, iteration) = tyler_estimator_covariance(𝐗[:,:,t], tol, iter_max)

        # Computing determinant add adding it to log_denominator_determinant_terms
        log_denominator_determinant_terms = log_denominator_determinant_terms + \
                                            N*np.log(np.abs(np.linalg.det(𝚺_t)))

        # Computing texture estimation
        𝛕_0 =  𝛕_0 + np.diagonal(𝐗[:,:,t].conj().T@i𝚺_0@𝐗[:,:,t]) / T
        log𝛕_t = log𝛕_t + np.log(np.diagonal(𝐗[:,:,t].conj().T@np.linalg.inv(𝚺_t)@𝐗[:,:,t]))

    # Computing quadratic terms
    log_numerator_quadtratic_terms = T*p*np.sum(np.log(𝛕_0))
    log_denominator_quadtratic_terms = p*np.sum(log𝛕_t)

    # Final expression of the statistic
    if scale=='linear':
        λ = np.exp(np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms))
    else:
        λ = np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms)

    return λ

def shape_equality_robust_statistic(𝐗, args):
    """ GLRT test for testing a change in the shape of 
        a deterministic SIRV model.
        Inputs:
            * 𝐗 = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * args = tol, iter_max for Tyler, scale
        Outputs:
            * the statistic given the observations in input"""

    tol, iter_max, scale = args
    (p, N, T) = 𝐗.shape

    # Estimating 𝚺_0 using all the observations
    (𝚺_0, δ, niter) = tyler_estimator_covariance(𝐗.reshape((p,T*N)), tol, iter_max)
    i𝚺_0 = np.linalg.inv(𝚺_0)

    # Some initialisation
    log_numerator_determinant_terms = T*N*np.log(np.abs(np.linalg.det(𝚺_0)))
    log_denominator_determinant_terms = 0
    log𝛕_0 = 0
    log𝛕_t = 0
    # Iterating on each date to compute the needed terms
    for t in range(0,T):
        # Estimating 𝚺_t
        (𝚺_t, δ, iteration) = tyler_estimator_covariance(𝐗[:,:,t], tol, iter_max)

        # Computing determinant add adding it to log_denominator_determinant_terms
        log_denominator_determinant_terms = log_denominator_determinant_terms + \
                                            N*np.log(np.abs(np.linalg.det(𝚺_t)))

        # Computing texture estimation
        log𝛕_0 =  log𝛕_0 + np.log(np.diagonal(𝐗[:,:,t].conj().T@i𝚺_0@𝐗[:,:,t]))
        log𝛕_t = log𝛕_t + np.log(np.diagonal(𝐗[:,:,t].conj().T@np.linalg.inv(𝚺_t)@𝐗[:,:,t]))

    # Computing quadratic terms
    log_numerator_quadtratic_terms = p*np.sum(log𝛕_0)
    log_denominator_quadtratic_terms = p*np.sum(log𝛕_t)

    # Final expression of the statistic
    if scale=='linear':
        λ = np.exp(np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms))
    else:
        λ = np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms)

    return λ


def scale_equality_robust_statistic(𝐗, args):
    """ GLRT test for testing a change in the scale of 
        a deterministic SIRV model.
        Inputs:
            * 𝐗 = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * args = tol, iter_max for Tyler, scale
        Outputs:
            * the statistic given the observations in input"""

    tol, iter_max, scale = args
    (p, N, T) = 𝐗.shape

    # Estimating 𝚺_t under H0 regime using all the observations
    (𝚺_0, 𝛅, niter) = tyler_estimator_covariance_text(𝐗, tol, iter_max)

    # Some initialisation
    log_numerator_determinant_terms = 0
    log_denominator_determinant_terms = 0
    𝛕_0 = 0
    log𝛕_t = 0
    # Iterating on each date to compute the needed terms
    for t in range(0,T):

        # Estimating 𝚺_t under H1 regime
        (𝚺_t, δ, iteration) = tyler_estimator_covariance(𝐗[:,:,t], tol, iter_max)

        # Computing determinant add adding it to log_numerator_determinant_terms
        log_numerator_determinant_terms = log_numerator_determinant_terms + \
                                        N*np.log(np.abs(np.linalg.det(𝚺_0[:,:,t])))

        # Computing determinant add adding it to log_denominator_determinant_terms
        log_denominator_determinant_terms = log_denominator_determinant_terms + \
                                            N*np.log(np.abs(np.linalg.det(𝚺_t)))

        # Computing texture estimation
        𝛕_0 =  𝛕_0 + np.diagonal(𝐗[:,:,t].conj().T@np.linalg.inv(𝚺_0[:,:,t])@𝐗[:,:,t]) / T
        log𝛕_t = log𝛕_t + np.log(np.diagonal(𝐗[:,:,t].conj().T@np.linalg.inv(𝚺_t)@𝐗[:,:,t]))

    # Computing quadratic terms
    log_numerator_quadtratic_terms = T*p*np.sum(np.log(𝛕_0))
    log_denominator_quadtratic_terms = p*np.sum(log𝛕_t)

    # Final expression of the statistic
    if scale=='linear':
        λ = np.exp(np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms))
    else:
        λ = np.real(log_numerator_determinant_terms - log_denominator_determinant_terms + \
        log_numerator_quadtratic_terms - log_denominator_quadtratic_terms)

    return λ

##############################################################################
# Some Functions
##############################################################################
def tyler_estimator_covariance_matandtext(𝐗, tol=0.0001, iter_max=20):
    """ A function that computes the Modified Tyler Fixed Point Estimator for 
    covariance matrix estimation under problem MatAndText.
        Inputs:
            * 𝐗 = a matrix of size p*N*T with each saptial observation along column dimension and time
                observation along third dimension.
            * tol = tolerance for convergence of estimator
            * iter_max = number of maximum iterations
        Outputs:
            * 𝚺 = the estimate
            * δ = the final distance between two iterations
            * iteration = number of iterations til convergence """

    (p, N, T) = 𝐗.shape
    δ = np.inf # Distance between two iterations
    𝚺 = np.eye(p) # Initialise estimate to identity
    iteration = 0

    # Recursive algorithm
    while (δ>tol) and iteration < iter_max:

        # Compute the textures for each pixel using all the dates avalaibe
        τ = 0
        i𝚺 = np.linalg.inv(𝚺)
        for t in range(0, T):
            τ = τ + np.diagonal(𝐗[:,:,t].conj().T@i𝚺@𝐗[:,:,t])

        # Computing expression of the estimator
        𝚺_new = 0
        for t in range(0, T):
            𝐗_bis = 𝐗[:,:,t] / np.sqrt(τ)
            𝚺_new = 𝚺_new + (p/N) * 𝐗_bis@𝐗_bis.conj().T

        # Imposing trace constraint: Tr(𝚺) = p
        𝚺_new = p*𝚺_new/np.trace(𝚺_new)

        # Condition for stopping
        δ = np.linalg.norm(𝚺_new - 𝚺, 'fro') / np.linalg.norm(𝚺, 'fro')

        # Updating 𝚺
        𝚺 = 𝚺_new
        iteration = iteration + 1

    if iteration == iter_max:
        warnings.warn('Recursive algorithm did not converge')

    return (𝚺, δ, iteration)


def tyler_estimator_covariance_text(𝐗, tol=0.0001, iter_max=20):
    """ A function that computes the Modified Tyler Fixed Point Estimator for
    covariance matrix estimation under problem TextGen.
        Inputs:
            * 𝐗 = a matrix of size p*N*T with each saptial observation along column dimension and time
                observation along third dimension.
            * tol = tolerance for convergence of estimator
            * iter_max = number of maximum iterations
        Outputs:
            * 𝚺 = array of size (p,p,T) to the different estimates
            * 𝛅 = the final distance between two iterations for each estimate
            * iteration = number of iterations til convergence """

    (p, N, T) = 𝐗.shape
    𝛅 = np.inf*np.ones(T) # Distance between two iterations for each t
    𝚺 = np.tile(np.eye(p).reshape(p,p,1), (1,1,T)) # Initialise all estimates to identity
    iteration = 0

    # Recursive algorithm
    while (np.max(𝛅)>tol) and iteration < iter_max:

        # Compute the textures for each pixel using all the dates avalaibe
        τ = 0
        for t in range(0, T):
            i𝚺_t = np.linalg.inv(𝚺[:,:,t])
            τ = τ + np.diagonal(𝐗[:,:,t].conj().T@i𝚺_t@𝐗[:,:,t])

        # Computing expression of the estimator
        𝚺_new = np.zeros((p,p,T)).astype(complex)
        for t in range(0, T):
            𝐗_bis = 𝐗[:,:,t] / np.sqrt(τ)
            𝚺_new[:,:,t] = (T*p/N) * 𝐗_bis@𝐗_bis.conj().T

            # Imposing trace constraint: Tr(𝚺) = p
            𝚺_new[:,:,t] = p*𝚺_new[:,:,t]/np.trace(𝚺_new[:,:,t])

            # Condition for stopping
            𝛅[t] = np.linalg.norm(𝚺_new[:,:,t] - 𝚺[:,:,t], 'fro') / \
                     np.linalg.norm(𝚺[:,:,t], 'fro')
        
        # Updating 𝚺
        𝚺 = 𝚺_new
        iteration = iteration + 1

    if iteration == iter_max:
        warnings.warn('Recursive algorithm did not converge')

    return (𝚺, 𝛅, iteration)