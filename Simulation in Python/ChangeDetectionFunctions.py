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
from Functions import *

##############################################################################
# Gaussian Statistics
##############################################################################
def covariance_equality_glrt_gaussian_statistic(X, *args):
    """ GLRT statistic for detecting a change of covariance matrix in a multivariate Gaussian Time Series.
        At each time, Ni.i.d samples are available. A description of the statistic can be found in:
        D. Ciuonzo, V. Carotenuto and A. De Maio, 
        "On Multiple Covariance Equality Testing with Application to SAR Change Detection," 
        in IEEE Transactions on Signal Processing, vol. 65, no. 19, pp. 5078-5091, 1 Oct.1, 2017.
        doi: 10.1109/TSP.2017.2712124
        Inputs:
            * X = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
        Outputs:
            * the GLRT statistic given the observations in input"""

    (p, N, T) = X.shape
    S = 0
    logDenominator = 0
    for t in range(0, T):
        St = SCM(X[:, :, t])
        logDenominator = logDenominator + N * np.log(np.abs(np.linalg.det(St)))
        S = S + St / T
    logNumerator = N * T * np.log(np.abs(np.linalg.det(S)))
    return np.exp(logNumerator - logDenominator)


def covariance_equality_t1_gaussian_statistic(X, *args):
    """ t1 statistic for detecting a change of covariance matrix in a multivariate Gaussian Time Series.
        At each time, Ni.i.d samples are available. A description of the statistic can be found in:
        D. Ciuonzo, V. Carotenuto and A. De Maio, 
        "On Multiple Covariance Equality Testing with Application to SAR Change Detection," 
        in IEEE Transactions on Signal Processing, vol. 65, no. 19, pp. 5078-5091, 1 Oct.1, 2017.
        doi: 10.1109/TSP.2017.2712124
        Inputs:
            * X = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
        Outputs:
            * the t1 statistic given the observations in input"""

    (N, K, M) = X.shape

    Sigma_10 = SCM(X.reshape((N, K*M)))
    iSigma_10 = np.linalg.inv(Sigma_10)
    t1 = 0
    for t in range(0, M):
        Sigma_m1 = SCM(X[:, :, t])
        S = (iSigma_10 @ Sigma_m1)
        t1 = t1 + np.trace( S @ S )/M;
    return np.real(t1)


def covariance_equality_Wald_gaussian_statistic(X, *args):
    """ Wald statistic for detecting a change of covariance matrix in a multivariate Gaussian Time Series.
        At each time, Ni.i.d samples are available. A description of the statistic can be found in:
        D. Ciuonzo, V. Carotenuto and A. De Maio, 
        "On Multiple Covariance Equality Testing with Application to SAR Change Detection," 
        in IEEE Transactions on Signal Processing, vol. 65, no. 19, pp. 5078-5091, 1 Oct.1, 2017.
        doi: 10.1109/TSP.2017.2712124
        Inputs:
            * X = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
        Outputs:
            * the Wald statistic given the observations in input"""

    (N, K, M) = X.shape
    L = 0;
    O = 0;
    Q = 0;
    Sigma_11 = SCM(X[:, :, 0])
    for m in range(0,M):
        Sigma_m1 = SCM(X[:, :, m])
        iSigma_m1 = np.linalg.inv(Sigma_m1)
        if m != 0:
            S = np.eye(N) - Sigma_11@iSigma_m1
            L = L + K*np.trace(S@S)
            Q = Q + K*(iSigma_m1 - iSigma_m1@Sigma_11@iSigma_m1)
        O = O + K*np.kron(iSigma_m1.T, iSigma_m1)
    
    return np.real(L - vec(Q).conj().T @ (np.linalg.inv(O)@vec(Q)))


##############################################################################
# Robust Statistics
##############################################################################
def scale_and_shape_equality_robust_statistic(X, args):
    """ GLRT test for testing a change in the scale or/and shape in 
        a deterministic SIRV model.
        Inputs:
            * X = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * args = tol, iterMax for Tyler
        Outputs:
            * the statistic given the observations in input"""
    tol, iterMax = args
    (p, N, T) = X.shape
    (Sigma_0, err, niter) = TylerFixedPointMatAndText(X, tol, iterMax)
    iSigma_0  = np.linalg.inv(Sigma_0)
    logDenominator = T*p*np.log(T)
    logNumerator = T*N*np.log(np.abs(np.linalg.det(Sigma_0)))
    numTemp = 0
    for t in range(0, T):
        xkt = X[:, :, t]
        (Sigma_t, err, niter) = TylerFixedPoint(xkt, tol, iterMax)
        numTemp = numTemp + np.diagonal(xkt.conj().T @ iSigma_0 @ xkt)
        logDenominator = logDenominator + N*np.log(np.abs(np.linalg.det(Sigma_t))) + \
                         p*np.sum(np.log(np.diagonal(xkt.conj().T @ np.linalg.inv(Sigma_t) @ xkt)))
    logNumerator = logNumerator + T*p*np.sum(np.log(numTemp))
    return np.exp(np.real(logNumerator - logDenominator))

def shape_equality_robust_statistic(X, args):
    """ GLRT test for testing a change in the shape in 
        a deterministic SIRV model.
        Inputs:
            * X = a (p, N, T) numpy array with:
                * p = dimension of vectors
                * N = number of Samples at each date
                * T = length of time series
            * args = tol, iterMax for Tyler
        Outputs:
            * the statistic given the observations in input"""

    tol, iterMax = args
    (p, N, T) = X.shape
    (Sigma_0, err, niter) = TylerFixedPoint(np.reshape(X, (p, N*T)), tol, iterMax)
    iSigma_0 = np.linalg.inv(Sigma_0)
    logDenominator = 0
    logNumerator = T * N * np.log(np.abs(np.linalg.det(Sigma_0)))
    lognumTemp = 0
    for t in range(0, T):
        xkt = X[:, :, t]
        (Sigma_t, err, niter) = TylerFixedPoint(xkt, tol, IterMax)
        lognumTemp = lognumTemp + p * np.sum(np.log(np.diagonal(xkt.conj().T @ iSigma_0 @ xkt)))
        logDenominator = logDenominator + N * np.log(np.abs(np.linalg.det(Sigma_t))) + \
                         p * np.sum(np.log(np.diagonal(xkt.conj().T @ np.linalg.inv(Sigma_t) @ xkt)))
    logNumerator = logNumerator + lognumTemp
    return np.real(logNumerator - logDenominator)

##############################################################################
# Some Functions
##############################################################################
def TylerFixedPointMatAndText( X, tol=0.0001, iterMax=20):
    """ A function that computes the Modified Tyler Fixed Point Estimator for covariance matrix estimation
        under problem MatAndText.
        Inputs:
            * X = a matrix of size p*N*T with each saptial observation along column dimension and time
                observation along third dimension.
            * tol = tolerance for convergence of estimator
            * iterMax = number of maximum iterations
        Outputs:
            * Sigma = the estimate
            * error = the final error between two iterations
            * iteration = number of iterations til convergence """

    (p, N, T) = X.shape
    error = np.inf
    converged = False
    iteration = 0
    Sigma = np.eye(p)
    while not converged and iteration < iterMax:
        iteration = iteration + 1

        # Compute the textures estimate for each pixel
        a = 0
        R = np.linalg.cholesky(Sigma).conj().T
        iR = np.linalg.inv(R)
        for t in range(0, T):
            x = X[:, :, t]
            vx = np.dot(np.linalg.inv(iR), x)
            ax = np.mean(vx * np.conj(vx), axis=0)
            a = a + ax
        a = np.tile(a, (p, T))
        # Compute the estimate for next iteration
        x = np.reshape(X, (p, N*T))
        xbis = x / np.sqrt(a)
        Sigma_new = np.dot(xbis, xbis.conj().T) / N
        Sigma_new = p * Sigma_new / np.trace(Sigma_new)

        # Compute error
        error = np.linalg.norm(Sigma_new - Sigma, 'fro') / np.linalg.norm(Sigma, 'fro')
        converged = error < tol
        Sigma = Sigma_new

    if iteration == iterMax:
        print("Waring: TylerMT Fixed point algorithm has not converged")

    return (Sigma, error, iteration)