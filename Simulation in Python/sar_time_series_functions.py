##############################################################################
# Some functions to generate synthetic data to model SAR time series
# Authored by Ammar Mian, 04/11/2018
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
from generic_functions import *


def generate_time_series_Gaussian_distribution_texture_equality(Args):
    """ A specific function aimed at generating Gaussian-distributed samples with texture equality constraint.
            * Args compromising:
                * p = dimension of samples
                * N = number of samples at each date sharing the same distribution parameters
                * T = length of time series
                * N = number of Samples
                * 𝛍 = mean
                * 𝚺 = covariance matrix
        Outputs:
            * an array of shape (p, N, T) corresponding to the time series"""

    p, N, T, 𝛍, 𝚺 = Args
    𝐗 = np.zeros((p, N, T)).astype(complex)
    for t in range(0, T):
        𝐗[:, :, t] = multivariate_complex_normal_samples(𝛍, 𝚺, N)  
    return 𝐗


def generate_time_series_K_distribution_texture_equality(Args):
    """ A specific function aimed at generating K-distributed samples with texture equality constraint.
            * Args compromising:
                * p = dimension of samples
                * N = number of samples at each date sharing the same distribution parameters
                * T = length of time series
                * mu = Shape parameter
                * b = Scale parameter
                * N = number of Samples
                * 𝛍 = mean
                * 𝚺 = covariance matrix
        Outputs:
            * an array of shape (p, N, T) corresponding to the time series"""

    p, N, T, mu, b, N, 𝛍, 𝚺 = Args
    𝐗 = np.zeros((p, N, T)).astype(complex)
    τ = np.random.gamma(mu, 2/(b**2), N)
    for t in range(0, T):
        𝐳 = multivariate_complex_normal_samples(np.zeros(𝛍.shape), 𝚺, N)
        𝐗[:, :, t] = np.tile(𝛍.reshape((len(𝛍),1)),(1,N)) + 𝐳*np.sqrt(τ)[None,:]  
    return 𝐗


def generate_time_series_t_distribution_texture_equality(Args):
    """ A specific function aimed at generating t-distributed samples with texture equality constraint.
            * Args compromising:
                * p = dimension of samples
                * N = number of samples at each date sharing the same distribution parameters
                * T = length of time series
                * df = degrees of freedom of the chi-squared distribution
                * N = number of Samples
                * 𝛍 = mean
                * 𝚺 = covariance matrix
        Outputs:
            * an array of shape (p, N, T) corresponding to the time series"""

    p, N, T, df, N, 𝛍, 𝚺 = Args
    𝐗 = np.zeros((p, N, T)).astype(complex)
    if df == np.inf:
        τ = 1
    else:
        τ = np.random.chisquare(df, N)/df
    for t in range(0, T):
        𝐳 = multivariate_complex_normal_samples(np.zeros(𝛍.shape), 𝚺, N)
        𝐗[:, :, t] = np.tile(𝛍.reshape((len(𝛍),1)),(1,N)) + 𝐳*np.sqrt(τ)[None,:]  
    return 𝐗


def generate_time_series_Cauchy_distribution_texture_equality(Args):
    """ A specific function aimed at generating Cauchy-distributed samples with texture equality constraint.
            * Args compromising:
                * p = dimension of samples
                * N = number of samples at each date sharing the same distribution parameters
                * T = length of time series
                * mu = Shape parameter
                * b = Scale parameter
                * N = number of Samples
                * 𝛍 = mean
                * 𝚺 = covariance matrix
        Outputs:
            * an array of shape (p, N, T) corresponding to the time series"""

    p, N, T, mu, b, N, 𝛍, 𝚺 = Args
    𝐗 = np.zeros((p, N, T)).astype(complex)
    τ = np.random.gamma(mu, 2/(b**2), N)
    for t in range(0, T):
        𝐳 = multivariate_complex_normal_samples(np.zeros(𝛍.shape), 𝚺, N)
        𝐗[:, :, t] = np.tile(𝛍.reshape((len(𝛍),1)),(1,N)) + 𝐳/np.sqrt(τ)[None,:]  
    return 𝐗


def generate_time_series_Laplace_distribution_texture_equality(Args):
    """ A specific function aimed at generating Laplace-distributed samples with texture equality constraint.
            * Args compromising:
                * p = dimension of samples
                * N = number of samples at each date sharing the same distribution parameters
                * T = length of time series
                * b = Scale parameter
                * N = number of Samples
                * 𝛍 = mean
                * 𝚺 = covariance matrix
        Outputs:
            * an array of shape (p, N, T) corresponding to the time series"""

    p, N, T, b, N, 𝛍, 𝚺 = Args
    𝐗 = np.zeros((p, N, T)).astype(complex)
    τ = np.random.exponential(beta, N)
    for t in range(0, T):
        𝐳 = multivariate_complex_normal_samples(np.zeros(𝛍.shape), 𝚺, N)
        𝐗[:, :, t] = np.tile(𝛍.reshape((len(𝛍),1)),(1,N)) + 𝐳*np.sqrt(τ)[None,:]  
    return 𝐗

