##############################################################################
# Some functions to make Monte-Carlo simulations smoother
# Authored by Ammar Mian, 29/10/2018
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
from multiprocessing import Process, Queue
import numpy as np
from generic_functions import *
import time

def wrapper_multivariate_complex_normal_samples(data_args):
    """ A wrapper for the Gaussian data generation function multivariate_complex_normal_samples
        in order to have a generic form for Monte-Carlo function """
    mean, covariance, N, pseudo_covariance = data_args
    return multivariate_complex_normal_samples(mean, covariance, N, pseudo_covariance)


def wrapper_multivariate_complex_t_samples(data_args):
    """ A wrapper for the t distribution data generation function multivariate_complex_t_samples
        in order to have a generic form for Monte-Carlo function """
    mean, covariance, N, df, pseudo_covariance = data_args
    return multivariate_complex_t_samples(mean, covariance, N, df, pseudo_covariance)


def wrapper_multivariate_complex_K_samples(data_args):
    """ A wrapper for the K distribution data generation function wrapper_multivariate_complex_K_samples
        in order to have a generic form for Monte-Carlo function """
    mean, covariance, N, mu, b, pseudo_covariance = data_args
    return multivariate_complex_K_samples(mean, covariance, N, mu, b, pseudo_covariance)


def wrapper_multivariate_complex_Cauchy_samples(data_args):
    """ A wrapper for the Cauchy distribution data generation function multivariate_complex_Cauchy_samples
        in order to have a generic form for Monte-Carlo function """
    mean, covariance, N, mu, b, pseudo_covariance = data_args
    return multivariate_complex_Cauchy_samples(mean, covariance, N, mu, b, pseudo_covariance)


def wrapper_multivariate_complex_Laplace_samples(data_args):
    """ A wrapper for the Laplace distribution data generation function multivariate_complex_Laplace_samples
        in order to have a generic form for Monte-Carlo function """
    mean, covariance, N, beta, pseudo_covariance = data_args
    return multivariate_complex_Laplace_samples(mean, covariance, N, beta, pseudo_covariance)


def generate_time_series_multivariate_vector(Args):
    """ A function to generate a time series of random samples of dimension p, where at each date
        there is N independent observations having the same set of parameters. The length of the series is T.
        Inputs:
            * Args compromising:
                * p = dimension of samples
                * N = number of samples at each date sharing the same distribution parameters
                * T = length of time series
                * data_generation_function = a function to generate random samples for each date:
                    must generate an array of shape (p, N)
                * data_args = list of arguments corresponding to the values of the parameters of 
                the random distribution at each date
        Outputs:
            * an array of shape (p, N, T) corresponding to the time series"""
    p, N, T, generation_function, data_args_list = Args
    X = np.zeros((p, N, T)).astype(complex)
    for t in range(0, T):
        X[:, :, t] = generation_function(data_args_list[t])
    return X


def compute_several_statistics(X, Args):
    """ A function to compute and stack the results of several statistics on data X.
        Inputs:
            * X = the data
            * Args = list constitued of statistics_list and statistics_args
        Outputs:
            * a list correspond to the value of each statistic on the data X. """

    statistics_list, statistics_args = Args
    λ = []
    for i_statistic, statistic in enumerate(statistics_list):
        λ.append(statistic(X, statistics_args[i_statistic]))
    return λ 


def compute_monte_carlo(data_generation_function, data_generation_args, function_to_compute,
                        function_args, number_of_trials, multi=False, queue=0, jobno=0):
    """ A function that allowing to compute Monte-Carlo trials by generating random data and computing some
        function of these observations
        Inputs:
            * data_generation_function = a function to generate the random data
            * data_generation_args = arguments to pass to data_generation_function
            * function_to_compute = a function to compute the desired quantity
            * function_args = arguments to pass to function_to_compute
            * number_of_trials = number of Monte-Carlo trials to run
            * multi=False, queue=0, jobno=0 -> serves to parallelize
        Outputs:
            * a list containing the results at each trial """

    # To have a different seed for each job
    if multi:
        np.random.seed(int(time.time())+jobno)

    # Results container
    results = []
    for trial in np.arange(0, number_of_trials):

        # Generate Data
        𝐗 = data_generation_function(data_generation_args)

        # Compute the function of the observations
        result_at_this_trial = function_to_compute(𝐗, function_args)
        results.append(result_at_this_trial)

    if multi:
        queue.put(results)
    else:
        return results


def compute_monte_carlo_parallel(data_generation_function, data_generation_args, function_to_compute,
                                 function_args, number_of_trials, multi=False, number_of_threads=4):
    """ A function that is a prallelisation of compute_monte_carlo
        Inputs:
            * data_generation_function = a function to generate the random data
            * data_generation_args = arguments to pass to data_generation_function
            * function_to_compute = a function to compute the desired quantity
            * function_args = arguments to pass to function_to_compute
            * number_of_trials = number of Monte-Carlo trials to run 
                    (multiple of number_of_threads please 
                    or make sure that number_of_trials/number_of_threads is an 
                    integer at least)
            * multi = True if parallel computing, False if not
            * number_of_threads = number of thread to use (number of cores of the machine in general)
        Outputs:
            *  a list containing the results at each trial """

    if multi:
        results = [] # Results container
        queues = [Queue() for i in range(number_of_threads)] # Serves to obtain result for each thread
        args = [(data_generation_function, data_generation_args, function_to_compute,
                                 function_args, int(number_of_trials/number_of_threads), 
                                 True, queues[i], i) for i in range(number_of_threads)] 
        jobs = [Process(target=compute_monte_carlo, args=a) for a in args]
        # Starting parallel computation
        for j in jobs: j.start()
        # Obtaining result for each thread
        for q in queues: results = results + q.get()
        # Waiting for each thread to terminate
        for j in jobs: j.join()

    else:
        results = compute_monte_carlo(data_generation_function, data_generation_args, 
                    function_to_compute, function_args, number_of_trials, multi=False)
    return results
