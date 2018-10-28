##############################################################################
# Plotting ROC curves for change detection in a time series under Gaussian model
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
from ChangeDetectionFunctions import *

def compute_monte_carlo_Gaussian(p, N, T, Sigma_t, statistics_list, statistics_args, generation_args=None, ChunkSize=1, multi=False, queue=0, jobno=0):
    """ A function that will generate a Gaussian time series and compute statistics for change detection.
        Inputs:
            * p = dimension of vectors
            * N = number of Samples at each date
            * T = length of time series
            * Sigma_t = numpy array of size (p, p, T) corresponding to the covariance matrix at each date
            * statistics_list = list of function pointers to statistics of decision
            * statistics_args = list of tuples containing arguments to pass to the statistics
            * generation_args = unused here
            * ChunkSize = number of monte-carlo trials to run
            * multi=False, queue=0, jobno=0 -> serves to parallelize
        Outputs:
            * an array of size (len(statistics_list), ChunkSize) 
              containing the results at each trial for each statistic"""

    # To have a different seed for each job
    if multi:
        np.random.seed(jobno)

    # Results container
    statistics_results = np.zeros((ChunkSize, len(statistics_list)))
    for trial in tqdm(np.arange(0, ChunkSize)):

        # Generate Time Series
        X = np.zeros((p, N, T)).astype(complex)
        for t in range(0, T):
            X[:, :, t] = multivariate_complex_normal_samples(np.zeros(p), Sigma_t[:,:,t], N)

        # Compute each statistic
        for i_s in range(0, len(statistics_list)):
            statistic_function = statistics_list[i_s]
            args = statistics_args[i_s]
            statistics_results[trial, i_s] = statistic_function(X, *args)

    if multi:
        queue.put(statistics_results)
    else:
        return statistics_results


def compute_monte_carlo_parallel(p, N, T, Sigma_t, statistics_list, statistics_args, monte_carlo_function, generation_args=None, numberOfTrials=100, multi=False, CORES=4):
    """ Parallel implementation of compute_monte_carlo_xx functions where xx is either Gaussian, K, t or Laplace.
        monte_carlo_function is the function to compute in parallel. """
    if multi:
        statistics_results = []
        queues = [Queue() for i in range(CORES)]
        args = [(p, N, T, Sigma_t, statistics_list, statistics_args, generation_args, int(numberOfTrials/CORES), True, queues[i], i) for i in range(CORES)]
        jobs = [Process(target=monte_carlo_function, args=a) for a in args]
        for j in jobs: j.start()
        for q in tqdm(queues): statistics_results.append(q.get())
        for j in jobs: j.join()
        statistics_results = np.vstack(MSE)

    else:
        statistics_results = monte_carlo_function(p, N, T, Sigma_t, statistics_list, statistics_args, generation_args=generation_args, ChunkSize=numberOfTrials)
    return statistics_results


if __name__ == '__main__':

    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib import rc

    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    sns.set_style("darkgrid")
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"

    #########################################################################################
    # Simulation parameters
    #########################################################################################

    # General parameters
    p = 3
    N = 10
    T = 5

    # Change parameters
    t_C = 3 # change-point 
    Sigma_t = np.zeros((p,p,T)).astype(complex)
    rho_before = 0.1
    rho_after = 0.9
    for t in range(0, t_C):
        Sigma_t[:,:,t] = ToeplitzMatrix(rho_before, p)
    for t in range(t_C, T):
        Sigma_t[:,:,t] = ToeplitzMatrix(rho_after, p)

    # Monte-Carlo parameters
    numberOfTrials_pfa = 10000
    numberOfTrials_pd = 1000

    # Plot parameters
    numberOfPoints = 30
    Pfa = np.logspace(0, np.log10(100/numberOfTrials_pfa), numberOfPoints) # Pfa vector for plotting

    # Statistics to use
    statistics_list = [covariance_equality_glrt_gaussian_statistic, 
                        covariance_equality_t1_gaussian_statistic,
                        covariance_equality_Wald_gaussian_statistic]
    statistics_args = [(0,0), (0,0), (0,0)]
    statistics_names = [r'Gaussian GLRT statistic', r'$t_1$ statistic', r'Wald statistic']

    #########################################################################################
    # Computing Monte-Carlo trials in a Gaussian model to have the relation Pfa-treshold 
    # for each statistic
    #########################################################################################
    print('Computing PFa-thresholds relationship numerically and selectionning tresholds:')
    Pfa_temp = np.linspace(1, 1/numberOfTrials_pfa, numberOfTrials_pfa)
    rho = rho_before
    Sigma_t_pfa = np.zeros((p,p,T)).astype(complex)
    for t in range(0,T):
        Sigma_t_pfa[:,:,t] = ToeplitzMatrix(rho, p)
    thresholds_temp = compute_monte_carlo_Gaussian(p, N, T, Sigma_t_pfa, statistics_list, statistics_args, ChunkSize=numberOfTrials_pfa)
    # Obtaining a threshold for each point of the ROC curve by doin an interpolation from the results of the monte-carlo
    thresholds = np.zeros((numberOfPoints, len(statistics_list)))
    for i_pfa in range(0, len(Pfa)):
        pfa = Pfa[i_pfa]
        for i_s in range(0, len(statistics_list)):
            thresholds[i_pfa, i_s] = np.interp(pfa, Pfa_temp[::-1], np.sort(thresholds_temp[:, i_s])[::-1])

    #########################################################################################
    # Computing Monte-Carlo trials in Gaussian model to have the performance of detection
    #########################################################################################
    print('\nComputing ROC curves points:\n')
    Pd = np.zeros((numberOfPoints, len(statistics_list)))
    for i_pfa in range(0, len(Pfa)):
        print('Computing iteration %d of %d' % (i_pfa, len(Pfa)))
        results = compute_monte_carlo_Gaussian(p, N, T, Sigma_t, statistics_list, statistics_args, ChunkSize=numberOfTrials_pd)
        for i_s in range(0, len(statistics_list)):
            thresh = thresholds[i_pfa, i_s]
            Pd[i_pfa, i_s] = np.mean(results[:,i_s]>=thresh)

    #########################################################################################
    # Plotting
    #########################################################################################
    plt.figure()
    markers = ['o', 'd', 's', '*', '+']
    ax = plt.gca()
    for i_s in range(0, len(statistics_list)):
        ax.scatter(Pfa, Pd[:, i_s], label=statistics_names[i_s], marker=markers[i_s])
    plt.legend()
    ax.set_yscale('log')
    ax.set_xscale('log')