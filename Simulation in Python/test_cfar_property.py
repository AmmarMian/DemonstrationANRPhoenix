##############################################################################
# Plotting CFARness property of considered statistics
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
import numpy as np
import scipy as sp
from Functions import *
from ChangeDetectionFunctions import *
from bunny import bunny
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from monte_carlo_tools import *

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
sns.set_style("darkgrid")
__spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"


if __name__ == '__main__':

    #########################################################################################
    # Simulation parameters
    #########################################################################################

    # General parameters
    p = 3
    N = 10
    T = 2

    # Monte-Carlo parameters
    number_of_trials = 10000

    # Statistics to use
    statistics_list = [covariance_equality_glrt_gaussian_statistic,
                        covariance_equality_t1_gaussian_statistic,
                        covariance_equality_Wald_gaussian_statistic,
                        scale_and_shape_equality_robust_statistic]
    statistics_args = [None]*3+[(0.001, 30)]
    statistics_names = [r'Gaussian GLRT statistic', r'$t_1$ statistic',
                        r'Wald statistic', " RobustScale And Shape statistic"]
    statistics_scale = ['log', 'linear', 'linear', 'log']


    #########################################################################################
    # Covariance CFARness under Gaussian model
    #########################################################################################
    print("Computing Monte-Carlo simulation for testing the covariance CFAR property")

    # Select zero-mean, no pseudo-covariance, t-distribution
    data_generation_function = generate_time_series_multivariate_vector
    generation_function = wrapper_multivariate_complex_t_samples
    df = 1
    𝛍 = np.zeros(p)
    pseudo_𝚺 = 0  

    # Setting a function in order to compute all the statistics for each trial
    function_to_compute = compute_several_statistics
    function_args = [statistics_list, statistics_args]

    # Values for Toeplitz matrices to compute the distribution for
    ρ_vec = np.linspace(0.1,0.99,3)

    # Container to store the values of the statistics
    λ = np.zeros((number_of_trials, len(statistics_list), len(ρ_vec)))

    i_ρ = 0
    for ρ in bunny(ρ_vec): # Iterate for each value of ρ
        
        # Generating Covariance matrix with given ρ
        𝚺 = ToeplitzMatrix(ρ, p)

        # Generating parameters to pass to the Monte-Carlo function: since
        # the series is homogeneous, we repeat T times
        data_args_list = [[𝛍, 𝚺, N, df, pseudo_𝚺]]*T 
        data_generation_args = [p, N, T, generation_function, data_args_list]

        # Computing the Monte-Carlo simulation for this value of ρ
        λ[:,:,i_ρ] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                        function_to_compute, function_args, 
                                        number_of_trials, multi=True, number_of_threads=4))

        i_ρ = i_ρ + 1

    # Plotting
    for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
        plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
        ax = plt.gca()
        for i_ρ, ρ in enumerate(ρ_vec): # Plotting the different histograms for each value of ρ
            if statistics_scale[i_Λ] == 'log':
                plt.hist(np.log(λ[:,i_Λ, i_ρ]), label=r'$\rho=%.2f$' % ρ, alpha=0.9)
            else:
                plt.hist(λ[:,i_Λ, i_ρ], label=r'$\rho=%.2f$' % ρ, alpha=0.9)
        plt.legend()
        if statistics_scale[i_Λ] == 'log':
            plt.xlabel(r'$\log(\lambda)$')
        else:
            plt.xlabel(r'$\lambda$')
        plt.ylabel(r'Count')
        plt.title(r'CFARness property for %s' % Λ)
    plt.show()