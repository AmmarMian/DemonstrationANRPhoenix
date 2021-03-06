##############################################################################
# Plotting CFARness property of considered statistics: only K-dist
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
from generic_functions import *
from change_detection_functions import *
import seaborn as sns
import matplotlib.pyplot as plt
from monte_carlo_tools import *
from tqdm import tqdm


def covariance_cfar():
    #########################################################################################
    # Covariance CFARness
    #########################################################################################
    print("Computing Monte-Carlo simulation for testing the covariance CFAR property")

    # Select zero-mean, no pseudo-covariance, K-distribution
    data_generation_function = generate_time_series_multivariate_vector
    generation_function = wrapper_multivariate_complex_K_samples
    b = 0.3 # Scale of K-distribution texture
    mu = 0.1 # Shape of K-distribution texture
    𝛍 = np.zeros(p)
    pseudo_𝚺 = 0  
    # Values for Toeplitz matrices to compute the distribution for
    ρ_vec = np.linspace(0.01,0.99,3)

    # Container to store the values of the statistics
    λ = np.zeros((number_of_trials, len(statistics_list), len(ρ_vec)))

    i_ρ = 0
    for ρ in tqdm(ρ_vec): # Iterate for each value of ρ
        
        # Generating Covariance matrix with given ρ
        𝚺 = ToeplitzMatrix(ρ, p)

        # Generating parameters to pass to the Monte-Carlo function: since
        # the series is homogeneous, we repeat T times
        data_args_list = [[𝛍, 𝚺, N, mu, b, pseudo_𝚺]]*T 
        data_generation_args = [p, N, T, generation_function, data_args_list]

        # Computing the Monte-Carlo simulation for this value of ρ
        λ[:,:,i_ρ] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                        function_to_compute, function_args, 
                                        number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))

        i_ρ = i_ρ + 1

    # Plotting
    for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
        plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
        ax = plt.gca()
        for i_ρ, ρ in enumerate(ρ_vec): # Plotting the different histograms for each value of ρ
            # Sometimes infinite value can appear due to overflow, we skip them
            λ_to_plot = λ[:,i_Λ, i_ρ]
            λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
            if statistics_scale[i_Λ] == 'log':
                plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                         label=r'$\rho=%.2f$' % ρ, alpha=0.5, density=True,
                         edgecolor='black', linewidth=0.5)
            else:
                plt.hist(λ_to_plot, number_of_bins_histogram, 
                         label=r'$\rho=%.2f$' % ρ, alpha=0.5, density=True,
                         edgecolor='black', linewidth=0.5)
        plt.legend()
        if statistics_scale[i_Λ] == 'log':
            plt.xlabel(r'$\log(\lambda)$')
        else:
            plt.xlabel(r'$\lambda$')
        plt.ylabel(r'PDF')
        plt.title(r'Martix CFAR property for %s' % Λ)
    plt.show()
    
    return λ


def texture_non_equality_cfar():
    #########################################################################################
    # Texture CFARness: The textures are not equal between the dates
    #########################################################################################
   print("Computing Monte-Carlo simulation for testing the texture CFAR property")
   print("The textures is are equal between the dates here")

   # Select zero-mean, no pseudo-covariance, K-distribution
   data_generation_function = generate_time_series_multivariate_vector
   generation_function = wrapper_multivariate_complex_K_samples
   𝚺 = ToeplitzMatrix(0.5, p)
   𝛍 = np.zeros(p)
   pseudo_𝚺 = 0  
   # Scale parameter of for texture in K-distribution
   b = 0.5
   # Values for K-distribution shape parameter to compute the distribution for
   mu_vec = np.linspace(0.1,10,5) 

   # Container to store the values of the statistics
   λ = np.zeros((number_of_trials, len(statistics_list), len(mu_vec)))

   i_mu = 0
   for mu in tqdm(mu_vec): # Iterate for each value of b

       # Generating parameters to pass to the Monte-Carlo function: since
       # the series is homogeneous, we repeat T times
       data_args_list = [[𝛍, 𝚺, N, mu, b, pseudo_𝚺]]*T 
       data_generation_args = [p, N, T, generation_function, data_args_list]

       # Computing the Monte-Carlo simulation for this value of mu
       λ[:,:,i_mu] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                       function_to_compute, function_args, 
                                       number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))

       i_mu = i_mu + 1

   # Plotting
   for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
       plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
       ax = plt.gca()
       for i_mu, mu in enumerate(mu_vec): # Plotting the different histograms for each value of b
            # Sometimes infinite value can appear due to overflow, we skip them
            λ_to_plot = λ[:,i_Λ, i_mu]
            λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
            if statistics_scale[i_Λ] == 'log':
               plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                        label=r'$\mu=%.2f$' % mu, alpha=0.5, density=True,
                        edgecolor='black', linewidth=0.5)
            else:
               plt.hist(λ_to_plot, number_of_bins_histogram, 
                        label=r'$\mu=%.2f$' % mu, alpha=0.5, density=True,
                        edgecolor='black', linewidth=0.5)
       plt.legend()
       if statistics_scale[i_Λ] == 'log':
           plt.xlabel(r'$\log(\lambda)$')
       else:
           plt.xlabel(r'$\lambda$')
       plt.ylabel(r'PDF')
       plt.title(r'Texture CFAR property for %s. The textures are not equal between the dates' % Λ)
   plt.show()
   
   return λ

def texture_equality_cfar():
    #########################################################################################
    # Texture CFARness: The textures are equal between the dates
    #########################################################################################
   print("Computing Monte-Carlo simulation for testing the texture CFAR property")
   print("The textures are equal between the dates here")

   # Select zero-mean, no pseudo-covariance, K-distribution
   data_generation_function = generate_time_series_K_distribution_texture_equality
   𝚺 = ToeplitzMatrix(0.5, p)
   𝛍 = np.zeros(p)
   # Scale parameter of for texture in K-distribution
   b = 0.5
   # Values for K-distribution shape parameter to compute the distribution for
   mu_vec = np.linspace(0.1,10,5) 

   # Container to store the values of the statistics
   λ = np.zeros((number_of_trials, len(statistics_list), len(mu_vec)))

   i_mu = 0
   for mu in tqdm(mu_vec): # Iterate for each value of b

        # Setting arguments to pass to data_generation_function
        data_generation_args = [p, N, T, mu, b, N, 𝛍, 𝚺]

        # Computing the Monte-Carlo simulation for this value of mu
        λ[:,:,i_mu] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                       function_to_compute, function_args, 
                                       number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))

        i_mu = i_mu + 1

   # Plotting
   for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
       plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
       ax = plt.gca()
       for i_mu, mu in enumerate(mu_vec): # Plotting the different histograms for each value of b
            # Sometimes infinite value can appear due to overflow, we skip them
            λ_to_plot = λ[:,i_Λ, i_mu]
            λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
            if statistics_scale[i_Λ] == 'log':
               plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                        label=r'$\mu=%.2f$' % mu, alpha=0.5, density=True,
                        edgecolor='black', linewidth=0.5)
            else:
               plt.hist(λ_to_plot, number_of_bins_histogram, 
                        label=r'$\mu=%.2f$' % mu, alpha=0.5, density=True,
                        edgecolor='black', linewidth=0.5)
       plt.legend()
       if statistics_scale[i_Λ] == 'log':
           plt.xlabel(r'$\log(\lambda)$')
       else:
           plt.xlabel(r'$\lambda$')
       plt.ylabel(r'PDF')
       plt.title(r'Texture CFAR property for %s. The textures are equal between the dates' % Λ)
   plt.show()
   
   return λ


def generate_time_series_K_distribution_texture_equality(Args):
    """ A specific function aimed at geerating K-distributed samples with texture equality constraint.
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


if __name__ == '__main__':

    #########################################################################################
    # Simulation parameters
    #########################################################################################
    # Activate latex in figures (or not)
    latex_in_figures = False
    if latex_in_figures:
      enable_latex_infigures()

    # Enable parallel processing (or not)
    # In general the optimal parameters are obtained for 
    # number_of_threads = number of cores on the machine
    enable_multi = True
    number_of_threads = 48 # for parallel compuatation

    # General parameters
    p = 10
    N = 25
    T = 3

    # Monte-Carlo parameters
    number_of_trials = 1200
    # Statistics to use
    statistics_list = [covariance_equality_glrt_gaussian_statistic,
                        covariance_equality_t1_gaussian_statistic,
                        covariance_equality_Wald_gaussian_statistic,
                        shape_equality_robust_statistic,
                        scale_equality_robust_statistic,
                        scale_and_shape_equality_robust_statistic]
    
    # The Gaussian statistics need no arguments while the robust ones need 
    # convergence criterion and number max of iterations
    statistics_args = [None]*3+[(0.01, 100)]*3

    # Naming the statistics for plotting
    statistics_names = [r'Gaussian GLRT statistic', r'$t_1$ statistic',
                        r'Wald statistic',
                        r'Robust Shape statistic',
                        r'Robust Scale statistic',
                        r'Robust Scale and/or Shape statistic']

    # If the scale of the statistic is to be plotted in logarithmic or linear
    statistics_scale = ['log', 'linear', 'linear', 'log', 'log', 'log']
    
    # Setting a function in order to compute all the statistics for each trial
    function_to_compute = compute_several_statistics
    function_args = [statistics_list, statistics_args]

    # Plotting options
    number_of_bins_histogram = 100
    
    # Printing some information about simulation
    print("This simulation is aimed at testing CFAR properties of various statistics.")
    print("Parameters of the simulation:")
    print("    * p=%d, N=%d, T=%d" % (p,N,T))
    print("    * K distribution is selected")
    print("    * %d Monte-Carlo Trials will be done" % number_of_trials)
    if enable_multi:
        print("    * The simulation will be done using %d threads in parallel" % number_of_threads)
    print("The following statistics have been selected (Name: Arguments):")
    for i_Λ, Λ in enumerate(statistics_names):
        print("    * %s: %s" % (Λ, str(statistics_args[i_Λ])))

    # Letting the choice to test Covariance, texture with non-equality and texture with equality
    choice = 3
    while (choice>2) or (choice<0):
        print("Please chose an option:")
        print("0: Covariance CFAR")
        print("1: Texture CFAR with no equality between dates")
        print("2: Texture CFAR with equality between dates")
        choice = int(input())

    # Doing Simulation
    possibilities = {0 : covariance_cfar,
       1 : texture_non_equality_cfar,
       2 : texture_equality_cfar}
    possibilities[choice]()

