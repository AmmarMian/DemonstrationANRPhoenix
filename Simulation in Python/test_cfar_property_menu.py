##############################################################################
# Plotting CFARness property of considered statistics
# In this file, a menu allows to chose some parameters of the simulation
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
from sar_time_series_functions import *
from change_detection_functions import *
import seaborn as sns
sns.set_style("darkgrid")
import matplotlib.pyplot as plt
from matplotlib import rc
from monte_carlo_tools import *
from tqdm import tqdm
import ast



def covariance_cfar():
    #########################################################################################
    # Covariance CFARness
    #########################################################################################
    # Activate latex in figures (or not)
    latex_in_figures = False
    if latex_in_figures:
      enable_latex_infigures()

    # Select zero-mean, no pseudo-covariance, chosen distrbution
    data_generation_function = generate_time_series_multivariate_vector
    𝛍 = np.zeros(p)
    pseudo_𝚺 = 0 

    # Treating the different possibilities of distributions
    data_functions_per_distribution = {0: wrapper_multivariate_complex_normal_samples,
        1: wrapper_multivariate_complex_K_samples,
        2: wrapper_multivariate_complex_t_samples,
        3: wrapper_multivariate_complex_Cauchy_samples,
        4: wrapper_multivariate_complex_Laplace_samples,}
    generation_function = data_functions_per_distribution[choice_distrib]
    if choice_distrib==1:
        print("Enter shape and scale argument for K-distribution (as a list \"[mu, b]\"):")
        mu, b = ast.literal_eval(input())
    if choice_distrib==2:
        print("Enter degree of freedoms argument for t-distribution:")
        df = float(input())
    if choice_distrib==3:
        print("Enter shape and scale argument for Cauchy-distribution (as a list \"[mu, b]\"):")
        mu, b = ast.literal_eval(input())
    if choice_distrib==4:
        print("Enter scale argument for Laplace-distribution:")
        b = float(input())

    # Asking for values of ρ (Toeplitz matrices)
    print("Enter values of rho (the Toeplitz parameter of covariance) as a list:")
    ρ_vec = ast.literal_eval(input())

    # Container to store the values of the statistics
    λ = np.zeros((number_of_trials, len(statistics_list), len(ρ_vec)))

    i_ρ = 0
    print("Computing Monte-Carlo simulation for testing the covariance CFAR property")
    print_bunny()
    for ρ in tqdm(ρ_vec): # Iterate for each value of ρ
        
        # Generating Covariance matrix with given ρ
        𝚺 = ToeplitzMatrix(ρ, p)

        # Generating parameters to pass to the Monte-Carlo function: since
        # the series is homogeneous, we repeat T times
        if choice_distrib==0:
            data_args_list = [[𝛍, 𝚺, N, pseudo_𝚺]]*T 
        if choice_distrib==1:
            data_args_list = [[𝛍, 𝚺, N, mu, b, pseudo_𝚺]]*T 
        if choice_distrib==2:
            data_args_list = [[𝛍, 𝚺, N, df, pseudo_𝚺]]*T 
        if choice_distrib==3:
            data_args_list = [[𝛍, 𝚺, N, mu, b, pseudo_𝚺]]*T 
        if choice_distrib==4:
            data_args_list = [[𝛍, 𝚺, N, b, pseudo_𝚺]]*T     
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

    # Select zero-mean, no pseudo-covariance
    data_generation_function = generate_time_series_multivariate_vector
    𝚺 = ToeplitzMatrix(0.5, p)
    𝛍 = np.zeros(p)
    pseudo_𝚺 = 0

    # Treating the different possibilities of distributions
    data_functions_per_distribution = {0: wrapper_multivariate_complex_normal_samples,
        1: wrapper_multivariate_complex_K_samples,
        2: wrapper_multivariate_complex_t_samples,
        3: wrapper_multivariate_complex_Cauchy_samples,
        4: wrapper_multivariate_complex_Laplace_samples,}
    generation_function = data_functions_per_distribution[choice_distrib]

    if choice_distrib==0:
        print("Texture cfar cannot be tested when data is Gaussian")
        return 0

    # Since K-dist and Cauchy-dist share the same parameters
    if choice_distrib==1 or choice_distrib==3: 
        choice_parameter = None
        while (choice_parameter!=0) and (choice_parameter!=1):
            print("Select which parameter of texture vary:")
            print("0. Scale")
            print("1. Shape")
            choice_parameter = int(input())
            print("choice:%s"%choice_parameter)

        # Case Scale vary
        if choice_parameter==0:
            # Get the parameters from the user
            print("Enter a value for the shape:")
            mu = float(input())
            print("Enter several values for the scale as a list:")
            b_vec = ast.literal_eval(input())

            # Container to store the values of the statistics
            λ = np.zeros((number_of_trials, len(statistics_list), len(b_vec)))

            i_b = 0
            print("Computing Monte-Carlo simulation for testing the texture CFAR property")
            print_bunny()
            for b in tqdm(b_vec): # Iterate for each value of b
               # Generating parameters to pass to the Monte-Carlo function: since
               # the series is homogeneous, we repeat T times
               data_args_list = [[𝛍, 𝚺, N, mu, b, pseudo_𝚺]]*T 
               data_generation_args = [p, N, T, generation_function, data_args_list]

               # Computing the Monte-Carlo simulation for this value of b
               λ[:,:,i_b] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                               function_to_compute, function_args, 
                                               number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))
               i_b = i_b + 1
               
            # Plotting
            for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
               plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
               ax = plt.gca()
               for i_b, b in enumerate(b_vec): # Plotting the different histograms for each value of b
                    # Sometimes infinite value can appear due to overflow, we skip them
                    λ_to_plot = λ[:,i_Λ, i_b]
                    λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
                    if statistics_scale[i_Λ] == 'log':
                       plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
                    else:
                       plt.hist(λ_to_plot, number_of_bins_histogram, 
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
               plt.legend()
               if statistics_scale[i_Λ] == 'log':
                   plt.xlabel(r'$\log(\lambda)$')
               else:
                   plt.xlabel(r'$\lambda$')
               plt.ylabel(r'PDF')
               plt.title(r'Texture CFAR property for %s. The textures are not equal between the dates' % Λ)
            plt.show()

        # Case shape vary
        else:
            print("Enter a value for the scale:")
            b = float(input())
            print("Enter several values for the shape as a list:")
            mu_vec = ast.literal_eval(input())

            # Container to store the values of the statistics
            λ = np.zeros((number_of_trials, len(statistics_list), len(mu_vec)))

            i_mu = 0
            print("Computing Monte-Carlo simulation for testing the texture CFAR property")
            print_bunny()
            for mu in tqdm(mu_vec): # Iterate for each value of mu

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
               for i_mu, mu in enumerate(mu_vec): # Plotting the different histograms for each value of mu
                    # Sometimes infinite value can appear due to overflow, we skip them
                    λ_to_plot = λ[:,i_Λ, i_mu]
                    λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
                    if statistics_scale[i_Λ] == 'log':
                       plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
                    else:
                       plt.hist(λ_to_plot, number_of_bins_histogram, 
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
               plt.legend()
               if statistics_scale[i_Λ] == 'log':
                   plt.xlabel(r'$\log(\lambda)$')
               else:
                   plt.xlabel(r'$\lambda$')
               plt.ylabel(r'PDF')
               plt.title(r'Texture CFAR property for %s. The textures are not equal between the dates' % Λ)
            plt.show()

    # Either t-distrib or Cauchy: they have the same nmber of parameters
    else:
        if choice_distrib==2:
            # Get the parameters from the user
            print("Enter several values for the degree of freedom as a list:")
            df_or_b_vec = ast.literal_eval(input())
            name_param = 'df'
        else:
            # Get the parameters from the user
            print("Enter several values for the scale as a list:")
            df_or_b_vec = ast.literal_eval(input())
            name_param = 'b'

        # Container to store the values of the statistics
        λ = np.zeros((number_of_trials, len(statistics_list), len(df_or_b_vec)))

        i_df_or_b = 0
        print("Computing Monte-Carlo simulation for testing the texture with no equality between dates CFAR property")
        print_bunny()
        for df_or_b in tqdm(df_or_b_vec): # Iterate for each value of df_or_b

           # Generating parameters to pass to the Monte-Carlo function: since
           # the series is homogeneous, we repeat T times
           data_args_list = [[𝛍, 𝚺, N, df_or_b, pseudo_𝚺]]*T 
           data_generation_args = [p, N, T, generation_function, data_args_list]

           # Computing the Monte-Carlo simulation for this value of df_or_b
           λ[:,:,i_df_or_b] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                           function_to_compute, function_args, 
                                           number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))
           i_df_or_b = i_df_or_b + 1

        # Plotting
        for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
           plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
           ax = plt.gca()
           for i_df_or_b, df_or_b in enumerate(df_or_b_vec): # Plotting the different histograms for each value of df_or_b_vec
                # Sometimes infinite value can appear due to overflow, we skip them
                λ_to_plot = λ[:,i_Λ, i_df_or_b]
                λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
                if statistics_scale[i_Λ] == 'log':
                   plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                            label=r'$%s=%.2f$' % (name_param, df_or_b), alpha=0.5, density=True,
                            edgecolor='black', linewidth=0.5)
                else:
                   plt.hist(λ_to_plot, number_of_bins_histogram, 
                            label=r'$%s=%.2f$' % (name_param, df_or_b), alpha=0.5, density=True,
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

    # Select zero-mean, no pseudo-covariance
    𝚺 = ToeplitzMatrix(0.5, p)
    𝛍 = np.zeros(p)

    # Treating the different possibilities of distributions
    data_functions_per_distribution = {0: generate_time_series_Gaussian_distribution_texture_equality,
        1: generate_time_series_K_distribution_texture_equality,
        2: generate_time_series_t_distribution_texture_equality,
        3: generate_time_series_Cauchy_distribution_texture_equality,
        4: generate_time_series_Laplace_distribution_texture_equality,}
    data_generation_function = data_functions_per_distribution[choice_distrib]

    if choice_distrib==0:
        print("Texture cfar cannot be tested when data is Gaussian")
        return 0

    # Since K-dist and Cauchy-dist share the same parameters
    if choice_distrib==1 or choice_distrib==3: 
        choice_parameter = None
        while (choice_parameter!=0) and (choice_parameter!=1):
            print("Select which parameter of texture vary:")
            print("0. Scale")
            print("1. Shape")
            choice_parameter = int(input())
            print("choice:%s"%choice_parameter)

        # Case Scale vary
        if choice_parameter==0:
            # Get the parameters from the user
            print("Enter a value for the shape:")
            mu = float(input())
            print("Enter several values for the scale as a list:")
            b_vec = ast.literal_eval(input())

            # Container to store the values of the statistics
            λ = np.zeros((number_of_trials, len(statistics_list), len(b_vec)))

            i_b = 0
            print("Computing Monte-Carlo simulation for testing the texture CFAR property")
            print_bunny()
            for b in tqdm(b_vec): # Iterate for each value of b

               # Generating parameters to pass to the Monte-Carlo function: 
               data_generation_args = [p, N, T, mu, b, N, 𝛍, 𝚺]

               # Computing the Monte-Carlo simulation for this value of b
               λ[:,:,i_b] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                               function_to_compute, function_args, 
                                               number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))
               i_b = i_b + 1
               
            # Plotting
            for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
               plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
               ax = plt.gca()
               for i_b, b in enumerate(b_vec): # Plotting the different histograms for each value of b
                    # Sometimes infinite value can appear due to overflow, we skip them
                    λ_to_plot = λ[:,i_Λ, i_b]
                    λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
                    if statistics_scale[i_Λ] == 'log':
                       plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
                    else:
                       plt.hist(λ_to_plot, number_of_bins_histogram, 
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
               plt.legend()
               if statistics_scale[i_Λ] == 'log':
                   plt.xlabel(r'$\log(\lambda)$')
               else:
                   plt.xlabel(r'$\lambda$')
               plt.ylabel(r'PDF')
               plt.title(r'Texture CFAR property for %s. The textures are equal between the dates' % Λ)
            plt.show()

        # Case shape vary
        else:
            print("Enter a value for the scale:")
            b = float(input())
            print("Enter several values for the shape as a list:")
            mu_vec = ast.literal_eval(input())

            # Container to store the values of the statistics
            λ = np.zeros((number_of_trials, len(statistics_list), len(mu_vec)))

            i_mu = 0
            print("Computing Monte-Carlo simulation for testing the texture CFAR property")
            print_bunny()
            for mu in tqdm(mu_vec): # Iterate for each value of mu

               # Generating parameters to pass to the Monte-Carlo function: 
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
               for i_mu, mu in enumerate(mu_vec): # Plotting the different histograms for each value of mu
                    # Sometimes infinite value can appear due to overflow, we skip them
                    λ_to_plot = λ[:,i_Λ, i_mu]
                    λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
                    if statistics_scale[i_Λ] == 'log':
                       plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
                    else:
                       plt.hist(λ_to_plot, number_of_bins_histogram, 
                                label=r'$\mu=%.2f, b=%.2f$' % (mu, b), alpha=0.5, density=True,
                                edgecolor='black', linewidth=0.5)
               plt.legend()
               if statistics_scale[i_Λ] == 'log':
                   plt.xlabel(r'$\log(\lambda)$')
               else:
                   plt.xlabel(r'$\lambda$')
               plt.ylabel(r'PDF')
               plt.title(r'Texture CFAR property for %s. The textures are equal between the dates' % Λ)
            plt.show()

    # Either t-distrib or Cauchy: they have the same nmber of parameters
    else:
        if choice_distrib==2:
            # Get the parameters from the user
            print("Enter several values for the degree of freedom as a list:")
            df_or_b_vec = ast.literal_eval(input())
            name_param = 'df'
        else:
            # Get the parameters from the user
            print("Enter several values for the scale as a list:")
            df_or_b_vec = ast.literal_eval(input())
            name_param = 'b'

        # Container to store the values of the statistics
        λ = np.zeros((number_of_trials, len(statistics_list), len(df_or_b_vec)))

        i_df_or_b = 0
        print("Computing Monte-Carlo simulation for testing the texture with equality between dates CFAR property")
        print_bunny()
        for df_or_b in tqdm(df_or_b_vec): # Iterate for each value of df_or_b

           # Generating parameters to pass to the Monte-Carlo function:
           data_generation_args = [p, N, T, df_or_b, N, 𝛍, 𝚺]

           # Computing the Monte-Carlo simulation for this value of df_or_b
           λ[:,:,i_df_or_b] = np.array(compute_monte_carlo_parallel(data_generation_function, data_generation_args, 
                                           function_to_compute, function_args, 
                                           number_of_trials, multi=enable_multi, number_of_threads=number_of_threads))
           i_df_or_b = i_df_or_b + 1

        # Plotting
        for i_Λ, Λ in enumerate(statistics_names): # A figure by statistic
           plt.figure(figsize=(12, 7), dpi=80, facecolor='w')
           ax = plt.gca()
           for i_df_or_b, df_or_b in enumerate(df_or_b_vec): # Plotting the different histograms for each value of df_or_b_vec
                # Sometimes infinite value can appear due to overflow, we skip them
                λ_to_plot = λ[:,i_Λ, i_df_or_b]
                λ_to_plot = λ_to_plot[λ_to_plot!=np.inf]
                if statistics_scale[i_Λ] == 'log':
                   plt.hist(np.log(λ_to_plot), number_of_bins_histogram,
                            label=r'$%s=%.2f$' % (name_param, df_or_b), alpha=0.5, density=True,
                            edgecolor='black', linewidth=0.5)
                else:
                   plt.hist(λ_to_plot, number_of_bins_histogram, 
                            label=r'$%s=%.2f$' % (name_param, df_or_b), alpha=0.5, density=True,
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


def print_bunny():
    print( '|￣￣￣￣￣￣￣￣|')
    print( '|   COMPUTING   |') 
    print( '|   in progress |')
    print( '|               |' )  
    print( '| ＿＿＿_＿＿＿＿|') 
    print( ' (\__/) ||') 
    print( ' (•ㅅ•) || ')
    print( ' / 　 づ')

if __name__ == '__main__':

    #########################################################################################
    # Simulation parameters
    #########################################################################################

    # Enable parallel processing (or not)
    # In general the optimal parameters are obtained for 
    # number_of_threads = number of cores on the machine
    enable_multi = True
    number_of_threads = 6 # for parallel compuatation


    # General parameters
    p = 3
    N = 10
    T = 3

    # Monte-Carlo parameters
    number_of_trials = 6000
    if (number_of_trials%number_of_threads) !=0:
        number_of_trials = int(number_of_trials/number_of_threads) * number_of_threads

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
                        r'Wald statistic (Still buggy)',
                        r'Robust Shape statistic',
                        r'Robust Scale statistic',
                        r'Robust Scale and/or Shape statistic']

    # If the scale of the statistic is to be plotted in logarithmic or linear
    statistics_scale = ['log', 'linear', 'linear', 'log', 'log', 'log']
    
    # Plotting options
    number_of_bins_histogram = 50
    
    #########################################################################################
    # Printing some information about simulation
    #########################################################################################
    
    print("This simulation is aimed at testing CFAR properties of various statistics.")
    print("Parameters of the simulation:")
    print("    * p=%d, N=%d, T=%d" % (p,N,T))
    print("    * %d Monte-Carlo Trials will be done" % number_of_trials)
    if enable_multi:
        print("    * The simulation will be done using %d threads in parallel" % number_of_threads)


    #########################################################################################
    # Menu to select a type of simulation
    #########################################################################################
    # Letting the choice of the statistics to use
    print("The following statistics are available (Name: Arguments):")
    for i_Λ, Λ in enumerate(statistics_names):
        print("    %d. %s: %s" % (i_Λ, Λ, str(statistics_args[i_Λ])))
    print("Select which to use (as a list of format: \"[0,1]\"):")
    list_chosen = ast.literal_eval(input())
    statistics_list = [statistics_list[x] for x in list_chosen]
    statistics_args = [statistics_args[x] for x in list_chosen]
    statistics_names = [statistics_names[x] for x in list_chosen]
    statistics_scale = [statistics_scale[x] for x in list_chosen]
    # Setting a function in order to compute all the chosen statistics for each trial
    function_to_compute = compute_several_statistics
    function_args = [statistics_list, statistics_args]

    # Letting the choice of the type of distribution
    possibilities_distribution = {0 : 'Gaussian distribution',
       1 : 'K-distribution',
       2 : 't-distribution',
       3 : 'Cauchy distribution (Still buggy)',
       4 : 'Laplace distribution'}
    choice_distrib = None
    while choice_distrib not in possibilities_distribution.keys():
        print("Select the distribution among:")
        for i_dis in possibilities_distribution.keys():
            print("%d. %s" %(i_dis, possibilities_distribution[i_dis]))
        choice_distrib = int(input())


    # Letting the choice to test Covariance, texture with non-equality and texture with equality
    choice = None
    possibilities_cfar = {0 : covariance_cfar,
       1 : texture_non_equality_cfar,
       2 : texture_equality_cfar}
    while choice not in possibilities_cfar.keys():
        print("Please chose an option:")
        print("0. Covariance CFAR")
        print("1. Texture CFAR with no equality between dates")
        print("2. Texture CFAR with equality between dates")
        choice = int(input())

    #########################################################################################
    # Doing simulation
    #########################################################################################
    possibilities_cfar[choice]()

