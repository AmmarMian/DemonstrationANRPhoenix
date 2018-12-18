##############################################################################
# Some test on Snjoaq data UAVSAR
# WARNING: will download about 28 Go of data
# Authored by Ammar Mian, 12/11/2018
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
from generic_functions import *
import matplotlib.pyplot as plt
from monte_carlo_tools import *
from multivariate_images_tools import *
from change_detection_functions import *
from read_sar_data import *
import os
import time
import seaborn as sns
sns.set_style("darkgrid")


if __name__ == '__main__':

    # Activate latex in figures (or not)
    latex_in_figures = True
    if latex_in_figures:
      enable_latex_infigures()

    # Enable parallel processing (or not)
    enable_multi = True
    # These two variables serves to split the original image into sub-images to be treated in parallel
    # In general the optimal parameters are obtained for 
    # number_of_threads_rows*number_of_threads_columns = number of cores on the machine
    number_of_threads_rows = 8
    number_of_threads_columns = 6

    
    # Reading data using the class
    print( '|￣￣￣￣￣￣￣￣|')
    print( '|   READING     |') 
    print( '|   dataset     |')
    print( '|               |' )  
    print( '| ＿＿＿_＿＿＿＿|') 
    print( ' (\__/) ||') 
    print( ' (•ㅅ•) || ')
    print( ' / 　 づ')
    # Reading data using the class
    data_class = uavsar_slc_stack_1x1('D:/UAVSAR/Snjoaq/')
    data_class.read_data(polarisation=['HH', 'HV', 'VV'], segment=2, crop_indexes=[27000,30000,3000,4000])
   


    # Parameters
    image = data_class.data
    n_r, n_rc, p, T = image.shape
    windows_mask = np.ones((11,11))
    m_r, m_c = windows_mask.shape
    function_to_compute = compute_several_statistics
#    function_args = [[covariance_equality_glrt_gaussian_statistic,
#                        covariance_equality_t1_gaussian_statistic,
#                        scale_and_shape_equality_robust_statistic,
#                        shape_equality_robust_statistic,
#                        scale_equality_robust_statistic], 
#                     ["log", None] + [(0.001,30,"log")]*3]
#    if latex_in_figures:
#        statistic_names = ['$\hat{\Lambda}_{\mathrm{G}}$', '$\hat{\Lambda}_{t_1}$', 
#        '$\hat{\Lambda}_{\mathrm{MT}}$', '$\hat{\Lambda}_{\mathrm{Mat}}$', '$\hat{\Lambda}_{\mathrm{Tex}}$']
#    else:
#        statistic_names = ['Gaussian GLRT', 'Gaussian t1', 'Robust Scale and Shape', 'Robust Shape', 'Robust Scale']
        
        
    function_args = [[covariance_equality_glrt_gaussian_statistic,
                      scale_and_shape_equality_robust_statistic,
                      shape_equality_robust_statistic], 
                     ["log"]+ [(0.001,30,"log")]*2]
    if latex_in_figures:
        statistic_names = ['$\hat{\Lambda}_{\mathrm{G}}$']
    else:
        statistic_names = ['Gaussian GLRT']


    # Computing statistics on both images
    print( '|￣￣￣￣￣￣￣￣|')
    print( '|   COMPUTING   |') 
    print( '|   in progress |')
    print( '|               |' )  
    print( '| ＿＿＿_＿＿＿＿|') 
    print( ' (\__/) ||') 
    print( ' (•ㅅ•) || ')
    print( ' / 　 づ')
#    t_beginning = time.time()
#    results = sliding_windows_treatment_image_time_series_parallel(image, windows_mask, function_to_compute, 
#                    function_args, multi=enable_multi, number_of_threads_rows=number_of_threads_rows,
#                    number_of_threads_columns=number_of_threads_columns)
#    print("Elpased time: %d s" %(time.time()-t_beginning))


    # Plotting the time series in Grayscale representation
    dynamic = 50
    for t, image in enumerate(data_class.unique_identifiers_time_list):
        plt.figure(figsize=(9, 8), dpi=80, facecolor='w')
        Span = np.sum(np.sqrt(np.abs(data_class.data[:,:,:,t])**2), axis=2)
        plt.imshow(20*np.log10(Span), aspect='auto', cmap='gray',
                   vmin=20*np.log10(Span).max()-dynamic, vmax=20*np.log10(Span).max())
        plt.axis('off')
        plt.colorbar()
        plt.title('Image at time %d' % t)
    Span = None

    # Showing statistics results raw
#    for i_s, statistic in enumerate(statistic_names):
#        plt.figure(figsize=(9, 8), dpi=80, facecolor='w')
#        plt.imshow(results[:,:,i_s], aspect='auto', cmap='inferno',
#                   interpolation='spline36')
#        plt.title(statistic)
#        plt.colorbar()
#        plt.grid(False)
    plt.show()