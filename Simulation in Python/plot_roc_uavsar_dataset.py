##############################################################################
# ROC plots on real data UAVSAR
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
import wget
import os
import time
import seaborn as sns
sns.set_style("darkgrid")

def download_uavsar_cd_dataset(parth='./Data/'):

    # if directory exists just catch error
    try:
        os.mkdir(path)
    except:
        path

    # Links to UAVSAR datasets
    list_of_files = ['http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_09014_007_090423_L090HH_03_BC_s4_1x1.slc',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_09014_007_090423_L090HV_03_BC_s4_1x1.slc',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_09014_007_090423_L090VV_03_BC_s4_1x1.slc',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_09014_007_090423_L090HH_03_BC.ann',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_09014_007_090423_L090HV_03_BC.ann',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_09014_007_090423_L090VV_03_BC.ann',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_15059_006_150511_L090HH_03_BC_s4_1x1.slc',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_15059_006_150511_L090HV_03_BC_s4_1x1.slc',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_15059_006_150511_L090VV_03_BC_s4_1x1.slc',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_15059_006_150511_L090HH_03_BC.ann',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_15059_006_150511_L090HV_03_BC.ann',
                    'http://downloaduav2.jpl.nasa.gov/Release25/SanAnd_26524_03/SanAnd_26524_15059_006_150511_L090VV_03_BC.ann']

    for file_url in list_of_files:
        if not os.path.exists(path + file_url.split('/')[-1]):
            print("File %s not found, downloading it" % file_url.split('/')[-1])
            wget.download(url=file_url, out=path+file_url.split('/')[-1])


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
    number_of_threads_rows = 6
    number_of_threads_columns = 4

    # Downloading data if needed
    download_uavsar_cd_dataset()
    
    # Reading data using the class
    print( '|￣￣￣￣￣￣￣￣|')
    print( '|   READING     |') 
    print( '|   dataset     |')
    print( '|               |' )  
    print( '| ＿＿＿_＿＿＿＿|') 
    print( ' (\__/) ||') 
    print( ' (•ㅅ•) || ')
    print( ' / 　 づ')
    data_class = uavsar_slc_stack_1x1('./Data/')
    data_class.read_data(polarisation=['HH', 'HV', 'VV'], segment=4, crop_indexes=[25600,27900,3235,3835])
    image = data_class.data

    # Parameters
    n_r, n_rc, p, T = image.shape
    windows_mask = np.ones((11,11))
    m_r, m_c = windows_mask.shape
    function_to_compute = compute_several_statistics
    function_args = [[covariance_equality_glrt_gaussian_statistic,
                        covariance_equality_t1_gaussian_statistic,
                        scale_and_shape_equality_robust_statistic,
                        shape_equality_robust_statistic,
                        scale_equality_robust_statistic], 
                     [None, None] + [(0.001,30)]*3]
    if latex_in_figures:
        statistic_names = ['$\hat{\Lambda}_{\mathrm{G}}$', '$\hat{\Lambda}_{t_1}$', 
        '$\hat{\Lambda}_{\mathrm{MT}}$', '$\hat{\Lambda}_{\mathrm{Mat}}$', '$\hat{\Lambda}_{\mathrm{Tex}}$']
    else:
        statistic_names = ['Gaussian GLRT', 'Gaussian t1', 'Robust Scale and Shape', 'Robust Shape', 'Robust Scale']


    # Computing statistics on both images
    print( '|￣￣￣￣￣￣￣￣|')
    print( '|   COMPUTING   |') 
    print( '|   in progress |')
    print( '|               |' )  
    print( '| ＿＿＿_＿＿＿＿|') 
    print( ' (\__/) ||') 
    print( ' (•ㅅ•) || ')
    print( ' / 　 づ')
    t_beginning = time.time()
    results = sliding_windows_treatment_image_time_series_parallel(image, windows_mask, function_to_compute, 
                    function_args, multi=enable_multi, number_of_threads_rows=number_of_threads_rows,
                    number_of_threads_columns=number_of_threads_columns)
    print("Elpased time: %d s" %(time.time()-t_beginning))


    # Computing ROC curves
    number_of_points = 100
    ground_truth = np.load('./Data/ground_truth_uavsar_scene2.npy')
    ground_truth = ground_truth[int(m_r/2):-int(m_r/2), int(m_c/2):-int(m_c/2)]
    pfa_array = np.zeros((number_of_points, len(function_args[0])))
    pd_array = np.zeros((number_of_points, len(function_args[0])))
    for i_s, statistic in enumerate(statistic_names):

        # Sorting values of statistic
        λ_vec = np.sort(vec(results[:,:,i_s]), axis=0)
        λ_vec = λ_vec[np.logical_not(np.isinf(λ_vec))]

        # Selectionning number_of_points values from beginning to end
        indices_λ = np.floor(np.logspace(0, np.log10(len(λ_vec)-1), num=number_of_points))
        λ_vec = np.flip(λ_vec, axis=0)
        λ_vec = λ_vec[indices_λ.astype(int)]

        # Thresholding and summing for each value
        for i_λ, λ in enumerate(λ_vec):
            good_detection = (results[:,:,i_s] >= λ) * ground_truth
            false_alarms = (results[:,:,i_s] >= λ) * np.logical_not(ground_truth)
            pd_array[i_λ, i_s] = good_detection.sum() / (ground_truth==1).sum()
            pfa_array[i_λ, i_s] = false_alarms.sum() / (ground_truth==0).sum()

    # Showing images 
    for t in range(T):
        Span = np.sum(np.abs(image[:,:,:,t])**2, axis=2)/p
        plt.figure(figsize=(10, 12), dpi=80, facecolor='w')
        plt.imshow(20*np.log10(Span), aspect='auto', cmap='gray')
        plt.title('Image %d' %(t+1))
        plt.grid(False)

    # Showing ground truth
    plt.figure(figsize=(10, 12), dpi=80, facecolor='w')
    plt.imshow(ground_truth, aspect='auto')
    plt.title('Ground Truth')
    plt.grid(False)

    # Showing statistics results raw
    for i_s, statistic in enumerate(statistic_names):
        plt.figure(figsize=(10, 12), dpi=80, facecolor='w')
        plt.imshow(np.log10(results[:,:,i_s]), aspect='auto', cmap='jet',
                   interpolation='spline36')
        plt.title(statistic)
        plt.colorbar()
        plt.grid(False)

    # Showing statistics results ROC
    markers = ['o', 's', 'd', '*', '+']
    plt.figure(figsize=(16, 9), dpi=80, facecolor='w')
    for i_s, statistic in enumerate(statistic_names):
        plt.semilogx(pfa_array[:,i_s], pd_array[:,i_s], linestyle='--', label=statistic,
            marker=markers[i_s])
    plt.title('ROC curve scene 2')
    plt.legend()
    plt.show()
