##############################################################################
# Some data reading functions
# Authored by Ammar Mian, 09/11/2018
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
import os, fnmatch
import matplotlib.pyplot as plt 

class uavsar_slc_stack_1x1():
    """A class to store data corresponding to a SLC stack of UAVSAR data
        * path = path to a folder containing the files obtained from UAVSAR (.slc, .ann, .llh)"""
    def __init__(self, path):
        super(uavsar_slc_stack_1x1, self).__init__()
        self.path = path
        self.meta_data = {}
        self.llh_grid = {}
        self.slc_data = {}


    def read_data(self, polarisation=['HH', 'HV', 'VV'], segment=1, crop_indexes=None):
        """ A method to read UAVSAR SLC 1x1 data stack
            Inputs:
                * polarisation = a list of polarisations to read
                * crop_indexes = if we want to read a portion of the image, a list of the form
                    [lowerIndex axis 0, UpperIndex axis 0, lowerIndex axis 1, UpperIndex axis 1]"""

        # Obtain list of all files in directory
        listOfFiles = os.listdir(self.path)

        # Iterate on those files to search for an annotation file
        for entry in listOfFiles:  
            if fnmatch.fnmatch(entry, "*.ann"):

                # Read the ann file to obtain metadata
                self.meta_data[entry.split('.')[0]] = {} # initialise dict for the currect file
                with open(self.path + entry, 'r') as f:
                    for line in f: # Iterate on each line
                        # Discard commented lines
                        line = line.strip().split(';')[0]
                        if not (line == ''):
                            category = ' '.join(line.split()[:line.split().index('=')-1])
                            value = ' '.join(line.split()[line.split().index('=')+1:])
                            self.meta_data[entry.split('.')[0]][category] = value

        # Read slc file corresponding to the segment of interest and crop it

        # First, we obtain a list containing the different filenames for each date
        # We put POL at the place of the polarisation and SEGMENT at the place of the segment in order to replace it after
        self.unique_identifiers_time_list = []
        for entry in list(self.meta_data.keys()):
            unique_identifiers_time = '_'.join(entry.split('_')[:-2])[:-2] + "POL_" + \
                                      '_'.join(entry.split('_')[-2:]) + '_sSEGMENT'
            if unique_identifiers_time not in self.unique_identifiers_time_list:
                self.unique_identifiers_time_list.append(unique_identifiers_time)

        # Then we read the files one by one for each polarisation and time
        if crop_indexes is not None:
            self.data = np.zeros((crop_indexes[1]-crop_indexes[0], crop_indexes[3]-crop_indexes[2], len(polarisation), 
                        len(self.unique_identifiers_time_list)), dtype='complex64')
            for t, entry_time in enumerate(self.unique_identifiers_time_list):
                for i_pol, pol in enumerate(polarisation):
                    # Read slc file in a temporary array an then crop it
                    file_name = entry_time.replace('POL', pol).replace('SEGMENT', str(segment))
                    print("Reading %s" % (self.path+file_name))
                    shape = (int(self.meta_data['_'.join(file_name.split('_')[:-1])]['slc_1_1x1 Rows']),
                             int(self.meta_data['_'.join(file_name.split('_')[:-1])]['slc_1_1x1 Columns']))               
                    temp_array = np.fromfile( self.path + file_name + '_1x1.slc', dtype='complex64').reshape(shape)
                    temp_array = temp_array[crop_indexes[0]:crop_indexes[1], crop_indexes[2]:crop_indexes[3]]
                    self.data[:,:,i_pol,t] = temp_array
        else:
            for t, entry_time in enumerate(self.unique_identifiers_time_list):
                for i_pol, pol in enumerate(polarisation):
                    # Read slc file in a temporary array an then crop it
                    file_name = entry_time.replace('POL', pol).replace('SEGMENT', str(segment))
                    print("Reading %s" % (self.path+file_name))
                    shape = (int(self.meta_data['_'.join(file_name.split('_')[:-1])]['slc_1_1x1 Rows']),
                             int(self.meta_data['_'.join(file_name.split('_')[:-1])]['slc_1_1x1 Columns']))               
                    temp_array = np.fromfile( self.path + file_name + '_1x1.slc', dtype='complex64').reshape(shape)
                    if t==0 and i_pol==0:
                        self.data = temp_array
                    else:
                        self.data = np.dtack(self.data, temp_array)




if __name__ == '__main__':
    
    # Reading data using the class
    data_class = uavsar_slc_stack_1x1('D:/UAVSAR/SanAnd_26524_03/')
    data_class.read_data(polarisation=['HH', 'HV', 'VV'], segment=4, crop_indexes=[25600,27900,3235,3835])
   
    # Plotting the time series in Grayscale representation
    dynamic = 50
    for t, image in enumerate(data_class.unique_identifiers_time_list):
        plt.figure(figsize=(6, 7), dpi=80, facecolor='w')
        Span = np.sum(np.sqrt(np.abs(data_class.data[:,:,:,t])**2), axis=2)
        plt.imshow(20*np.log10(Span), aspect='auto', cmap='gray',
                   vmin=20*np.log10(Span).max()-dynamic, vmax=20*np.log10(Span).max())
        plt.axis('off')
        plt.colorbar()
        plt.title('Image at time %d' % t)
    Span = None

