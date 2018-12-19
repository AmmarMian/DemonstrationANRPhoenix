# Demonstration of Robust Statistics for Change Detection

This folder contains codes for a demonstration of the work I have done in the context of [ANR Phoenix](http://am.atto.free.fr/index_phoenix.htm).
The works consists in the development of robust statistics for Change Detection in SAR Images.

The description of the work can be found in:

	"New Robust Statistics for Change Detection in Time Series of Multivariate SAR Images",
	Ammar Mian , Guillaume Ginolhac , Jean-Philippe Ovarlez , Abdourahmane M. Atto,
	in Transactions on Signal Processing
	URL: https://ieeexplore.ieee.org/document/8552453
	Preprint available at: https://ammarmian.github.io/publication/tsp-2018/

If you use any code of the repository, please consider citing the above mentionned reference.

## Files' organisation

This folder is organised as follows:
 - **Simulation in Matlab/** contains code for the simulation in Matlab (2017a).
 
 	**WARNING**: This folder correspond to an older version of code not well commented. If you can, prefer the Python version which was specifically
 	developped to be shared.
 - **Simulation in Python/** contains code for simulation done in Python (3.7).


## Requirements for Python
	The code provided was developped and tested using Python 3.7. The following packages must be installed 
	in order to run everything smoothly:
	- Scipy/numpy
	- matplotlib
	- seaborn
	- wget (if you do not have the UAVSAR data yet)
	- tqdm

The code use parallel processing which can be disabled by putting the boolean enable_multi to False in each file.
The figures can be plotted using Latex formatting by putting the boolean latex_in_figures to True in each file (must have a latex distribution installed).

## Files' organisation in Simulation in Python folder

This folder is organised as follows:
 - **test_cfar_property.py** contains a code to test Matrix and CFAR properties of statistics.
 - **test_cfar_property_menu.py** is an interactive version of test_cfar_property.py, where a menu allows to chose the parameters of the simulation.
 - **plot_roc_uavsar_dataset.py** contains a code to compare the results of statistics on a real UAVSAR dataset.

 	 **WARNING**: If the data is not already in the path specified at the beginning of the file, it will download it automatically.
 	 		  The data is approximately 28 Go in size.
 - **generic_functions.py** contains some general use functions including random vectors generation ones.
 - **monte_carlo_tools.py** contains some functions used in order to compute Monte-Carlo simulations efficiently.
 - **multivariate_images_tools.py** contains some general use functions in order to compute statistics on a time series 
 	of image using a sliding windows.
 - **plot_roc_uavsar_dataset.py** contains a class used to read UAVSAR dataset.
 - **sar_time_series_functions.py.py** contains functions to generate time series of vectors according to the three models presented in subsection III.B of:

	"New Robust Statistics for Change Detection in Time Series of Multivariate SAR Images",
	Ammar Mian , Guillaume Ginolhac , Jean-Philippe Ovarlez , Abdourahmane M. Atto,
	in Transactions on Signal Processing
	URL: https://ieeexplore.ieee.org/document/8552453
	Preprint available at: https://ammarmian.github.io/publication/tsp-2018/

## Files' organisation in Simulation in Matlab folder
 - **ChangeDetection/** contains a code to compute change detection over real dataset UAVSAR. You must specify the good path to the data since there is no
 automatic download like in Python.
 - **Theoretical/** contains several codes for testing CFARness, plotting theoretical ROC, testing convergence properties.
 - **Detectors/** contains functions that computes the statistics for Change detection

## Credits
**Author:** Ammar Mian, Ph.d student at SONDRA, CentraleSupélec
 - **E-mail:** ammar.mian@centralesupelec.fr
 - **Web:** https://ammarmian.github.io/
 
 Acknowledgements to:
 - [**Guillaume Ginolhac**](https://www.listic.univ-smb.fr/presentation/membres/enseignants-chercheurs/guillaume-ginolhac/), LISTIC, Université Savoie Mont-Blanc
 - [**Jean-Philippe Ovarlez**](http://www.jeanphilippeovarlez.com/), DEMR, ONERA , Université Paris-Saclay  & SONDRA, CentraleSupélec
 - [**Abdourahmane M. Atto**](http://am.atto.free.fr/), LISTIC, Université Savoie Mont-Blanc

 
## Copyright
 
 Copyright 2018 @CentraleSupelec

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

