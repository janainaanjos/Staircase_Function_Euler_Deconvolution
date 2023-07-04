# Variable regularization degrees in processing aeromagnetic data with first-order derivatives to improve geological mapping and automatic depth estimates

by Janaína A. Melo (IAG-USP), Carlos A. Mendonça (IAG-USP) and Yara R. Marangoni (IAG-USP) (2023)

Email: janaina.melo@usp.br (J.A. Melo); carlos.mendonca@iag.usp.br (C.A. Mendonça); yaramaran@usp.br. (Y.R. Marangoni)

This repository contains several Python codes and synthetic data examples for a paper about the open-source package staircase function and Euler deconvolution. The synthetic example data reproduce the results and figures shown in the publication. The codes 'filtering.py', 'synthetic_data.py', 'euler.py', and 'plot_figure.py' are compatible with both Python 2.7 and Python 3.7 programming languages. To run the main program "synthetic_data.py" and "euler.py", the numpy library is required. In addition, the numpy and sklearn.linear_model libraries are necessary to run the complementary code "filtering.py" and the matplotlib and numpy libraries are required to run the script "plot_figure.py".


1 - Abstract

2 - Content

3 - Prerequisites

4 - Parameterization

5 - Running the files

6 - License

===============================================================================


1 - Abstract
----------------------
The evaluation of numerical derivatives serves as the basis of several techniques applied to magnetic data processing to map structural trends and determine depth estimates for the associated magnetic sources. This differential operation is unstable, since amplifying the high-frequency content of the observed anomaly. Regularized derivatives based on Tikhonov regularization have been proposed to balance the oscillatory property of the numerical derivatives at the cost of incorporating a certain smoothing degree as determined by a particular choice for the associated regularization parameter. We apply a graphical procedure to select the regularization parameter for different transformations requiring first-order derivatives, by normalizing the L2-norm of transformed fields for trial regularization parameters. This approach produces a characteristic staircase function ranging from 0 to 1, respectively from subtle to over-regularized conditions, the ramp between the two end-term steps better trading-off noise amplification with smoothed outputs. As illustrated with synthetic and real data applications, a proper choice for the regularization parameter at the ramp interval can be made according to general properties of the transformed outputs (with oscillations or not) or inconsistency for depth to the top estimates, for example when inferring magnetic sources above the ground surface or at stratigraphic levels incompatible with a priori geological information. For grid operations, a mild regularization dose at the staircase ramp seems enough to produce reliable products, better centring fields over sources meanwhile preventing noise amplification. Real data application for the intracontinental Transbrasiliano shear-zone corridor, revealed that moderate regularization is necessary to obtain depth estimates consistent with independent geological and geophysical information.  


2 - Content
----------------------

Python codes:

	- filtering.py:
		General Python module containing the functions to compute the non-regularized and 
        	regularized directional first-order derivatives, S-function of the regularized 
		directional first-order derivatives, regularization parameters, analytical signal 
		amplitude and tilt derivative.
	
	- synthetic_data.py:
		Python script to generate the synthetic results. The script loads the total-field 
		anomaly of a synthetic model without and with noise of 0.1% e 1% from the files 
		"nonoise_synthetic_data.dat", "noise01_synthetic_data.dat", and "noise1_synthetic_
		data.dat", respectively. The script computes the S-function of the regularized 
		directional first-order derivatives, regularization parameters, non-regularized 
		and regularized directional first-order derivatives, analytical signal amplitude, 
		and tilt derivative using the functions in "filtering.py". The figures are generated 
		using the function "plot_figure.py".

	- euler.py: 
		Python script of Melo & Barbosa (2020), avaliable at https://github.com/ffigura/Euler
		-deconvolution-python, was modified to perform the regularized Euler deconvolution on 
		gridded data. In the script of Melo & Barbosa (2020), we only incorporated the function 
		"regularized_deriv" which calculates the regularized derivatives, and the function 
		"euler_deconv_regularized" which computes the Euler deconvolution using regularized 
		derivatives. 	
	
	- plot_figure.py:
		Python script to generate the figures of the synthetic data.
	
Outputs (folders): 
 
	- figures - figures of the synthetic data will be saved in this folder. The figures 
		    generated may differ from the publication and can be adapted in the script 
		    'plot_figure.py'.
						 
	- results - txt files of the solutions [x, y, depth, base level] of the regularized and
		    non-regularized Euler deconvolution for corrupted synthetic data with 1% noise.

Test data:

	- nonoise_synthetic_data.dat:
		Synthetic total-field anomaly data without noise

	- noise01_synthetic_data.dat:
		Synthetic total-field anomaly data corrupted with 0.1% noise

	- noise1_synthetic_data.dat:
		Synthetic total-field anomaly data corrupted with 1% noise	


3 - Prerequisites
----------------------
The Python program "filtering.py" requires the Python packages "numpy" and "sklearn", the 
scripts "synthetic_data.py" and "euler.py" require the Python package "numpy", and the script 
"plot_figure.py" requires the Python packages "numpy" and "matplotlib". 
The easier way to get Python and all libraries installed is through the Anaconda Python 
distribution (https://www.anaconda.com/distribution/). After installed Anaconda, install the libraries 
by running the following command in your terminal:

	conda install numpy matplotlib sklearn

 
4 - Parameterization
----------------------
To construct the S-function, the user must configure the parameter 'alpha_test' that corresponds 
to the variation range of the trial regularization parameters to calculate the Euclidean norm of 
the regularized directional derivatives.

To determine the regularization parameter associated with a particular Euclidean norm value, 
'value_norm', located in the slope of the S-function, the user must configure two parameters: 
'inferior_limit' and 'upper_limit'. These parameters represent the interval limits in which the 
S-function presents a linear variation.

To run the Euler deconvolution, the user must define three parameters: 1) "winsize" - moving data 
window size; 2) "filt" - percentage of the solutions with the higher vertical derivatives that 
will be kept; 3) "SI" -structural index.


5 - Running the files
----------------------
After Anaconda and the libraries are installed, you need to run 
the code 'synthetic_data.py' to generate the synthetic data results.
There are two easy ways to run:

1. Via command line: 

	1.1 On linux based systems open the command prompt windows (terminal), and run the command:
	
	On Python 2.7: python synthetic_data.py 

	On Python 3.7: python3 synthetic_data.py 

	1.2 On windows based systems open the Iphyton console (Python 2.7) or 
	Anaconda Navigator and then Jupyter QtConsole (Python 3.7), and run the command:

	run synthetic_data.py

2. Open the PyCharm program that comes with the Ananconda installation. In the PyCharm, open
the example script "synthetic_data.py" and click in "Run" inside the PyCharm
to execute it.

If everything runs ok, the figures and results will be saved, respectively, in the folders: 
'/figures' and '/results'.


6 - License
----------------------
The following legal note is restricted solely to the content of the named files. It cannot
overrule licenses from the Python standard distribution modules, which are imported and
used therein.

The 'filtering.py', 'synthetic_data.py', 'euler.py', and 'plot_figure.py' files are 
distributed under the following license agreement:

MIT license

Copyright (c) 2023 Janaína A. Melo, Carlos A. Mendonça and Yara R. Marangoni.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the names of the copyright holders nor the names of any contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
