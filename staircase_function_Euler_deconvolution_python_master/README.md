# Variable regularization degrees in processing aeromagnetic data with first-order derivatives to improve geological mapping and automatic depth estimates

by
Janaína A. Melo, Carlos A. Mendonça and Yara R. Marangoni 

## About

This paper has been submitted in the *Geophysical Journal International*. Melo, J.A., Mendonça, C.A, Marangoni, Y.R., 2023. Variable regularization degrees in processing aeromagnetic data with first-order derivatives to improve geological mapping and automatic depth estimates.

This repository contains the source code to perform the synthetic test presented. The codes 'filtering.py', 'synthetic_data.py', 'euler.py', and 'plot_figure.py' generate the results related to our methodology.

These programs are compatible with both Python 2.7 and Python 3.7 programming language.


## Abstract

The evaluation of numerical derivatives serves as the basis of several techniques applied to magnetic data processing to map structural trends and determine depth estimates for the associated magnetic sources. This differential operation is unstable, since amplifying the high-frequency content of the observed anomaly. Regularized derivatives based on Tikhonov regularization have been proposed to balance the oscillatory property of the numerical derivatives at the cost of incorporating a certain smoothing degree as determined by a particular choice for the associated regularization parameter. We apply a graphical procedure to select the regularization parameter for different transformations requiring first-order derivatives, by normalizing the L2-norm of transformed fields for trial regularization parameters. This approach produces a characteristic staircase function ranging from 0 to 1, respectively from subtle to over-regularized conditions, the ramp between the two end-term steps better trading-off noise amplification with smoothed outputs. As illustrated with synthetic and real data applications, a proper choice for the regularization parameter at the ramp interval can be made according to general properties of the transformed outputs (with oscillations or not) or inconsistency for depth to the top estimates, for example when inferring magnetic sources above the ground surface or at stratigraphic levels incompatible with a priori geological information. For grid operations, a mild regularization dose at the staircase ramp seems enough to produce reliable products, better centring fields over sources meanwhile preventing noise amplification.  Real data application for the intracontinental Transbrasiliano shear-zone corridor, revealed that moderate regularization is necessary to obtain depth estimates consistent with independent geological and geophysical information.  


## Content

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
	
Test data:

- nonoise_synthetic_data.dat:
		Synthetic total-field anomaly data without noise

- noise01_synthetic_data.dat:
		Synthetic total-field anomaly data corrupted with 0.1% noise 

- noise1_synthetic_data.dat:
		Synthetic total-field anomaly data corrupted with 1% noise	


## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/janainaanjos/Staircase_Function_Euler_Deconvolution.git

or [download a zip archive](https://github.com/janainaanjos/Staircase_Function_Euler_Deconvolution/archive/master.zip).


## Dependencies

The Python program "filtering.py" requires the Python packages "numpy" and "sklearn", the 
scripts "synthetic_data.py" and "euler.py" require the Python package "numpy", and the script 
"plot_figure.py" requires the Python packages "numpy" and "matplotlib". 
The easier way to get Python and all libraries installed is through the Anaconda Python 
distribution (https://www.anaconda.com/distribution/). After installed Anaconda, install the libraries 
by running the following command in your terminal:

	conda install numpy matplotlib sklearn


## Reproducing the results

For the synthetic data, the results are reproducible from the folder '/results' and the figures 
are found in the folder '/figures'. Running the code 'synthetic_data.py' will allow the reproduction 
of the results of our methodology. For more information read the file 'README.md'or 
'README.txt' in the folder '/code'.


## License

All source code is made available under a MIT license. You can freely use 
and modify the code, without warranty, so long as you provide attribution
to the authors. See 'LICENSE.md' for the full license text.

The manuscript text is not open source. The authors reserve the rights to 
the article content, which is currently submitted for publication in the
*Geophysical Journal International*.
