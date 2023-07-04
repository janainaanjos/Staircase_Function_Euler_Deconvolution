"""
Synthetic data test 

Python script to generate the synthetic results. The script loads the total-field anomaly of a synthetic model without and with noise 
of 0.1% e 1% from the files "nonoise_synthetic_data.dat", "noise01_synthetic_data.dat", and "noise1_synthetic_data.dat", respectively.
The script computes the S-function of the regularized directional first-order derivatives, regularization parameters, non-regularized 
and regularized directional first-order derivatives, analytical signal amplitude, and tilt derivative using the functions in "filtering.py". 
The Euler deconvolution is performed by functions in "euler.py" of the program of Melo & Barbosa (2020). In the script of Melo & Barbosa 
(2020), we only incorporated the function "regularized_deriv" which calculates the regularized derivatives, and the function 
"euler_deconv_regularized" which computes the Euler deconvolution using regularized derivatives. The figures are generated using the functions 
in "plot_figure.py".

This code is released from the paper: "Variable regularization degrees in processing aeromagnetic data with first-order derivatives to 
improve geological mapping and automatic depth estimates".

The program is under the conditions terms in the file README.txt.

authors:Janaína A. Melo (IAG-USP), Carlos A. Mendonça (IAG-USP) and Yara R. Marangoni (IAG-USP) (2023)
email: janaina.melo@usp.br (J.A. Melo); carlos.mendonca@iag.usp.br (C.A. Mendonça); yaramaran@usp.br. (Y.R. Marangoni)
"""


"""
Input:

- input/nonoise_synthetic_data.dat: 2D-array with "n" rows by 4 columns, where "n" rows correspond to the size of the data.
            x-coordinate, y-coordinate, z-coordinate, total-field anomaly without noise

- input/noise01_synthetic_data.dat: 2D-array with "n" rows by 4 columns, where "n" rows correspond to the size of the data.
            x-coordinate, y-coordinate, z-coordinate, total-field anomaly corrupted with 0.1% noise

- input/noise1_synthetic_data.dat: 2D-array with "n" rows by 4 columns, where "n" rows correspond to the size of the data.
            x-coordinate, y-coordinate, z-coordinate, total-field anomaly corrupted with 1% noise
            

Parameters:

- Variation range of the trial regularization parameters: 
            alpha_test - 1D-array 

- Reference value of the S-function located on the slope (linear variation) of the S-function:
            value_norm - float
            
- Linear variation interval limits of the S-function:
            inferior_limit - float
            upper_limit - float

- Moving data window size:
            winsize - float

- Percentage of the solutions with higher vertical derivatives that will be kept:
            filt - float

- Structural index (SI):
            SI - integer

"""


import numpy as np
from filtering import *
from euler import *
from plot_figure import *




# Input data
true_data = np.loadtxt("input\nonoise_synthetic_data.dat")
data2 = np.loadtxt("input\noise01_synthetic_data.dat")
data = np.loadtxt("input\noise1_synthetic_data.dat")

x = data[:,0]                         # x-coordinates (m)
y = data[:,1]                         # y-coordinates (m)
z = data[:,2]                         # z-coordinates (m)
true_tfa = true_data[:,3]             # total-field anomaly without noise (nT)
tfa2 = data2[:,3]                     # total-field anomaly with noise of 0.1% (nT)
tfa = data[:,3]                       # total-field anomaly with noise of 1% (nT)

inc, dec = 45, -5                     # geomagnetic inclination and declination (degrees)

area = (0, 20000, 0, 20000)           # (x1, x2, y1, y2) - mesh boundaries
nx, ny = 200, 200                     # number of points in the x and y axis
shape = (nx, ny)                      # number of points in the x and y axis




'''
STAIRCASE FUNCTION OF THE DIRECTIONAL FIRST-ORDER DERIVATIVES
'''

# The user establishes the interval of the trial regularization parameters
l = np.arange(-6,14.5,0.5)
alpha_test = 10**(l[:])

# Calculates the Euclidean norm of the regularized directional derivatives of the total-field anomaly corrupted with 1% of noise 
norm_sol_dx, norm_sol_dy, norm_sol_dz = s_function_derivative(x, y, tfa, shape, alpha_test)


'''The user establishes the interval limits in which the S-function presents a linear variation to determine the regularization parameter 
associated with the Euclidean norm value equal to 0.50'''
value_norm = 0.50
upper_limit = 0.60
inferior_limit = 0.40

# Determines the regularization parameters of the directional derivatives to S=0.50
alpha_x_05 = regularization_parameter(norm_sol_dx, alpha_test, upper_limit, inferior_limit, value_norm)
alpha_y_05 = regularization_parameter(norm_sol_dy, alpha_test, upper_limit, inferior_limit, value_norm)
alpha_z_05 = regularization_parameter(norm_sol_dz, alpha_test, upper_limit, inferior_limit, value_norm)

alpha_vector = [alpha_x_05, alpha_y_05, alpha_z_05]

# Grid regularization parameter
alpha_grid = np.mean(alpha_vector)


'''The user establishes the interval limits in which the S-function presents a linear variation to determine the regularization parameter 
associated with the Euclidean norm value equal to 0.75'''
value_norm2 = 0.75
upper_limit2 = 0.80
inferior_limit2 = 0.60

# Determines the regularization parameters of the directional derivatives to S=0.75
alpha_x_075 = regularization_parameter(norm_sol_dx, alpha_test, upper_limit2, inferior_limit2, value_norm2)
alpha_y_075 = regularization_parameter(norm_sol_dy, alpha_test, upper_limit2, inferior_limit2, value_norm2)
alpha_z_075 = regularization_parameter(norm_sol_dz, alpha_test, upper_limit2, inferior_limit2, value_norm2)

alpha_vector2 = [alpha_x_075, alpha_y_075, alpha_z_075]

# Euler deconvolution regularization parameter
alpha_euler075 = np.mean(alpha_vector2)


'''The user establishes the interval limits in which the S-function presents a linear variation to determine the regularization parameter 
associated with the Euclidean norm value equal to 0.83'''
value_norm4 = 0.83
upper_limit4 = 0.90
inferior_limit4 = 0.70

# Determines the regularization parameters of the directional derivatives
alpha_x_083 = regularization_parameter(norm_sol_dx, alpha_test, upper_limit4, inferior_limit4, value_norm4)
alpha_y_083 = regularization_parameter(norm_sol_dy, alpha_test, upper_limit4, inferior_limit4, value_norm4)
alpha_z_083 = regularization_parameter(norm_sol_dz, alpha_test, upper_limit4, inferior_limit4, value_norm4)

alpha_vector4 = [alpha_x_083, alpha_y_083, alpha_z_083]
alpha_euler083 = np.mean(alpha_vector4)


'''The user establishes the interval limits in which the S-function presents a linear variation to determine the regularization parameter 
associated with the Euclidean norm value equal to 0.90'''
value_norm3 = 0.90
upper_limit3 = 0.95
inferior_limit3 = 0.80

# Determines the regularization parameters of the directional derivatives to S=0.90
alpha_x_090 = regularization_parameter(norm_sol_dx, alpha_test, upper_limit3, inferior_limit3, value_norm3)
alpha_y_090 = regularization_parameter(norm_sol_dy, alpha_test, upper_limit3, inferior_limit3, value_norm3)
alpha_z_090 = regularization_parameter(norm_sol_dz, alpha_test, upper_limit3, inferior_limit3, value_norm3)

alpha_vector3 = [alpha_x_090, alpha_y_090, alpha_z_090]

# Euler deconvolution regularization parameter
alpha_euler090 = np.mean(alpha_vector3)


# Print the exponents of the regularization parameters
print(np.round(alpha_vector[0], 1))
print(np.round(alpha_vector[1], 1))
print(np.round(alpha_vector[2], 1))
print(np.round(alpha_grid, 1))
print(np.round(alpha_euler075, 1))
print(np.round(alpha_euler083, 1))
print(np.round(alpha_euler090, 1))




'''
ANALYTIC SIGNAL AMPLITUDE (ASA) AND TILT DERIVATIVE (TDR)
'''

# First-order non-regularized derivatives (nT/m) of the total-field anomaly corrupted with 1% of noise 
dy_tfa, dx_tfa, dz_tfa = nonregularized_derivative(x, y, tfa, shape, order=1)
true_dy_tfa, true_dx_tfa, true_dz_tfa = nonregularized_derivative(x, y, true_tfa, shape, order=1)

# First-order regularized derivatives (nT/m) of the total-field anomaly corrupted with 1% of noise
reg_dy_tfa, reg_dx_tfa, reg_dz_tfa = regularized_derivative(x, y, tfa, shape, alpha=10**(alpha_grid))
reg2_dy_tfa, reg2_dx_tfa, reg2_dz_tfa = regularized_derivative(x, y, tfa, shape, alpha=10**(alpha_euler075))

# Non-regularized ASA (nT/m) and TDR (rad)
asa, tdr = asa_tdr(dx_tfa, dy_tfa, dz_tfa)
true_asa, true_tdr = asa_tdr(true_dx_tfa, true_dy_tfa, true_dz_tfa)

# Regularized ASA (nT/m) and TDR (rad)
reg_asa, reg_tdr = asa_tdr(reg_dx_tfa, reg_dy_tfa, reg_dz_tfa)
reg2_asa, reg2_tdr = asa_tdr(reg2_dx_tfa, reg2_dy_tfa, reg2_dz_tfa)




'''
EULER DECONVOLUTION - MELO & BARBOSA (2020)
'''

# The user sets the parameters

winsize = 6             # moving data window size
filt = 0.035            # percentage of 3.5% of the solutions with the higher vertical derivatives that will be kept
SI = 1                  # define the SI 

# Non-regularized Euler solutions [x, y, depth, base level]
euler_sol = euler_deconv(tfa, x, y, z, shape, area, SI, winsize, filt)
euler2_sol = euler_deconv(tfa2, x, y, z, shape, area, SI, winsize, filt)

# Regularized Euler solutions [x, y, depth, base level] to total-field anomaly corrupted with 1% of noise
reg_euler_sol = euler_deconv_regularized(tfa, x, y, z, shape, area, SI, winsize, filt, alpha=10**(alpha_euler083))
reg_euler_sol1 = euler_deconv_regularized(tfa, x, y, z, shape, area, SI, winsize, filt, alpha=10**(alpha_euler075))
reg_euler_sol2 = euler_deconv_regularized(tfa, x, y, z, shape, area, SI, winsize, filt, alpha=10**(alpha_euler090))

sol_depth = np.array([euler_sol[:,2], reg_euler_sol[:,2], reg_euler_sol1[:,2], reg_euler_sol2[:,2]])

# Regularized Euler solutions [x, y, depth, base level] to total-field anomaly corrupted with 0.1% of noise
reg_euler2_sol = euler_deconv_regularized(tfa2, x, y, z, shape, area, SI, winsize, filt, alpha=10**(alpha_euler075))
reg_euler2_sol1 = euler_deconv_regularized(tfa2, x, y, z, shape, area, SI, winsize, filt, alpha=10**(alpha_euler090))

sol_depth2 = np.array([euler2_sol[:,2], reg_euler2_sol[:,2], reg_euler2_sol1[:,2]])

# Saves estimates [x, y, depth, base level] in a txt file
np.savetxt('results\euler_solutions_synthetic.txt', euler_sol, delimiter="\t")
np.savetxt('results\reg_euler_solutions_synthetic.txt', reg_euler_sol, delimiter="\t")

# Depth ranges
xrange1 = []
xrange2 = []
xrange3 = []
yrange1 = []
yrange2 = []
yrange3 = []

x2range1 = []
x2range2 = []
x2range3 = []
y2range1 = []
y2range2 = []
y2range3 = []

reg_xrange1 = []
reg_xrange2 = []
reg_xrange3 = []
reg_yrange1 = []
reg_yrange2 = []
reg_yrange3 = []

reg_x2range1 = []
reg_x2range2 = []
reg_x2range3 = []
reg_y2range1 = []
reg_y2range2 = []
reg_y2range3 = []


for i in range (len(euler_sol[:,2])):
    if (euler_sol[i,2]>=195 and euler_sol[i,2]<=205):
        xrange1.append(euler_sol[i,0])
        yrange1.append(euler_sol[i,1])
    if (euler_sol[i,2]>=95 and euler_sol[i,2]<=105):
        xrange2.append(euler_sol[i,0])
        yrange2.append(euler_sol[i,1])
    else:
        xrange3.append(euler_sol[i,0])
        yrange3.append(euler_sol[i,1])


for i in range (len(reg_euler_sol[:,2])):
    if (reg_euler_sol[i,2]>=195 and reg_euler_sol[i,2]<=205):
        reg_xrange1.append(reg_euler_sol[i,0])
        reg_yrange1.append(reg_euler_sol[i,1])
    if (reg_euler_sol[i,2]>=95 and reg_euler_sol[i,2]<=105):
        reg_xrange2.append(reg_euler_sol[i,0])
        reg_yrange2.append(reg_euler_sol[i,1])
    else:
        reg_xrange3.append(reg_euler_sol[i,0])
        reg_yrange3.append(reg_euler_sol[i,1])


for i in range (len(euler2_sol[:,2])):
    if (euler2_sol[i,2]>=195 and euler2_sol[i,2]<=205):
        x2range1.append(euler2_sol[i,0])
        y2range1.append(euler2_sol[i,1])
    if (euler2_sol[i,2]>=95 and euler2_sol[i,2]<=105):
        x2range2.append(euler2_sol[i,0])
        y2range2.append(euler2_sol[i,1])
    else:
        x2range3.append(euler2_sol[i,0])
        y2range3.append(euler2_sol[i,1])


for i in range (len(reg_euler2_sol[:,2])):
    if (reg_euler2_sol[i,2]>=195 and reg_euler2_sol[i,2]<=205):
        reg_x2range1.append(reg_euler2_sol[i,0])
        reg_y2range1.append(reg_euler2_sol[i,1])
    if (reg_euler2_sol[i,2]>=95 and reg_euler2_sol[i,2]<=105):
        reg_x2range2.append(reg_euler2_sol[i,0])
        reg_y2range2.append(reg_euler2_sol[i,1])
    else:
        reg_x2range3.append(reg_euler2_sol[i,0])
        reg_y2range3.append(reg_euler2_sol[i,1])


xrange1 = np.array(xrange1)
xrange2 = np.array(xrange2)
xrange3 = np.array(xrange3)
yrange1 = np.array(yrange1)
yrange2 = np.array(yrange2)
yrange3 = np.array(yrange3)

x2range1 = np.array(x2range1)
x2range2 = np.array(x2range2)
x2range3 = np.array(x2range3)
y2range1 = np.array(y2range1)
y2range2 = np.array(y2range2)
y2range3 = np.array(y2range3)

reg_xrange1 = np.array(reg_xrange1)
reg_xrange2 = np.array(reg_xrange2)
reg_xrange3 = np.array(reg_xrange3)
reg_yrange1 = np.array(reg_yrange1)
reg_yrange2 = np.array(reg_yrange2)
reg_yrange3 = np.array(reg_yrange3)

reg_x2range1 = np.array(reg_x2range1)
reg_x2range2 = np.array(reg_x2range2)
reg_x2range3 = np.array(reg_x2range3)
reg_y2range1 = np.array(reg_y2range1)
reg_y2range2 = np.array(reg_y2range2)
reg_y2range3 = np.array(reg_y2range3)

xrange = [xrange1,xrange2,xrange3]
yrange = [yrange1,yrange2,yrange3]
reg_xrange = [reg_xrange1,reg_xrange2,reg_xrange3]
reg_yrange = [reg_yrange1,reg_yrange2,reg_yrange3]

x2range = [x2range1,x2range2,x2range3]
y2range = [y2range1,y2range2,y2range3]
reg_x2range = [reg_x2range1,reg_x2range2,reg_x2range3]
reg_y2range = [reg_y2range1,reg_y2range2,reg_y2range3]




'''
PLOT THE FIGURES
'''

# Plot the total-field anomaly, non-regularized and regularized ASA, and non-regularized and regularized TDR - Figure 1
plot_figure1(x, y, tfa, asa, reg_asa, tdr, reg_tdr, true_asa, reg2_asa)

# Plot regularization parameters to directional derivatives of the total-field anomaly - Figure 2
plot_figure2(alpha_test, norm_sol_dx, norm_sol_dy, norm_sol_dz, alpha_vector)

# Plot the histograms of the non-regularized and regularized Euler solutions to total-field anomaly corrupted with 0.1% of noise - Figure 3
plot_figure3(x, y, tfa2, x2range, y2range, reg_x2range, reg_y2range, sol_depth2)

# Plot the maps of the non-regularized and regularized Euler solutions total-field anomaly corrupted with 1% of noise - Figure 4
plot_figure4(x, y, tfa, xrange, yrange, reg_xrange, reg_yrange, sol_depth)

