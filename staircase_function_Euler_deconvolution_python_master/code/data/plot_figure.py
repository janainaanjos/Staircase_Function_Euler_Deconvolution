"""
Plot functions 

A Python program to plot the total-field anomaly, non-regularized and regularized ASA, non-regularized and regularized TDR, the S-function, depth
histograms, and map of the solutions Euler.

This code plot the figures 1, 2, 3, and 4 of the synthetic data in the folder 'figures'.

This code is released from the paper: "Variable regularization degrees in processing aeromagnetic data with first-order derivatives to 
improve geological mapping and automatic depth estimates".

The program is under the conditions terms in the file README.txt.

authors:Janaína A. Melo (IAG-USP), Carlos A. Mendonça (IAG-USP) and Yara R. Marangoni (IAG-USP) (2023)
email: janaina.melo@usp.br (J.A. Melo); carlos.mendonca@iag.usp.br (C.A. Mendonça); yaramaran@usp.br. (Y.R. Marangoni)
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



def plot_figure1(x, y, tfa, asa, reg_asa, tilt, reg_tilt, true_asa, reg2_asa):

    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(7, 10))

    true_asa = true_asa * 1000
    asa = asa * 1000
    reg_asa = reg_asa * 1000
    reg2_asa = reg2_asa * 1000

    xx = np.arange(5, 15, 0.1)

    ax[0][1].plot(xx, true_asa[30000:30100],'-o', color='black', markersize=2, label='sem ruído')
    ax[0][1].plot(xx, asa[30000:30100], color='lime', label='$S=1$')
    ax[0][1].plot(xx, reg2_asa[30000:30100], color='blueviolet', markersize=2, label='$S=0.75$')
    ax[0][1].plot(xx, reg_asa[30000:30100], color='red', markersize=2, label='$S=0.50$')

    v1_ = np.linspace(min(tfa), max(tfa), 5, endpoint=True)
    v2 = np.linspace(min(reg_asa), max(reg_asa), 15, endpoint=True)
    v2_ = np.linspace(min(reg_asa), max(reg_asa), 5, endpoint=True)
    v3 = np.linspace(min(reg_tilt), max(reg_tilt), 15, endpoint=True)
    v3_ = np.linspace(min(reg_tilt), max(reg_tilt), 5, endpoint=True)

    tmp1 = ax[0][0].tricontourf(y/1000, x/1000, tfa, 30, cmap='gist_ncar')
    tmp2 = ax[1][0].tricontourf(y/1000, x/1000, asa, 30, cmap='gist_ncar', levels=v2)
    tmp3 = ax[1][1].tricontourf(y/1000, x/1000, reg_asa, 30, cmap='gist_ncar', levels=v2)
    tmp4 = ax[2][0].tricontourf(y/1000, x/1000, tilt, 30, cmap='gist_ncar', levels=v3)
    tmp5 = ax[2][1].tricontourf(y/1000, x/1000, reg_tilt, 30, cmap='gist_ncar', levels=v3)


    plt.colorbar(tmp1, ax=ax[0][0], fraction=0.030, aspect=20, spacing='uniform', format='%.f', orientation='vertical',
                      ticks=v1_).set_label('(nT)', fontsize=10, labelpad=-20, y=-0.10, rotation=0)

    plt.colorbar(tmp2, ax=ax[1][0], fraction=0.030, aspect=20, spacing='uniform', format='%.f', orientation='vertical',
                      ticks=v2_).set_label('(nT/km)', fontsize=10, labelpad=-10, y=-0.10, rotation=0)

    plt.colorbar(tmp3, ax=ax[1][1], fraction=0.030, aspect=20, spacing='uniform', format='%.f', orientation='vertical',
                      ticks=v2_).set_label('(nT/km)', fontsize=10, labelpad=-10, y=-0.10, rotation=0)

    plt.colorbar(tmp4, ax=ax[2][0], fraction=0.030, aspect=20, spacing='uniform', orientation='vertical', format='%.1f',
                      ticks=v3_).set_label('(rad)', fontsize=10, labelpad=-20, y=-0.10, rotation=0)

    plt.colorbar(tmp5, ax=ax[2][1], fraction=0.030, aspect=20, spacing='uniform', orientation='vertical', format='%.1f',
                      ticks=v3_).set_label('(rad)', fontsize=10, labelpad=-20, y=-0.10, rotation=0)

    
    ax[1][0].vlines(11, 5, 15, colors = 'black', linestyles = 'dashed', lw=2)
    ax[1][0].text(12, 5, 'S', fontsize=10, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[1][0].text(12, 15, 'N', fontsize=10, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[0][1].text(4.5, 330, 'S', fontsize=10, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[0][1].text(15, 330, 'N', fontsize=10, fontweight='bold', horizontalalignment='center', verticalalignment='center')

    ax[0][0].set_xlabel('y (km)')
    ax[0][1].set_xlabel('x (km)')
    ax[1][0].set_xlabel('y (km)')
    ax[1][1].set_xlabel('y (km)')
    ax[2][0].set_xlabel('y (km)')
    ax[2][1].set_xlabel('y (km)')
    ax[0][0].set_ylabel('x (km)')
    ax[0][1].set_ylabel('ASA (nT/km)')
    ax[1][0].set_ylabel('x (km)')
    ax[2][0].set_ylabel('x (km)')
    ax[1][1].set_ylabel('x (km)')
    ax[2][1].set_ylabel('x (km)')

    ax[0][0].yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax[0][1].yaxis.set_ticks([0, 80, 160, 240, 320])
    ax[1][0].yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax[1][1].yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax[2][0].yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax[2][1].yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax[0][0].xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))
    ax[0][1].xaxis.set_ticks([5,7,9,11,13,15])
    ax[1][0].xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))
    ax[1][1].xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))
    ax[2][0].xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))
    ax[2][1].xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))

    ax[0][0].text(3, 16, 'A', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[1][0].text(3, 16, 'A', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[1][1].text(3, 16, 'A', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[2][0].text(3, 16, 'A', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[2][1].text(3, 16, 'A', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')

    ax[0][0].text(9, 11, 'B', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[1][0].text(9, 11, 'B', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[1][1].text(9, 11, 'B', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[2][0].text(9, 11, 'B', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax[2][1].text(9, 11, 'B', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')

    ax[0][0].text(-4, 22, 'a)', fontsize=14, horizontalalignment='center', verticalalignment='center')
    ax[0][1].text(2.5, 352, 'd)', fontsize=14, horizontalalignment='center', verticalalignment='center')
    ax[1][0].text(-4, 22, 'b)', fontsize=14, horizontalalignment='center', verticalalignment='center')
    ax[1][1].text(-4, 22, 'e)', fontsize=14, horizontalalignment='center', verticalalignment='center')
    ax[2][0].text(-4, 22, 'c)', fontsize=14, horizontalalignment='center', verticalalignment='center')
    ax[2][1].text(-4, 22, 'f)', fontsize=14, horizontalalignment='center', verticalalignment='center')

    plt.subplots_adjust(wspace=0.6, hspace=0.5)
    ax[0][1].legend(edgecolor='black',loc='upper right', fontsize=7)

    plt.savefig('figures\FIG1.png', bbox_inches='tight', dpi=600)
    plt.close('all')

    return



def plot_figure2(alpha, norm_sol_dx, norm_sol_dy, norm_sol_dz, alpha_vector):

    fig, ax = plt.subplots(figsize=(4, 3))

    ax.plot(np.log10(alpha), norm_sol_dx, '-', color='purple', linewidth=1.5, label='S$_{x}(\u03B1)$')
    ax.plot(np.log10(alpha), norm_sol_dy, color='red', linewidth=1.5, label='S$_{y}(\u03B1)$')
    ax.plot(np.log10(alpha), norm_sol_dz, color='limegreen', linewidth=1.5, label='S$_{z}(\u03B1)$')

    ax.set_ylabel('S$_{\mu}(\u03B1)$', color='black', fontsize=7)
    ax.set_xlabel('log$_{10}$($\u03B1$)', fontsize=7)
    ax.xaxis.set_ticks(np.arange(min(np.log10(alpha)), max(np.log10(alpha)), step=3))

    ax.text(5.3, 0.02, '$\u03B1_{x}=10^{4.3}$', fontsize=5, color='purple', horizontalalignment='center', verticalalignment='center')
    ax.text(5.2, -0.026, '$\u03B1_{y}=10^{4.2}$', fontsize=5, color='red', horizontalalignment='center', verticalalignment='center')
    ax.text(2.5, 0.02, '$\u03B1_{z}=10^{3.9}$', fontsize=5, color='limegreen', horizontalalignment='center', verticalalignment='center')
    ax.text(-4.1, 0.78, 'S$_{g}(10^{3.1})$=0.75', fontsize=5, color='black', horizontalalignment='center',
            verticalalignment='center')
    ax.text(-4.22, 0.53, 'S$_{g}(10^{4.1})$=0.5', fontsize=5, color='black', horizontalalignment='center',
            verticalalignment='center')
    ax.text(-4.1, 0.86, 'S$_{g}(10^{2.9})$=0.83', fontsize=5, color='black', horizontalalignment='center',
            verticalalignment='center')
    ax.text(-4.22, 0.93, 'S$_{g}(10^{2.5})$=0.9', fontsize=5, color='black', horizontalalignment='center',
            verticalalignment='center')

    ax.hlines(0.5, -6, 14, colors='black', linestyles='dashed', linewidth=0.75)
    ax.hlines(0.75, -6, 14, colors='black', linestyles='dashed', linewidth=0.75)
    ax.hlines(0.83, -6, 14, colors='black', linestyles='dashed', linewidth=0.75)
    ax.hlines(0.90, -6, 14, colors='black', linestyles='dashed', linewidth=0.75)

    ax.vlines(np.round(alpha_vector[0], 2), 0, 0.5, colors = 'purple', linestyles = 'dashed', linewidth=0.75)
    ax.vlines(np.round(alpha_vector[1], 1), 0, 0.5, colors = 'red', linestyles = 'dashed', linewidth=0.75)
    ax.vlines(np.round(alpha_vector[2], 1), 0, 0.5, colors = 'limegreen', linestyles = 'dashed', linewidth=0.75)

    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.legend(loc='lower left', fontsize=5, edgecolor='black')

    plt.savefig('figures\FIG2.png', bbox_inches='tight', dpi=600)
    plt.close('all')

    return



def plot_figure3(x, y, tfa2, x2range, y2range, reg_x2range, reg_y2range, sol_depth2):

    fig = plt.figure(figsize=(6, 10))
    ax1 = plt.subplot2grid((3, 2), (0, 0))
    ax2 = plt.subplot2grid((3, 2), (0, 1))
    ax3 = plt.subplot2grid((3, 2), (1, 0), colspan=2)
    ax4 = plt.subplot2grid((3, 2), (2, 0), colspan=2)

    ax1.tricontourf(y/1000, x/1000, tfa2, 5, cmap='gist_gray')
    ax2.tricontourf(y/1000, x/1000, tfa2, 5, cmap='gist_gray')

    ax1.scatter(y2range[2]/1000, x2range[2]/1000, s=6, marker='o', c='yellow')
    ax1.scatter(y2range[0]/1000, x2range[0]/1000, s=6, marker='o', c='blue')
    ax1.scatter(y2range[1]/1000, x2range[1]/1000, s=6, marker='o', c='red')
    ax2.scatter(reg_y2range[2]/1000, reg_x2range[2]/1000, s=6, marker='o', c='yellow')
    ax2.scatter(reg_y2range[0]/1000, reg_x2range[0]/1000, s=6, marker='o', c='blue')
    ax2.scatter(reg_y2range[1]/1000, reg_x2range[1]/1000, s=6, marker='o', c='red')

    N1, bins1, patches1 = ax3.hist(sol_depth2[0], density=True, color='lightgrey', edgecolor='black', bins=65, alpha=1, label='$S_{g}=1$')
    N1, bins1, patches1 = ax4.hist(sol_depth2[0], density=True, color='lightgrey', edgecolor='black', bins=65, alpha=1, label='$S_{g}=1$')
    N2, bins2, patches2 = ax3.hist(sol_depth2[2], density=True, color='lime', edgecolor='black', bins=65, alpha=1, label='$S_{g}=$0.9')
    N3, bins3, patches3 = ax4.hist(sol_depth2[1], density=True, color='lime', edgecolor='black', bins=65, alpha=1, label='$S_{g}=$0.75')

    line = [Line2D([0], [0], marker='o', color='w', label='100 $\pm$ 5', markerfacecolor='red', markersize=7,
                    markeredgecolor='black'),
             Line2D([0], [0], marker='o', color='w', label='200 $\pm$ 5', markerfacecolor='blue', markersize=7,
                    markeredgecolor='black'),
             Line2D([0], [0], marker='o', color='w', label='other', markerfacecolor='yellow', markersize=7,
                    markeredgecolor='black')]

    ax2.legend(handles=line, loc=4, edgecolor='black', title='Depth (m)', title_fontsize=8,  numpoints=1, ncol=1, fontsize=7)

    ax1.set_xlabel('y (km)', fontsize=10)
    ax2.set_xlabel('y (km)',fontsize=10)
    ax1.set_ylabel('x (km)', fontsize=10)
    ax2.set_ylabel('x (km)', fontsize=10)

    ax3.set_xlabel('Depth (m)', fontsize=10)
    ax4.set_xlabel('Depth (m)', fontsize=10)
    ax3.set_ylabel('Probability density', fontsize=10)
    ax4.set_ylabel('Probability density', fontsize=10)

    ax1.tick_params(axis='both',which='major', labelsize=9)
    ax2.tick_params(axis='both',which='major', labelsize=9)
    ax3.tick_params(axis='both',which='major', labelsize=9)
    ax4.tick_params(axis='both',which='major', labelsize=9)

    ax1.yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax2.yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax1.xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))
    ax2.xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))

    ax3.vlines(x=100, ymin=0,ymax= 0.03, color='red', linestyles='dashed')
    ax3.vlines(x=200, ymin=0,ymax= 0.03, color='red', linestyle='dashed')
    ax4.vlines(x=100, ymin=0, ymax=0.03, color='red', linestyles='dashed')
    ax4.vlines(x=200, ymin=0, ymax=0.03, color='red', linestyle='dashed')

    ax3.text(100, 0.033, '100 (A)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')
    ax3.text(200, 0.033, '200 (B)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')
    ax4.text(100, 0.033, '100 (A)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')
    ax4.text(200, 0.033, '200 (B)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')

    ax3.yaxis.set_ticks([0, 0.01, 0.02, 0.03, 0.04,0.05])
    ax4.yaxis.set_ticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])
   
    ax3.set_xlim([60, 240])
    ax4.set_xlim([60, 240])

    ax3.xaxis.set_ticks([60, 80, 100, 120, 140, 160, 180, 200, 220, 240])
    ax4.xaxis.set_ticks([60, 80, 100, 120, 140, 160, 180, 200, 220, 240])

    ax1.text(3, 16, 'A', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax2.text(3, 16, 'A', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax1.text(9, 11, 'B', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax2.text(9, 11, 'B', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')

    ax3.legend(edgecolor='black', fontsize=9, loc='upper left')
    ax4.legend(edgecolor='black', fontsize=9, loc='upper left')

    ax1.text(-4.2, 20, 'a)', fontsize=13, horizontalalignment='center', verticalalignment='center')
    ax2.text(-4.2, 20, 'b)', fontsize=13, horizontalalignment='center', verticalalignment='center')
    ax3.text(40, 0.053, 'c)', fontsize=13, horizontalalignment='center', verticalalignment='center')
    ax4.text(40, 0.053, 'd)', fontsize=13, horizontalalignment='center', verticalalignment='center')

    plt.subplots_adjust(hspace=0.36, wspace=0.32)

    plt.savefig('figures\FIG3.png', bbox_inches='tight', dpi=600)
    plt.close('all')

    return



def plot_figure4(x, y, tfa, xrange, yrange, reg_xrange, reg_yrange, sol_depth):

    fig = plt.figure(figsize=(6, 10))
    ax1 = plt.subplot2grid((3, 2), (0, 0))
    ax2 = plt.subplot2grid((3, 2), (0, 1))
    ax3 = plt.subplot2grid((3, 2), (1, 0), colspan=2)
    ax4 = plt.subplot2grid((3, 2), (2, 0), colspan=2)

    ax1.tricontourf(y/1000, x/1000, tfa, 5, cmap='gist_gray')
    ax2.tricontourf(y/1000, x/1000, tfa, 5, cmap='gist_gray')

    ax1.scatter(yrange[2]/1000, xrange[2]/1000, s=6, marker='o', c='yellow')
    ax1.scatter(yrange[0]/1000, xrange[0]/1000, s=6, marker='o', c='blue')
    ax1.scatter(yrange[1]/1000, xrange[1]/1000, s=6, marker='o', c='red')
    ax2.scatter(reg_yrange[2]/1000, reg_xrange[2]/1000, s=6, marker='o', c='yellow')
    ax2.scatter(reg_yrange[0]/1000, reg_xrange[0]/1000, s=6, marker='o', c='blue')
    ax2.scatter(reg_yrange[1]/1000, reg_xrange[1]/1000, s=6, marker='o', c='red')

    N1, bins1, patches1 = ax3.hist(sol_depth[0], density=True, color='lightgrey', edgecolor='black', bins=80, alpha=1, label='$S_{g}=1$')
    N1, bins1, patches1 = ax4.hist(sol_depth[0], density=True, color='lightgrey', edgecolor='black', bins=80, alpha=1, label='$S_{g}=1$')
    N2, bins2, patches2 = ax3.hist(sol_depth[3], density=True, color='lime', edgecolor='black', bins=80, alpha=1, label='$S_{g}=$0.9')
    N3, bins3, patches3 = ax4.hist(sol_depth[1], density=True, color='lime', edgecolor='black', bins=80, alpha=1, label='$S_{g}=$0.83')

    line = [Line2D([0], [0], marker='o', color='w', label='100 $\pm$ 5', markerfacecolor='red', markersize=7,
                    markeredgecolor='black'),
             Line2D([0], [0], marker='o', color='w', label='200 $\pm$ 5', markerfacecolor='blue', markersize=7,
                    markeredgecolor='black'),
             Line2D([0], [0], marker='o', color='w', label='others', markerfacecolor='yellow', markersize=7,
                    markeredgecolor='black')]

    ax2.legend(handles=line, loc=4, edgecolor='black', title='Depth (m)', title_fontsize=8,  numpoints=1, ncol=1, fontsize=7)

    ax1.set_xlabel('y (km)', fontsize=10)
    ax2.set_xlabel('y (km)',fontsize=10)
    ax1.set_ylabel('x (km)', fontsize=10)
    ax2.set_ylabel('x (km)', fontsize=10)

    ax3.set_xlabel('Depth (m)', fontsize=10)
    ax4.set_xlabel('Depth (m)', fontsize=10)
    ax3.set_ylabel('Probability density', fontsize=10)
    ax4.set_ylabel('Probability density', fontsize=10)

    ax1.tick_params(axis='both',which='major', labelsize=9)
    ax2.tick_params(axis='both',which='major', labelsize=9)
    ax3.tick_params(axis='both',which='major', labelsize=9)
    ax4.tick_params(axis='both',which='major', labelsize=9)

    ax1.yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax2.yaxis.set_ticks(np.arange(min(y)/1000, max(y)/1000 + 5, step=5))
    ax1.xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))
    ax2.xaxis.set_ticks(np.arange(min(x)/1000, max(x)/1000 + 5, step=5))

    ax3.vlines(x=100, ymin=0,ymax= 0.03, color='red', linestyles='dashed')
    ax3.vlines(x=200, ymin=0,ymax= 0.03, color='red', linestyle='dashed')
    ax4.vlines(x=100, ymin=0, ymax=0.03, color='red', linestyles='dashed')
    ax4.vlines(x=200, ymin=0, ymax=0.03, color='red', linestyle='dashed')

    ax3.text(100, 0.033, '100 (A)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')
    ax3.text(200, 0.033, '200 (B)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')
    ax4.text(100, 0.033, '100 (A)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')
    ax4.text(200, 0.033, '200 (B)', fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')

    ax3.yaxis.set_ticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax4.yaxis.set_ticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])
   
    ax3.set_xlim([0, 240])
    ax4.set_xlim([0, 240])

    ax3.xaxis.set_ticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240])
    ax4.xaxis.set_ticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240])

    ax1.text(3, 16, 'A', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax2.text(3, 16, 'A', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax1.text(9, 11, 'B', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')
    ax2.text(9, 11, 'B', color='white', fontsize=14, fontweight='bold', horizontalalignment='center', verticalalignment='center')

    ax3.legend(edgecolor='black', fontsize=9, loc='upper left')
    ax4.legend(edgecolor='black', fontsize=9, loc='upper left')

    ax1.text(-4.2, 20, 'a)', fontsize=13, horizontalalignment='center', verticalalignment='center')
    ax2.text(-4.2, 20, 'b)', fontsize=13, horizontalalignment='center', verticalalignment='center')
    ax3.text(-27, 0.053, 'c)', fontsize=13, horizontalalignment='center', verticalalignment='center')
    ax4.text(-27, 0.053, 'd)', fontsize=13, horizontalalignment='center', verticalalignment='center')

    plt.subplots_adjust(hspace=0.36, wspace=0.32)

    plt.savefig('figures\FIG4.png', bbox_inches='tight', dpi=600)
    plt.close('all')

    return
