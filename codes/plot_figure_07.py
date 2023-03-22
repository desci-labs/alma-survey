
# ======================== Import Packages ==========================

import os,sys,pdb
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join, MaskedColumn
from astropy import constants as const
import csv
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator


# ============================= Code ==================================

#### LOAD IN TABLES FROM PAPER SUPPLEMENTAL MATERIAL
TS = Table.read('../input/apjaa2846t1_mrf.txt', format='ascii.cds')
TG = Table.read('../input/apjaa2846t3_mrf.txt', format='ascii.cds')
T = join(TS, TG, join_type='inner')

### REMOVE SOURCES WITH UNKNOWN STELLAR MASSES
T = T[~T['Mass'].mask]

### DEFINE PLOT VARIABLES
x_data = np.ma.log10(T['Mass'])
x_err  = np.ma.array(0.434*(T['e_Mass'] / T['Mass']))
y_data = np.ma.log10(T['Mgas'])
y_max  = np.ma.log10(T['B_Mgas'])
y_min  = np.ma.log10(T['b_Mgas'])

### SETUP PLOT
mpl.rc('xtick', labelsize=15) 
mpl.rc('ytick', labelsize=15)
mpl.rc('xtick.major', size=7, pad=7, width=1)
mpl.rc('ytick.major', size=7, pad=7, width=1)
mpl.rc('axes', linewidth=1)
mpl.rc('lines', markersize=5)
fig = plt.figure(figsize = (8, 6))
ax = fig.add_subplot(111)
ax.set_xlim([-1.3, 0.7])
ax.set_ylim([-2.0, 2.0])
ax.yaxis.set_major_locator(MaxNLocator(5, integer=True))
ax.xaxis.set_major_locator(MaxNLocator(4))
ax.set_xlabel(r'$\mathregular{log(M_{\ast})}$' + ' ' + r'$\mathregular{[M_{\odot}]}$', fontsize=17)
ax.set_ylabel(r'$\mathregular{log(M_{gas})}$' + ' ' + r'$\mathregular{[M_{Jup}]}$', fontsize=17)
ax.tick_params(which='minor', axis='x', length=3, color='k', width=1)
ax.tick_params(which='minor', axis='y', length=3, color='k', width=1)
ax.minorticks_on()

### PLOT DATA 
for i, val in enumerate(T['Name']):

    ### PLOT NON-DETECTIONS
    if T['l_Mgas'][i] == '<':
        ax.errorbar(x_data[i], y_data[i], fmt='v', color='lightgray', ms=6, mec='black', mew=0.9, zorder=2)

    ### PLOT DETECTIONS WITH ONLY UPPER LIMITS
    elif T['b_Mgas'][i] is np.ma.masked:
        ax.scatter(x_data[i], y_data[i], marker='o', facecolor='lightblue', s=60, edgecolor='darkslategray', zorder=3, linewidth=1)
        ax.arrow(x_data[i], y_data[i], 0.0, -0.3, head_width=0.02, head_length=0.04, fc='lightgray', ec='lightgray', linewidth=1, zorder=1)
        ax.errorbar(x_data[i], y_data[i], yerr=[[0], [y_max[i] - y_data[i]]], xerr=[x_err[i]],
                    fmt='o', mfc='lightblue', ms=0, mec='black', mew=1, ecolor='lightgray', elinewidth=1, zorder=1, capsize=3)
        
    ### PLOT DETECTIONS WITH BOTH UPPER AND LOWER LIMITS
    else:        
        ax.errorbar(x_data[i], y_data[i], yerr=[[y_data[i] - y_min[i]], [y_max[i] - y_data[i]]], xerr=[x_err[i]],
                    fmt='o', mfc='lightblue', ms=0, mec='black', mew=1, ecolor='lightgray', elinewidth=1, zorder=1, capsize=3)
        ax.scatter(x_data[i], y_data[i], marker='o', facecolor='lightblue', s=60, edgecolor='darkslategray', zorder=3, linewidth=1)

### SAVE FIGURE
fig.savefig('../output/figure_07.pdf', bbox_inches='tight', dpi=100)
plt.close('all')