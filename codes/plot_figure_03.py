
# ======================== Import Packages ==========================

import os,sys,pdb
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join, MaskedColumn
from astropy import constants as const
import csv
import matplotlib as mpl
from astropy import constants as const


# ============================= Code ==================================

#### LOAD IN TABLES FROM PAPER SUPPLEMENTAL MATERIAL
Table_Dust = Table.read('../input/apjaa2846t2_mrf.txt', format='ascii.cds')
Table_Gas = Table.read('../input/apjaa2846t3_mrf.txt', format='ascii.cds')
Table_All = join(Table_Dust, Table_Gas, join_type='inner')

### ONLY KEEP DUST DETECTIONS, CHANGE UNITS TO SOLAR MASSES, AND SORT TABLE BY INCREASING DUST MASS
Table_All = Table_All[Table_All['FCont']/Table_All['e_FCont'] > 3.0]
Table_All['MDust'] = Table_All['MDust']*(const.M_earth.cgs/const.M_sun.cgs).value
Table_All['e_MDust']  = np.sqrt( (Table_All['e_MDust'])**2 + (0.1*Table_All['e_MDust'])**2 )*(const.M_earth.cgs/const.M_sun.cgs).value
Table_All.sort('MDust')

### CHANGE GAS MASS UNITS TO SOLAR MASSES
Table_All['Mgas'] = Table_All['Mgas'] * (const.M_jup.cgs/const.M_sun.cgs).value
Table_All['b_Mgas'] = Table_All['b_Mgas'] * (const.M_jup.cgs/const.M_sun.cgs).value
Table_All['B_Mgas'] = Table_All['B_Mgas'] * (const.M_jup.cgs/const.M_sun.cgs).value

### SETUP PLOT
mpl.rc('xtick', labelsize=15) 
mpl.rc('ytick', labelsize=15)
mpl.rc('xtick.major',size=7, pad=7, width=1)
mpl.rc('ytick.major',size=7, pad=7, width=1)
mpl.rc('axes', linewidth=1)
mpl.rc('lines', markersize=5)
fig = plt.figure(figsize = (16, 12))
    
### SET UP DUST MASS SUBPLOT
ax1 = fig.add_subplot(311)
ax1.set_ylim(1e-7, 1e-3)
ax1.set_yscale('log')
ax1.set_ylabel('Dust Mass ('+r'$\mathregular{M_{\odot}}$'+')', fontsize=17)
ax1.set_xlim(-1, len(Table_All))
ax1.xaxis.set_ticks(np.arange(-4, len(Table_All),1))
fig.canvas.draw()
ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax1.tick_params(axis='y', which='both', labelsize=12)
ax1.grid(color='gray',linewidth=1, linestyle=':', alpha=0.3)
plt.axhline(y=3e-06,color='gray',linewidth=1.0)
ax1.text(-3.5,4e-06,r'$\mathregular{1M_{\oplus}}$',fontsize=12,style='italic',color='gray')
plt.axhline(y=3e-05,color='gray',linewidth=1.0)
ax1.text(-3.5,4e-05,r'$\mathregular{10M_{\oplus}}$',fontsize=12,style='italic',color='gray')
ax1.errorbar(-1, 5.68e-05, fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax1.errorbar(-2, 1.20e-05, fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax1.errorbar(-3, 2.26e-07, fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax1.arrow(-3, 2.26e-07, 0.0, -0.4*2.26e-07, head_width=0.3, head_length=2.26e-07*0.1, fc='black', ec='black', linewidth=1, zorder=1)

### SET UP GAS MASS SUBPLOT
ax2 = fig.add_subplot(312)
ax2.set_ylim(1e-5, 1e-1)
ax2.set_yscale('log')
ax2.set_ylabel('Gas Mass ('+r'$\mathregular{M_{\odot}}$'+')', fontsize=17)
ax2.set_xlim(-1,len(Table_All))
ax2.xaxis.set_ticks(np.arange(-4, len(Table_All), 1))
ax2.tick_params(axis='y', which='both', labelsize=12)
fig.canvas.draw()
ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax2.grid(color='gray', linewidth=1, linestyle=':', alpha=0.3)
plt.axhline(y=1e-3, color='gray', linewidth=1.0)
plt.axhline(y=1e-2, color='gray', linewidth=1.0)
ax2.text(-3.5, 1.4e-3, r'$\mathregular{1M_{Jup}}$', fontsize=12, style='italic', color='gray')
ax2.text(-3.5, 1.2e-2,r'$\mathregular{MMSN}$', fontsize=12, style='italic', color='gray')
ax2.errorbar(-1, 3.8e-4, fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax2.errorbar(-2, 1.9e-4, fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax2.arrow(-2, 1.9e-4, 0.0, -0.5 * 1.9e-4, head_width=0.3, head_length=1.9e-4 * 0.1, fc='black', ec='black', linewidth=1, zorder=1)

### SET UP GAS-TO-DUST MASS SUBPLOT
ax3 = fig.add_subplot(313)
ax3.set_ylim(0.1, 1e4)
ax3.set_yscale('log')
ax3.set_ylabel('Gas-to-Dust Ratio', fontsize=17) 
ax3.set_xlim(-1,len(Table_All))
ax3.xaxis.set_ticks(np.arange(-4, len(Table_All), 1))
ax3.tick_params(axis='y', which='both', labelsize=12)
fig.canvas.draw()
labels = list(np.append([' ', 'Continuum non-det', r'$\mathregular{^{13}\!\hspace{0.1}CO}$' + ' & ' + r'$\mathregular{C^{18}\!\hspace{0.1}O}$' + 
                         ' non-det', r'$\mathregular{C^{18}\!\hspace{0.1}O}$' + ' non-det'], Table_All['Name']))
ax3.set_xticklabels(labels, rotation='vertical', size=9)
ax3.grid(color='gray',linewidth=1, linestyle=':', alpha=0.3)
plt.axhline(y=100, color='gray', linewidth=1.0)
ax3.text(-3.5, 1.3e2, 'ISM', fontsize=11, style='italic', color='gray')
ax3.errorbar(-1, 7, fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax3.errorbar(-2, 13,fmt='*', color='white', ms=12, mec='black', mew=1, zorder=999)
ax3.arrow(-2, 13, 0.0, -0.5 * 13, head_width=0.3, head_length=13 * 0.1, fc='black', ec='black', linewidth=1, zorder=1)

### PLOT DUST MASSES
ax1.errorbar(np.arange(len(Table_All)), Table_All['MDust'], yerr=Table_All['e_MDust'], 
             fmt='o', mfc='lightblue', ms=6, mec='black', mew=0.8, ecolor='black', elinewidth=1, zorder=999,capsize=2)

### PLOT GAS MASSES AND GAS-TO-DUST RATIOS
for i, val in enumerate(Table_All['Name']):

    ### PLOT NON-DETECTIONS
    if Table_All['l_Mgas'][i] == '<':

        ax2.errorbar(i, Table_All['Mgas'][i], fmt='v', color='lightgray', ms=7, mec='black', mew=0.9, zorder=999)
        ax3.errorbar(i, int(round(Table_All['Mgas'][i]/Table_All['MDust'][i])), fmt='v', color='lightgray', ms=7, mec='black', mew=0.9, zorder=999)

    ### PLOT DETECTIONS WITH ONLY UPPER LIMITS
    elif Table_All['b_Mgas'][i] is np.ma.masked:

        ax2.arrow(i, Table_All['Mgas'][i], 0.0, -0.5 * Table_All['Mgas'][i], head_width=0.3, head_length=Table_All['Mgas'][i] * 0.1,
                 fc='black', ec='black', linewidth=0.9, zorder=1)
        ax2.errorbar(i, Table_All['Mgas'][i], yerr=[[0],[Table_All['B_Mgas'][i] - Table_All['Mgas'][i]]],
                    fmt='o', mfc='lightblue', ms=7, mec='black', mew=1, ecolor='black', elinewidth=0.9, zorder=999, capsize=3)

        g2d = int(round(Table_All['Mgas'][i]/Table_All['MDust'][i]))
        g2d_max = int(round(Table_All['B_Mgas'][i]/Table_All['MDust'][i]))
        ax3.arrow(i, g2d, 0.0, -0.5 * g2d, head_width=0.3, head_length=g2d * 0.1, fc='black', ec='black', linewidth=1, zorder=1)
        ax3.errorbar(i, g2d,yerr=[[0],[g2d_max - g2d]], fmt='o', mfc='lightblue', ms=7, mec='black', mew=1, 
                     ecolor='black', elinewidth=1, zorder=999, capsize=3)

    ### PLOT DETECTIONS WITH BOTH UPPER AND LOWER LIMITS
    else:
        
        ax2.errorbar(i, Table_All['Mgas'][i], yerr=[[Table_All['Mgas'][i] - Table_All['b_Mgas'][i]], [Table_All['B_Mgas'][i] - Table_All['Mgas'][i]]],
                    fmt='o', mfc='lightblue', ms=6, mec='black', mew=1, ecolor='black', elinewidth=0.9, zorder=999, capsize=3)
        
        g2d = int(round(Table_All['Mgas'][i]/Table_All['MDust'][i]))
        g2d_max = int(round(Table_All['B_Mgas'][i]/Table_All['MDust'][i]))
        g2d_min = int(round(Table_All['b_Mgas'][i]/Table_All['MDust'][i]))
        ax3.errorbar(i, g2d, yerr=[[g2d - g2d_min],[g2d_max - g2d]], fmt='o', mfc='lightblue', 
                     ms=6, mec='black', mew=0.9, ecolor='black', elinewidth=1, zorder=999, capsize=3)
        
### SAVE PLOT
fig.savefig('../output/figure_03.pdf', bbox_inches='tight', dpi=100)
plt.close('all')