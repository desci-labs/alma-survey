"""
NOTE: units for "e_Mass" in "t1_mrf.txt" should be changed to "solMass"
      for compatibility with astropy.Table.read "ascii.cds" format

"""

# ======================== Import Packages ==========================

import sys, pdb, glob
import numpy as np
from astropy.table import Table


# ========================== Code ==========================

### READ IN FULL TABLE FROM PAPER SUPPLEMENTAL MATERIAL
T = Table.read('../input/apjaa2846t1_mrf.txt', format='ascii.cds')

### WRITE HEADER INFO
f = open('../output/table_01.tex', 'w')
f.write(r'\capstartfalse'                                     + ' \n')
f.write(r'\begin{deluxetable}{lrrrr}'                         + ' \n')
f.write(r'\tabletypesize{\footnotesize}'                      + ' \n')
f.write(r'\centering'                                         + ' \n')
f.write(r'\tablewidth{240pt}'                                 + ' \n')
f.write(r'\tablecaption{Stellar Properties \label{tab-star}}' + ' \n')
f.write(r'\tablecolumns{5} '                                  + ' \n')
f.write(r'\tablehead{'                                        + ' \n')
f.write(r' \colhead{Source}'                                  + ' \n')
f.write(r'&\colhead{$d$ (pc)}'                                + ' \n')
f.write(r'&\colhead{SpT}'                                     + ' \n')
f.write(r'&\colhead{$M_{\ast}$/$M_{\odot}$}'                  + ' \n')
f.write(r'&\colhead{Ref}'                                     + ' \n')
f.write(r'}'                                                  + ' \n')
f.write(r'\startdata'                                         + ' \n')

### WRITE FIRST 10 LINES
lim = 10
for i, val in enumerate(np.zeros(lim)):

    source = T['Name'][i]
    dist = str(T['Dis'][i])

    if np.ma.is_masked(T['Ref'][i]) is True:
        ref = '...'
    else:
        ref = str(T['Ref'][i])

    if np.ma.is_masked(T['SpType'][i]) is True:
        spt = '...'
    else:
        spt = T['SpType'][i]

    if np.ma.is_masked(T['Mass'][i]) is True:
        mstar = '...'
    else:
        mstar  = str(T['Mass'][i]) + r' $\pm$ ' + str(T['e_Mass'][i])

    ### END OF LINE
    if (i < lim - 1):
        end = r' \\' + '\n'
    else:
        end = '\n'

    f.write(source + ' & ' + dist + ' & ' + spt + ' & ' + mstar + ' & ' + ref + end)

### WRITE FOOTER INFO
f.write(r'\enddata' + ' \n')
f.write(r'\tablenotetext{}{References:' +
        r' (1) \cite{2014A&A...561A...2A},' +
        r' (2) Alcal\'{a} et al. (in prep),' +
        r' (3) \cite{2013MNRAS.429.1001A},' +
        r' (4) \cite{2011MNRAS.418.1194M},' +
        r' (5) \cite{2008ApJS..177..551M},' +
        r' (6) Cleeves et al. (in prep),' +
        r' (7) \cite{2015A&A...578A..23B},' +
        r' (8) \cite{2008hsf2.book..295C}.' +
        r' (This table is available in its entirety in machine-readible form.)}' + ' \n')
f.write(r'\end{deluxetable}' + ' \n')
f.write(r'\capstartfalse' + ' \n')
f.close()