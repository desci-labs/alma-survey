
# ======================== Import Packages ==========================

import sys, pdb, glob
import numpy as np
from astropy.table import Table


# ========================== Code ==========================

### READ IN FULL TABLE FROM PAPER SUPPLEMENTAL MATERIAL
T = Table.read('../input/apjaa2846t2_mrf.txt', format='ascii.cds')

### WRITE HEADER INFO
f = open('../output/table_02.tex', 'w')
f.write(r'\capstartfalse'                                  + ' \n')
f.write(r'\begin{deluxetable*}{lrrrccccr}'                 + ' \n')
f.write(r'\tabletypesize{\footnotesize}'                   + ' \n')
f.write(r'\centering'                                      + ' \n')
f.write(r'\tablewidth{500pt}'                              + ' \n')
f.write(r'\tablecaption{$890~\mu$m Continuum Properties \label{tab-cont}}' + ' \n')
f.write(r'\tablecolumns{9} '                               + ' \n')
f.write(r'\tablehead{'                                     + ' \n')
f.write(r' \colhead{Source}'                               + ' \n')
f.write(r'&\colhead{RA$_{\rm J2000}$}'                     + ' \n')
f.write(r'&\colhead{Dec$_{\rm J2000}$}'                    + ' \n')
f.write(r'&\colhead{$F_{\rm cont}$}'                       + ' \n')
f.write(r'&\colhead{rms}'                                  + ' \n')
f.write(r'&\colhead{$a$}'                                  + ' \n')
f.write(r'&\colhead{$i$}'                                  + ' \n')
f.write(r'&\colhead{PA}'                                   + ' \n')
f.write(r'&\colhead{$M_{\rm dust}$} \\'                    + ' \n')
f.write(r' \colhead{}'                                     + ' \n')
f.write(r'&\colhead{}'                                     + ' \n')
f.write(r'&\colhead{}'                                     + ' \n')
f.write(r'&\colhead{(mJy)}'                                + ' \n')
f.write(r'&\colhead{(mJy beam$^{-1}$)}'                    + ' \n')
f.write(r'&\colhead{arcsec}'                               + ' \n')
f.write(r'&\colhead{deg}'                                  + ' \n')
f.write(r'&\colhead{(deg)}'                                + ' \n')
f.write(r'&\colhead{($M_{\oplus}$)}'                       + ' \n')
f.write(r'}'                                               + ' \n')
f.write(r'\startdata'                                      + ' \n')

### WRITE FIRST 10 LINES
lim = 10
for i, val in enumerate(np.zeros(lim)):

    source = T['Name'][i]

    de = str(T['DE-'][i]) + "{0:02}".format(T['DEd'][i]) + ':' + "{0:02}".format(T['DEm'][i]) + ':' + "{0:05.2f}".format(T['DEs'][i])
    ra = "{0:02}".format(T['RAh'][i]) + ':' + "{0:02}".format(T['RAm'][i]) + ':' + "{0:05.2f}".format(T['RAs'][i])

    flx = "{0:.2f}".format(T['FCont'][i]) + r' $\pm$ ' + "{0:.2f}".format(T['e_FCont'][i])
    rms = "{0:.2f}".format(T['rms'][i])
    mdust = str(T['MDust'][i]) + r' $\pm$ ' + str(T['e_MDust'][i])

    if np.ma.is_masked(T['a'][i]) is True:
        bmaj = '...'
    else:
        bmaj  = str(T['a'][i]) + r' $\pm$ ' + str(T['e_a'][i])

    if np.ma.is_masked(T['i'][i]) is True:
        incl = '...'
    elif np.ma.is_masked(T['e_i'][i]) is True:
        incl = str(T['i'][i])
    else:
        incl  = str(T['i'][i]) + r' $\pm$ ' + str(T['e_i'][i])

    if np.ma.is_masked(T['PosAng'][i]) is True:
        pang = '...'
    elif np.ma.is_masked(T['e_PosAng'][i]) is True:
        pang = str(T['PosAng'][i])
    else:
        pang  = str(T['PosAng'][i]) + r' $\pm$ ' + str(T['e_PosAng'][i])

    ### END OF LINE
    if (i < lim - 1):
        end = r' \\' + '\n'
    else:
        end = '\n'

    f.write(source + ' & ' + ra + ' & ' + de + ' & ' + flx + ' & ' + rms + ' & ' + 
            bmaj + ' & ' + incl + ' & ' + pang + ' & ' + mdust + end)

f.write(r'\enddata' + ' \n')
f.write(r'\tablenotetext{}{(This table is available in its entirety in machine-readible form.)}' + ' \n')
f.write(r'\end{deluxetable*}' + ' \n')
f.write(r'\capstartfalse' + ' \n')
f.close()

