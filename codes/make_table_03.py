# ======================== Import Packages ==========================

import sys, pdb, glob
import numpy as np
from astropy.table import Table


# ========================== Code ==========================

### READ IN FULL TABLE FROM PAPER SUPPLEMENTAL MATERIAL
T = Table.read('../input/apjaa2846t3_mrf.txt', format='ascii.cds')

### WRITE HEADER INFO
f = open('../output/table_03.tex', 'w')
f.write(r'\capstartfalse'                                  + ' \n')
f.write(r'\begin{deluxetable*}{lrrrrrrr}'                    + ' \n')
f.write(r'\tabletypesize{\footnotesize}'                   + ' \n')
f.write(r'\centering'                                      + ' \n')
f.write(r'\tablewidth{500pt}'                              + ' \n')
f.write(r'\tablecaption{Gas Properties \label{tab-gas}}'   + ' \n')
f.write(r'\tablecolumns{8} '                               + ' \n')
f.write(r'\tablehead{'                                     + ' \n')
f.write(r' \colhead{Source}'                               + ' \n')
f.write(r'&\colhead{$F_{\rm 13CO}$}'                       + ' \n')
f.write(r'&\colhead{$E_{\rm 13CO}$}'                       + ' \n')
f.write(r'&\colhead{$F_{\rm C18O}$}'                       + ' \n')
f.write(r'&\colhead{$E_{\rm C18O}$}'                       + ' \n')
f.write(r'&\colhead{$M_{\rm gas}$}'                        + ' \n')
f.write(r'&\colhead{$M_{\rm gas,min}$}'                    + ' \n')
f.write(r'&\colhead{$M_{\rm gas,max}$} \\'                 + ' \n')
f.write(r' \colhead{}'                                     + ' \n')
f.write(r'&\colhead{(mJy~km~s$^{-1}$)}'                    + ' \n')
f.write(r'&\colhead{(mJy~km~s$^{-1}$)}'                    + ' \n')
f.write(r'&\colhead{(mJy~km~s$^{-1}$)}'                    + ' \n')
f.write(r'&\colhead{(mJy~km~s$^{-1}$)}'                    + ' \n')
f.write(r'&\colhead{($M_{\rm Jup}$)}'                      + ' \n')
f.write(r'&\colhead{($M_{\rm Jup}$)}'                      + ' \n')
f.write(r'&\colhead{($M_{\rm Jup}$)}'                      + ' \n')
f.write(r'}'                                               + ' \n')
f.write(r'\startdata'                                      + ' \n')


lim = 10
for i, val in enumerate(np.zeros(lim)):

    source = T['Name'][i]

    if np.ma.is_masked(T['l_Mgas'][i]) is True:
        mgas = "{0:01}".format(T['Mgas'][i])

    else:
        mgas = r'$' + str(T['l_Mgas'][i]) + "{0:01}".format(T['Mgas'][i]) + r'$'

    if np.ma.is_masked(T['b_Mgas'][i]) is True:
        b_mgas = '...'
    else:
        b_mgas = r'$' + "{0:01}".format(T['b_Mgas'][i]) + r'$'

    if np.ma.is_masked(T['B_Mgas'][i]) is True:
        B_mgas = '...'
    else:
        B_mgas = r'$' + "{0:01}".format(T['B_Mgas'][i]) + r'$'


    if np.ma.is_masked(T['l_F13CO'][i]) is True:
        f13 = "{0:0}".format(T['F13CO'][i])
        e13 = "{0:0}".format(T['e_F13CO'][i])
    else:
        f13 = r'$' + str(T['l_F13CO'][i]) + "{0:0}".format(T['F13CO'][i]) + r'$'
        e13 = '...'


    if np.ma.is_masked(T['l_F18CO'][i]) is True:
        f18 = "{0:0}".format(T['F18CO'][i])
        e18 = "{0:0}".format(T['e_F18CO'][i])
    else:
        f18 = r'$' + str(T['l_F18CO'][i]) + "{0:0}".format(T['F18CO'][i]) + r'$'
        e18 = '...'


    if (i < lim - 1):
        end = r' \\' + '\n'
    else:
        end = '\n'

    f.write(source + ' & ' + f13 + ' & ' + e13 + ' & ' + f18 + ' & ' + e18 + ' & ' + 
            mgas + ' & ' + b_mgas + ' & ' + B_mgas + end)

f.write(r'\enddata' + ' \n')
f.write(r'\tablenotetext{}{(This table is available in its entirety in machine-readible form.)}' + ' \n')
f.write(r'\end{deluxetable*}' + ' \n')
f.write(r'\capstartfalse' + ' \n')
f.close()