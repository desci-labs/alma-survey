
# ======================== Import Packages ==========================

import sys, pdb, glob
import numpy as np
from astropy.table import Table
import csv

# ========================== Code ==========================

### READ IN TABLE FROM PAPER SUPPLEMENTAL MATERIAL
with open('../input/apjaa2846t4_ascii.txt', newline='') as f:
    name, ra, de, cont, gas13, gas18 = [],[],[],[],[],[]
    reader = csv.reader(f, delimiter='\t')
    for i, row in enumerate(reader):            
            if i >= 6:
                name.append(row[0])
                ra.append(row[1])
                de.append(row[2])
                cont.append(row[3])
                gas13.append(row[4])
                gas18.append(row[5])


### WRITE HEADER INFO
f = open('../output/table_04.tex', 'w')
f.write(r'\begin{table}[!htb]'                        + ' \n')
f.write(r'\caption{Rejected Targets}'                 + ' \n')
f.write(r'\label{tab-rejected}'                       + ' \n')
f.write(r'\centering  '                               + ' \n')
f.write(r'\begin{tabular}{lccrcc}'                    + ' \n')
f.write(r'\hline\hline'                               + ' \n')
f.write(r'Source & RA$_{\rm J2000}$ & Dec$_{\rm J2000}$ & $F_{\rm cont}$ & $F_{\rm 13CO}$} & $F_{\rm C18O}$ \\' + ' \n')
f.write(r' &  &  & (mJy) & (mJy~km~s$^{-1}$) & (mJy~km~s$^{-1}$) \\' + ' \n')
f.write(r'\hline' + ' \n')

### WRITE DATA
for i, val in enumerate(name):

    ### END OF LINE
    if (i < len(name)):
        end = r' \\' + '\n'
    else:
        end = '\n'

    f.write(name[i] + ' & ' + ra[i] + ' & ' + de[i] + ' & ' + cont[i].replace('+or-', '$\pm$') + ' & ' +  
            gas13[i].replace('\\lt', '<') + ' & ' + gas18[i].replace('\\lt', '<') + end)


f.write(r'\hline'         + ' \n')
f.write(r'\end{tabular}'  + ' \n')
f.write(r'\end{table}'   + ' \n')
f.close()

