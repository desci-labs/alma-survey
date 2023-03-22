# ======================== Import Packages ==========================

import sys, os, pdb
import numpy as np
from astropy.table import Table, vstack, Column, join
import calc_dust_masses
import csv
from lifelines import KaplanMeierFitter


# ========================== Define Fuctions ==========================

def get_usc(f1, f2):

    """
    PURPOSE:    Create table for Upper Sco with info needed to run the 
                two-sample tests and the KME for comparing regions

                Calculates dust masses in same way as in Ansdell+2016
                Replaces non-detections with 3-sigma upper limits
    INPUT:      f1 = path to table 1 from Barenfeld+2016 (2016ApJ...827..142B)
                f2 = path to table 4 from Barenfeld+2016 (2016ApJ...827..142B)

    OUTPUT:     t = table that contains:
                  Source names (str)
                  Detection flags (int; 1=detection, 0=non-detection)
                  Dust masses with non-detections set to 3-sigma upper limits
                  Continuum fluxes with non-detections set to 3-sigma upper limits
                  Stellar masses and spectral types

    """

    ### GRAB UPPER SCO DATA
    t1 = Table.read(f1, format='ascii.cds')
    t2 = Table.read(f2, format='ascii.cds')
    t = join(t1, t2, join_type='inner')
    t['Source'].name = 'Name'

    ### FLAG (NON-)DETECTIONS
    ### REPLACE NON-DETECTION FLUXES WITH 3-SIGMA UPPER LIMITS
    t['det'] = np.repeat(0, len(t))
    t['det'][t['Snu'] / t['e_Snu'] >= 3.0] = 1
    t['Snu'][np.where(t['det'] == 0)] = 3.0 * t['e_Snu'][np.where(t['det'] == 0)]
    
    ### ONLY KEEP MSTAR >= 0.1
    t = t[np.where(10**t['logM'] >= 0.1)]

    ### CALCULATE DUST MASSES USING SAME METHOD AS LUPUS
    mdust = []
    for i, val in enumerate(t):
        mdust.append(calc_dust_masses.get_dustmass(1.0, 880., t['Snu'][i], 145., 20.))
    t['MDust'] = mdust

    ### ONLY KEEP "PRIMORDIAL" DISKS 
    ### TO MATCH SAMPLES OF LUPUS & TAURUS
    t = t[np.where( (t['Type'] == 'Full') | (t['Type'] == 'Transitional') | (t['Type'] == 'Evolved'))] 

    ### CLEAN UP
    t['Snu'].name = 'FCont'
    t['Mstar'] = 10**t['logM']

    return t['Name', 'SpT', 'Mstar', 'FCont', 'MDust', 'det']


def get_lup(f1, f2):

    """
    PURPOSE:    Create table for Lupus with info needed to run the 
                two-sample tests and the KME for comparing regions

                Calculates dust masses in same way as in Ansdell+2016
                Replaces non-detections with 3-sigma upper limits

    INPUT:      f1 = path to table 1 from Ansdell+2016 (2016ApJ...828...46A)
                f2 = path to table 2 from Ansdell+2016 (2016ApJ...828...46A)

    OUTPUT:     t = table that contains:
                  Source names (str)
                  Detection flags (int; 1=detection, 0=non-detection)
                  Dust masses with non-detections set to 3-sigma upper limits
                  Continuum fluxes with non-detections set to 3-sigma upper limits
                  Stellar masses and spectral types
                  Distances
    """

    t1 = Table.read(f1, format='ascii.cds')
    t2 = Table.read(f2, format='ascii.cds')
    t = join(t1, t2, join_type='inner')

    ### ONLY KEEP KNOWN STELLAR MASSES
    t = t[~t['Mass'].mask]

    ### ONLY KEEP MSTAR >= 0.1
    t = t[t['Mass'] >= 0.1]

    ### FLAG (NON-)DETECTIONS
    ### REPLACE NON-DETECTION FLUXES WITH 3-SIGMA UPPER LIMITS
    t['det'] = np.repeat(0, len(t))
    t['det'][t['FCont'] / t['e_FCont'] >= 3.0] = 1
    t['MDust'][np.where(t['det'] == 0)] = 3.0 * t['e_MDust'][np.where(t['det'] == 0)]

    ### CLEAN UP
    t['Mass'].name = 'Mstar'
    t['SpType'].name = 'SpT'

    return t['Name', 'SpT', 'Mstar', 'FCont', 'MDust', 'det', 'Dis']


def get_tau(f1, f2, f3):

    """
    PURPOSE:    Create table for Taurus with info needed to run the 
                two-sample tests and the KME for comparing regions

                Calculates dust masses in same way as in Ansdell+2016
                Replaces non-detections with 3-sigma upper limits

    INPUT:      f1 = path to table 2 from Andrews+2013 (2013ApJ...771..129A)
                f2 = path to table 3 from Andrews+2013 (2013ApJ...771..129A)
                f2 = path to table 4 from Andrews+2013 (2013ApJ...771..129A)

    OUTPUT:     t = table that contains:
                  Source names (str)
                  Detection flags (int; 1=detection, 0=non-detection)
                  Dust masses with non-detections set to 3-sigma upper limits
                  Continuum fluxes with non-detections set to 3-sigma upper limits
                  Stellar masses and spectral types
                  Flag for observations at 890 or 1300 microns

    """

    ### GET STELLAR MASSES
    name1, mstar = [], []
    with open(f2, newline='') as f:

        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if (i >= 6) & (i <= 217):

                ### REMOVE FOOTNOTES FROM NAMES
                row[0] = row[0].split('^')[0]

                ### THESE SOURCES HAD DIFFERENT NAMES IN OTHER TABLE
                if row[0] == 'JH 112 Aa':
                    row[0] = 'JH 112 A'
                if row[0] == 'JH 112 Ab':
                    row[0] = 'JH 112 B'

                ### THESE SOURCES ARE NOT IN OTHER TABLES
                if row[0] in ['J04361030+2159364', 'J04263055+2443558', 'J04335245+2612548', 'J04290068+2755033']:
                    continue

                ### THIS SOURCE WAS SEPARATED IN A & B COMPONENTS IN OTHER TABLE
                if row[0] == 'CIDA 11 AB':
                    row[0] = ['CIDA 11 A', 'CIDA 11 B']
                    row[9] = [row[9], row[9]]
                else:
                    ### SO ALL OTHERS ARE ALSO LISTS; NEEDED FOR EXTEND FUNCTION BELOW
                    row[0] = row[0].split(';')
                    row[9] = row[9].split(';')

                ### ADD TO LIST
                name1.extend(row[0])
                mstar.extend(row[9])

    ### GET FLUXES AND SPECTRAL TYPES
    name2, spt, f890, f1300, det890, det1300, note = [], [], [], [], [], [], []
    with open(f1, newline='') as f:
        
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if (i >= 6) & (i <= 215):

                ### THESE SOURCES ARE NOT IN OTHER TABLE
                if row[0] in ['J04290068+2755033']:
                    continue

                ### GET SPECTRAL TYPE
                if '(' in row[1]:
                    row[1] = row[1].split('-')[0][1:]
                spt.append(row[1].split(', ')[0].split(' +or- ')[0])

                ### GET NAME
                name2.append(row[0].split('^')[0])

                ### GET FLUXES AND FLAG DETECTIONS
                f890.append(row[3].split(' +or- ')[0].split('<')[-1])
                if '<' in row[3]:
                    det890.append(0)
                else:
                    det890.append(1)
                f1300.append(row[4].split(' +or- ')[0].split('<')[-1])
                if '<' in row[4]:
                    det1300.append(0)
                else:
                    det1300.append(1)

                ### GET NOTE ON OBSERVED OR EXTRAPOLATED
                note.append(row[5])                

    ### CREATE TABLE
    t = Table()
    t['Name'] = name1
    t['F890'] = np.array(f890).astype(float) * 1e3
    t['F1300'] = np.array(f1300).astype(float) * 1e3
    t['det890'] = det890
    t['det1300'] = det1300
    t['Mstar'] = 10**np.array(mstar).astype(float)
    t['SpT'] = spt

    ### FIGURE OUT WHICH WERE MEASURED AT WHAT WAVELENGTH
    otype = np.empty(len(t), dtype='U10')
    for i, val in enumerate(note):
        test = val.split(", ")
        if (test[0] == 'm') or (test[0] == '(m'): 
            otype[i] = 890
        else: 
            otype[i] = 1300
    t['ObsType'] = otype

    ### ONLY KEEP MSTAR >= 0.1
    t = t[t['Mstar'] >= 0.1]

    # ### REMOVE BINARIES
    # ### SINCE LUPUS AND UPPER SCO SAMPLES DIDN'T INCLUDE SECONDARY SOURCES
    # ind_bin = []
    # for i, val in enumerate(t['Name']):
    #     if i == 0:
    #         continue
    #     if (' B' in t['Name'][i]) & (' A' in t['Name'][i-1]):
    #         ind_bin.append(i)
    #     if val == 'UZ Tau Wb':
    #         ind_bin.append(i)
    # t.remove_rows(ind_bin)

    ### CALCULATE DUST MASSES USING SAME METHOD AS Ansdell+2016
    mdust, det, fcont = [], [], []
    for i, val in enumerate(t['ObsType']):
        if val == '890':
            mdust.append(calc_dust_masses.get_dustmass(1.0, 890., t['F890'][i], 140., 20.))
            det.append(t['det890'][i])
            fcont.append(t['F890'][i])
        elif val == '1300':
            mdust.append(calc_dust_masses.get_dustmass(1.0, 1300., t['F1300'][i], 140., 20.))
            det.append(t['det1300'][i])
            fcont.append(t['F1300'][i])
    t['MDust'] = mdust
    t['det'] = det
    t['FCont'] = fcont

    return t['Name', 'SpT', 'Mstar', 'FCont', 'MDust', 'det', 'ObsType']


# ============================= Code ==================================

#### LOAD IN TABLES 
TL = get_lup('../input/apjaa2846t1_mrf.txt', '../input/apjaa2846t2_mrf.txt')
TU = get_usc('../input/apjaa2b81t1_mrt.txt', '../input/apjaa2b81t4_mrt.txt')
TT = get_tau('../input/apj476413t2_ascii.txt', '../input/apj476413t3_ascii.txt', '../input/apj476413t4_ascii.txt')

### WRITE FILES
TL.write('../output/data_lup.txt', format='ascii.ipac')
TU.write('../output/data_usc.txt', format='ascii.ipac')
TT.write('../output/data_tau.txt', format='ascii.ipac')

