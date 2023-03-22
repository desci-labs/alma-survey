# ======================== Import Packages ==========================

import os,sys,pdb,glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, Column, join
import matplotlib as mpl
from astropy import constants as const
from astropy import units as u
from scipy.stats import spearmanr, linregress
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch


# ========================== Define Functions ==========================

def get_model_grid_old(f):

    ### LOAD FILE FROM WILLIAMS & BEST 2014 (2014ApJ...788...59W)
    gasgrid = Table.read(f, format='ascii.csv')

    ### LOAD RELEVANT MODEL OUTPUTS
    gasgrid['M_gas'].name = 'Mgas'
    gasgrid['f_3-2_13co'].name = 'F13CO32'
    gasgrid['f_3-2_c18o'].name = 'FC18O32'
    gasgrid['f_3-2_c18o_low'].name = 'FC18O32l'

    ### REMOVE MASKED VALUES (1 BAD END-OF-LINE IN GRID FROM JPW)
    gasgrid = gasgrid[np.where(~gasgrid['F13CO32'].mask)]

    return gasgrid


def get_model_grid(f):

    """
    PURPOSE:    Load gas model grid
                Remove unnecessary model grid points

    INPUT:      Table 3 from Williams & Best 2014 (2014ApJ...788...59W)

    OUTPUT:     Gas model grid

    """

    ### LOAD FILE FROM WILLIAMS & BEST 2014 (2014ApJ...788...59W)
    gg = Table.read(f, format='ascii.cds')

    ### ONLY KEEP RELEVANT MODEL OUTPUTS TO REDUCE SIZE FOR FITTING
    gg = gg[np.where(gg['gamma'] == 0.)]
    gg = gg[np.where(gg['Mgas'] <= 0.03)]

    return gg


def get_cal_err(f, e, e_cal=0.10):

    """
    PURPOSE:    Calculate multiplication factor for adding 
                calibration error to flux measurement error
                (assumed to be 10% unless otherwise specified)

    INPUT:      f = measured flux (float)
                e = measurement error (float)
                e_cal = calibration error fraction (float; optional)

    OUTPUT:     Multiplication factor for adding calibration error

    """

    mult = np.sqrt((f * e_cal)**2 + (e)**2)

    return mult


def get_model_idx(g13, g18_ism, g18_lo, f13, e13, f18, e18, d13, d18):

    """
    PURPOSE:    Index model grid points for a given flux measurement 

    INPUT:      g13 = model grid points for 13CO flux
                g18_hi = model grid points for C18O flux (ISM abundance)
                g18_lo = model grid points for C18O flux (low abundance)
                f13 = 13CO flux measurement (float)
                e13 = 13CO flux measurement error(float)
                f18 = C18O flux measurement (float)
                e18 = C18O flux measurement error(float)
                d13 = detection flags for 13CO flux (masked array)
                d18 = detection flags for C18O flux (masked array)

    OUTPUT:     ifit = indexes of model grid
                inote = note indicating type of detection

    """

    ### IF BOTH LINES DETECTED DETECTED
    ### INDEX GRID WITHIN ERRORS (MEASUREMENT + FLUX CAL)
    if (d13 == True) and (d18 == True):
        i13 = ( abs(g13 - f13)  < get_cal_err(f13, e13) )
        i18_ism = ( abs(g18_ism - f18)  < get_cal_err(f18, e18) )
        i18_lo = ( abs(g18_lo - f18) < get_cal_err(f18, e18) )
        inote = "GD"

    ### IF ONLY 13CO DETECTED
    ### INDEX GRID WITHIN UPPER LIMIT FOR C180
    elif (d13 == True) and (d18 == False):
        i13 = (abs(g13 - f13) < get_cal_err(f13, e13) )
        i18_ism = (g18_ism < f18)
        i18_lo = (g18_lo < f18)
        inote = "D13"

    ### IF BOTH LINES UNDETECTED
    ### INDEX GRID WITHIN UPPER LIMITS FOR BOTH
    elif (d13 == False) and (d18 == False):
        i13 = (g13 < f13)
        i18_ism = (g18_ism < f18)
        i18_lo = (g18_lo < f18)
        inote = "ND"

    ### STEP CODE IF UNKNOWN RESULT
    else:
        print("Unknown result")
        pdb.set_trace()

    ### COMBINING ISM & LOW C18O ABUNDANCES
    ### CONSERVATIVE SINCE CONSIDERING "ALL" POSSIBLE CO ABUNDANCES
    i18 = i18_ism + i18_lo

    ### KEEP ONLY THOSE IN BOTH LINES
    ### I.E., CREATING BOX AROUND MEASUREMENT OR UPPER LIMIT
    ifit = i13 & i18

    return ifit, inote


def get_model_gasmass(gm, ifit, inote):

    """
    PURPOSE:    Get gas mass from model grid 

    INPUT:      gm = all gas masses from model grid 
                ifit = model grid indexes for a given flux measurement
                inote = note indicating type of detection

    OUTPUT:     mgas_fit = gas mass based on model fit
                mgas_min = lower limti of gas mass based on model fit
                mgas_max = upper limit of gas mass based on model fit
                mgas_note = note indicating type of gas mass estimate
    """

    nfit = np.count_nonzero(ifit)
    if (nfit > 0):
        
        mgasfit = gm[ifit]

        ### ONLY 13CO DETECTED
        if (d13==True) and (d18==True):
            mgas_fit = 10**np.mean(np.log10(mgasfit))
            mgas_min, mgas_max = np.min(mgasfit), np.max(mgasfit)
            mgas_note = "GF"

        ### BOTH LINES DETECTED
        elif (d13==True) and (d18==False):
            mgas_fit = 10**np.mean(np.log10(mgasfit))
            mgas_min, mgas_max = 0.0, np.max(mgasfit)
            mgas_note = "GL"

        ### BOTH LINES UNDETECTED
        elif (d13==False) and (d18==False):
            mgas_fit = -99.0
            mgas_min, mgas_max = 0.0, np.max(mgasfit)
            mgas_note = "UL"

        ### STOP CODE IF UNKNOWN RESULT
        else:
            print("Unknown flux measurement result")
            pdb.set_trace()
                  
    else:

        ### STOP CODE IF NO MATCHES TO MODEL GRID
        print("No matches to model grid")
        pdb.set_trace()

    ### COMBINE NOTES
    mgas_note = (', ').join([inote, mgas_note])
                            
    return mgas_fit, mgas_min, mgas_max, mgas_note


def get_gasmass(g, f13, e13, d13, f18, e18, d18):

    """
    PURPOSE:    Calculate gas mass

    INPUT:      g = model grid 
                f13 = 13CO flux measurement (float)
                e13 = 13CO flux measurement error(float)
                d13 = detection flags for 13CO flux (masked array)
                f18 = C18O flux measurement (float)
                e18 = C18O flux measurement error(float)
                d18 = detection flags for C18O flux (masked array)

    OUTPUT:     mg_fit = gas mass based on model fit
                mgas_min = lower limti of gas mass based on model fit
                mgas_max = upper limit of gas mass based on model fit
                mgas_note = note indicating type of gas mass estimate
    """

    ### INDEX MODEL GRID FOR THIS FLUX MEASUREMENT
    i_fit, i_note = get_model_idx(g['F13CO32'], g['FC18O32'], g['FC18O32l'], 
                                  f13, e13, f18, e18, d13, d18)

    ### CALCULATE GAS MASS
    mgas_fit, mgas_lo, mgas_hi, mgas_note = get_model_gasmass(g['Mgas'], i_fit, i_note)
    
    return mgas_fit, mgas_lo, mgas_hi, mgas_note


# ============================= Code ==================================

### GET GAS MEASUREMENTS FROM ANSDELL+2016
TS = Table.read('../input/apjaa2846t1_mrt.txt', format='ascii.cds')
TG = Table.read('../input/apjaa2846t3_mrf.txt', format='ascii.cds')
T = join(TS, TG, join_type='inner')

### LOAD GAS MODEL GRID FROM WILLIAMS & BEST 2014 (2014ApJ...788...59W)
# G = get_model_grid('../input/apj495435t3_mrt.txt')
G = get_model_grid_old('../input/gasgrid.csv')

### GET GAS MASSES
mg_f, mg_m, mg_l, mg_h, mg_n = [], [], [], [], []
for i, val in enumerate(T['Name']):

    ### GET GAS FLUXES IN JY, SCALED TO 140 PC TO MATCH MODEL GRID
    f2l = (T['Dis'][i] / 140.)**2 / 1000.0
    f13, e13 = f2l * T['F13CO'][i], f2l * T['e_F13CO'][i]
    f18, e18 = f2l * T['F18CO'][i], f2l * T['e_F18CO'][i]

    ### FLAG (NON-)DETECTIONS 
    d13, d18 = T['l_F13CO'].mask[i], T['l_F18CO'].mask[i]

    ### CALCULATE GAS MASSES
    mgf, mgl, mgh, mgn = get_gasmass(G, f13, e13, d13, f18, e18, d18)

    ### SAVE GAS MASSES (M_JUP), LIMITS, AND FLAGS
    if mgn == 'ND, UL':
        mg_n.append('<')
        mg_f.append(round(mgh*(const.M_sun.cgs/const.M_jup.cgs).value, 3))
        mg_l.append(-99.)
        mg_h.append(-99.)

    elif mgn == 'D13, GL':
        mg_n.append('--')
        mg_f.append(round(mgf*(const.M_sun.cgs/const.M_jup.cgs).value, 3))
        mg_h.append(round(mgh*(const.M_sun.cgs/const.M_jup.cgs).value, 3))
        mg_l.append(-99.)

    elif mgn == 'GD, GF':
        mg_n.append('--')
        mg_f.append(round(mgf*(const.M_sun.cgs/const.M_jup.cgs).value, 3))
        mg_l.append(round(mgl*(const.M_sun.cgs/const.M_jup.cgs).value, 3))
        mg_h.append(round(mgh*(const.M_sun.cgs/const.M_jup.cgs).value, 3))

    else:
        print("Unknown gas mass result")
        pdb.set_trace()








    
