# ======================== Import Packages ===========================

import sys,os,pdb,glob
import numpy as np
from astropy.io import ascii
from scipy.interpolate import interp1d
from astropy.table import Table, join
from scipy.stats import norm 
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import scipy.ndimage
import csv


# ========================== Define Fuctions =========================

def norm_flux(flux, dist, freq):

    """
    PURPOSE:    Normalize fluxes to common distance and frequency

    INPUT:      flux = fluxes in mJy (array of floats)
                dist = distances in pc (array of floats)
                freq = frequency in GHz (array of floats)

    OUTPUT:     flux_norm = normalized fluxes (float)

    """

    ### REFERENCE PARAMETERS
    ### (SET FOR LUPUS)
    freq_ref = 335.8 # GHz
    dist_ref = 150.0 # pc
    beta_ref = 1.
    
    ### NORMALIZE TO REFERENCE FREQUENCY AND DISTANCE
    flux_norm = np.empty(len(flux), dtype=float)    
    for i, val in enumerate(flux):
        tmp = flux [i] * (freq_ref / freq[i])**(2.0 + beta_ref)
        tmp = tmp * ((dist[i] / dist_ref)**2)
        flux_norm[i] = float("{0:.4f}".format(tmp))

    return flux_norm


def get_nspt(sspt):

    """
    PURPOSE:    Converting string spectral type (e.g., M5) to 
                numerical spectral type (e.g., 66)

    INPUT:      sspt = string spectral type (str)

    OUTPUT:     nspt = numerical spectral type (float)
    
    """
    
    cl_conv = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L']
    nu_conv = [  0,  10,  20,  30,  40,  50,  60,  70]
    
    nspt = np.zeros(len(sspt))
    for i in np.arange(len(sspt)):
        nspt[i] = nu_conv[cl_conv.index(sspt[i][0])] + float(sspt[i][1:])

    return nspt


def censored_twosample(distA, censA, distB, censB):

    """
    PURPOSE:    Calculate a two-sample test for censored datasets following 
                methodology of Feigelson & Nelson 1985 (1985ApJ...293..192F)

    INPUT:      distA = distribution of Sample A
                censA = censoring of Sample A
                distB = distribution of Sample B
                censB = censoring of Sample B

    OUTPUT:     output = probability that the incomplete and/or biased comparison samples 
                are drawn from the same parent population as the complete reference sample

    """
    

    ### SWAP CENSORING DIRECTION
    vA = np.max(np.concatenate((distA, distB))) - distA
    vB = np.max(np.concatenate((distA, distB))) - distB

    ### SLICE INTO CENSORED AND UNCENSORED SUBSETS
    uncensored_vA = vA[~censA]
    censored_vA = vA[censA]
    uncensored_vB = vB[~censB]
    censored_vB = vB[censB]

    ### COMBINE UNCENSORED DISTRIBUTIONS
    uncensored_vAB = np.concatenate((uncensored_vA, uncensored_vB))
    xunc = np.unique(uncensored_vAB)

    ### DEFINE RANK ARRAYS
    nA = np.zeros(len(xunc), dtype='int')
    dA = np.zeros(len(xunc), dtype='int')
    mA = np.zeros(len(xunc), dtype='int')
    nB = np.zeros(len(xunc), dtype='int')
    dB = np.zeros(len(xunc), dtype='int')
    mB = np.zeros(len(xunc), dtype='int')

    ### COMPUTE RANKS FOR EACH DISTRIBUTION
    for i in np.arange(len(xunc)):
        xi = xunc[i]
        if (i < len(xunc)-1):
            xi1 = xunc[i+1]
        else:
            xi1 = xunc[len(xunc)-1]+1
        nA[i] = np.sum(vA >= xi)
        dA[i] = np.sum(uncensored_vA == xi)
        mA[i] = np.sum((censored_vA > xi) & (censored_vA < xi1))
        nB[i] = np.sum(vB >= xi)
        dB[i] = np.sum(uncensored_vB == xi)
        mB[i] = np.sum((censored_vB > xi) & (censored_vB < xi1))

    ### SUMMED RANKS
    nAB = nA+nB
    dAB = dA+dB
    mAB = mA+mB

    ### ESTIMATOR AND VARIANCE
    s = dA - 1. * dAB * nA / nAB
    var = 1. * dAB * nA * nB * (nAB - dAB) / (nAB[nAB > 2] - 1.)
    var[nAB == 1] = 0

    ### GEHAN
    L_gehan = np.sum(nAB * s)
    sig_gehan = np.sqrt(np.sum(var))
    P_gehan = 2 * (1 - norm.cdf(np.abs(L_gehan) / sig_gehan))

    ### LOGRANK
    L_logrank = np.sum(s)
    sig_logrank = np.sqrt(np.sum(var / nAB**2))
    P_logrank = 2 * (1 - norm.cdf(np.abs(L_logrank) / sig_logrank))

    ### PETO-PRENTICE
    w = np.zeros(len(xunc))
    b = np.zeros(len(xunc))
    for i in np.arange(len(w)):
        if (i == 0):
            w[i] = 1. * nAB[i] / (nAB[i] + 1)
            b[i] = 1. * (nAB[i] + 1) / (nAB[i] + 2)
        else:
            w[i] = w[i-1] * 1. * nAB[i] / (nAB[i] + 1)
            b[i] = b[i-1] * 1. * (nAB[i] + 1) / (nAB[i] + 2)
    g = 2 * dB + mB
    L_pp = np.sum(w * s)
    var = np.zeros_like(xunc)
    k = len(xunc) - 1
    for i in np.arange(k - 1):
        var[i] = w[i] * (1 - b[i]) * g[i] - (b[i] - w[i]) * g[i] * \
                 (w[i] * g[i] + 2 * np.sum(w[i+1:] * g[i+1:]))
    var[k] = w[k] * (1 - b[k]) * g[k] - w[k] * (b[k] - w[k]) * g[k]**2
    sig_pp = np.sqrt(np.sum(var))
    P_petoprentice = 2 * (1. - norm.cdf(np.abs(L_pp) / sig_pp))
    
    ### RETURN PROBABILITIES
    output = P_gehan, P_logrank, P_petoprentice

    return output


def do_mc(ref, comp, n_comp):

    """
    PURPOSE:    Perform Monte Carlo (MC) runs of two-sample tests
                following methodology of Andrews+2013 (2013ApJ...771..129A)

    INPUT:      ref = reference sample (table from get_lup())
                comp = comparison sample (table from get_tau() or get_usc())
                n_comp = comparison sample name (str)

    OUTPUT:     "twosample_test_[region].txt" = file with results of MC runs 
                of 3 types of two-sample tests calculted in censored_twosample()
    
    """

    ### SET VARIABLES
    nspt_comp  = np.array(comp['n_SpT'])
    nspt_ref   = np.array(ref['n_SpT'])
    Lmm_ref    = np.array(ref['FCont'])
    Lmm_comp   = np.array(comp['n_FCont'])
    flags_ref  = np.array(ref['det'] == 0)
    flags_comp = np.array(comp['det'] == 0)

    ### SET SPT BINS IN WHICH TO DISCRETIZE THE SAMPLES
    ### EVERY BIN MUST HAVE AT LEAST 1 SOURCE IN THE REFERENCE SAMPLE 
    blo = np.array([70, 65, 63, 62, 61, 58, 56, 53, 51, 45])
    bhi = np.array([65, 63, 62, 61, 58, 56, 53, 51, 45, 10])

    ### COMPUTE HOW MANY SOURCES ARE IN EACH BIN FOR THE COMPARISON SAMPLE
    nbin_comp = np.zeros(len(blo), dtype='int')
    for i in np.arange(len(blo)):
        nbin_comp[i] = np.sum( (nspt_comp <= blo[i]) & (nspt_comp > bhi[i]) )

    ### OPEN OUTPUT FILE
    out_file =  ('/').join(os.getcwd().split('/')[0:-1]) + '/output/twosample_test_' + n_comp + '.txt'
    os.system('rm -rf ' + out_file)

    ### PERFORM MC SAMPLING
    nmc = 10000
    print('\n Starting Monte Carlo sampling for ' + n_comp + '...')
    for im in np.arange(nmc):

        ### IN EACH SPT BIN, RANDOMLY DRAW (W/REPLACEMENT)
        ### A SUBSAMPLE FROM THE REFERENCE SAMPLE
        r_L = np.array([],dtype=int)
        r_f = np.array([],dtype=bool)
        for i in np.arange(len(blo)):
            
            if (nbin_comp[i] > 0):

                ### ALL REF FLUXES IN THIS SPT BIN (FLAG ND)
                t_L = Lmm_ref[(nspt_ref <= blo[i]) & (nspt_ref > bhi[i])]
                t_f = flags_ref[(nspt_ref <= blo[i]) & (nspt_ref > bhi[i])]

                ### RANDOMLY DRAW FROM REF SAMPLE WITH REPLACEMENT
                ### WHERE NUM OF DRAWS == NUM OF COMP SOURCES IN THIS SPT BIN
                rind = np.random.random_integers(0, len(t_L)-1, size=nbin_comp[i])

                ### STORE DRAWS IN ARRAYS                
                r_L = np.append(r_L, t_L[rind])
                r_f = np.append(r_f, t_f[rind])
    
        ### PRINT PROGRESS TO TERMINAL
        if (im % 1000 == 0 and im != 0): print("    {0} / {1} iterations".format(im,nmc))

        ### CALCULATE TWO-SAMPLE TESTS (p_null for Gehan, Logrank, Peto-Prentice)
        pg, pl, pp = censored_twosample(r_L, r_f, Lmm_comp, flags_comp)
        f = open(out_file, 'a')
        f.write("{0:e}   {1:e}   {2:e}\n".format(pg, pl, pp))
        f.close()


# ============================== Code ================================

#### LOAD IN TABLES 
TL = Table.read('../output/data_lup.txt', format='ascii.ipac')
TT = Table.read('../output/data_tau.txt', format='ascii.ipac')
TU = Table.read('../output/data_usc.txt', format='ascii.ipac')

### GET NUMERICAL SPECTRAL TYPE
TL['n_SpT'] = get_nspt(TL['SpT'])
TT['n_SpT'] = get_nspt(TT['SpT'])
TU['n_SpT'] = get_nspt(TU['SpT'])

### SET TAURUS FREQUENCIES
tt_freq = np.repeat(225.5, len(TT))
tt_freq[np.where(TT['ObsType'] == '890')] = 338.

### NORMALIZE
TL['n_FCont'] = norm_flux(TL['FCont'], TL['Dis'], np.repeat(335.8, len(TL)))
TU['n_FCont'] = norm_flux(TU['FCont'], np.repeat(145., len(TU)), np.repeat(341.1, len(TU)))
TT['n_FCont'] = norm_flux(TT['FCont'], np.repeat(140., len(TT)), tt_freq)

### PERFORM MC SAMPLING
do_mc(TL, TU,'usc')
do_mc(TL, TT,'tau')

