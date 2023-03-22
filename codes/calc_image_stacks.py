# ======================== Import Packages ==========================

import sys, os, pdb, glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
from astropy import units as u


# ===================== Define Functions ===================

def readfits(file):

    """
    PURPOSE:    Read in FITS file and header info

    INPUT:      Path to FITS file (str)

    OUTPUT:     Image (float arr)
                Image center coordinates in pixels (float)
                Image pixel width in deg/pix units (float)
                Beam major axis, minor axis, position angle (float)
                Image center coordinates in deg units (float)

    """

    ### READ IN FITS FILE
    hdulist = fits.open(file)
    data = hdulist[0].data[0, 0, :, :]
    head = hdulist[0].header
    hdulist.close()

    ### GET HEADER INFO
    xcen = head['CRPIX1']
    ycen = head['CRPIX2']
    xpix = head['CDELT1']
    ypix = head['CDELT2']
    xcen_ra = head['CRVAL1']
    xcen_de = head['CRVAL2']
    bmaj = head['bmaj']
    bmin = head['bmin']
    bpa  = head['bpa']

    return(data, xcen, ycen, xpix, ypix, bmaj, bmin, bpa, xcen_ra, xcen_de)


def write_fits(img, line, i, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all):

    """
    PURPOSE:    Write FITS file and header info

    INPUT:      Image to write (array)
                Line that was stacked (str)
                Pixel width in deg/pix units for all images in stack (array)
                Beam major axis, minor axis, position angle for all images in stack (array)

    OUTPUT:     Stacked image (FITS file)

    """

    os.system('rm ../output/stack_nd_'+line+'_'+str(i)+'.fits')
    hdu = fits.PrimaryHDU()
    hdu.data = img

    hdu.header['CRPIX1'] = hdu.header['NAXIS1']/2
    hdu.header['CRPIX2'] = hdu.header['NAXIS2']/2
    hdu.header['bmaj'] = bmaj_all.mean()
    hdu.header['bmin'] = bmin_all.mean()
    hdu.header['bpa'] = bpa_all.mean()
    hdu.header['cdelt1'] = xpix_all.mean()
    hdu.header['cdelt2'] = ypix_all.mean()

    hdu.writeto('../output/stack_nd_'+line+'_'+str(i)+'.fits')


def crop_img(file_img, hw_as, c_obj):

    ### LOAD IMAGE AND GET CENTER COORDINATES
    img, xcen_img, ycen_img, xpix_img, ypix_img, bmaj_img, bmin_img, bpa_img, xcen_ra_img, ycen_de_img = readfits(file_img)
    c_img = SkyCoord(xcen_ra_img, ycen_de_img, frame='icrs', unit='deg')
    
    ### CENTER IMAGE ON OBJECT LOCATION 
    dra, ddec = c_img.spherical_offsets_to(c_obj)
    width_pix = int(round(hw_as / (ypix_img * 3600.0)))
    xctr = xcen_img + dra.value / xpix_img
    yctr = ycen_img + ddec.value / ypix_img

    ### CROP IMAGE
    img = img[int(round(yctr - width_pix)):int(round(yctr + width_pix)),
              int(round(xctr - width_pix)):int(round(xctr + width_pix))]
    
    return img, width_pix, xpix_img, ypix_img, bmaj_img, bmin_img, bpa_img


def stackme(t, line):

    """
    PURPOSE:    Stack image

    INPUT:      Table of sources to be stacked (AstroPy Table)
                Line name (str; must be 'cont', '13CO', C18O')

    OUTPUT:     Stacked image (array)
                Pixel width in deg/pix units for all images in stack (array)
                Beam major axis, minor axis, position angle for all images in stack (array)

    """

    xpix_all, ypix_all = np.empty(len(t)), np.empty(len(t))
    bmaj_all, bmin_all, bpa_all = np.empty(len(t)), np.empty(len(t)), np.empty(len(t))
    
    for i,val in enumerate(t['Name']):

        if (line == 'C18O'): suffix = '_c18o32.mom0.fits'
        if (line == '13CO'): suffix = '_13co32.mom0.fits'
        if (line == 'cont'):   suffix = '_cont.fits'

        file_img = '../data/' + val + suffix
        file_img = file_img.replace(' ', '_')
        if os.path.isfile(file_img) is False:
            print('missing FITS file for ' + val, line)
            pdb.set_trace()

        ### GET COORDINATES OF OBJECT FROM PAPER TABLE
        de_obj = str(t['DE-'][i]) + str(t['DEd'][i]) + 'd' + str(t['DEm'][i]) + 'm' + str(t['DEs'][i]) + 's'
        ra_obj = str(t['RAh'][i]) + 'h' + str(t['RAm'][i]) + 'm' + str(t['RAs'][i]) + 's'
        c_obj = SkyCoord(ra_obj, de_obj, frame='icrs')

        img_cont, width_pix, xpix_img, ypix_img, bmaj_img, bmin_img, bpa_img = crop_img(file_img, 8.0, c_obj)
        
        ### SCALE WITH DISTANCE AND PUT INTO MJY UNITS
        img_cont = 1e3 * img_cont * ((t['Dis'][i] / 200.)**2)


        xpix_all[i], ypix_all[i] = xpix_img, ypix_img
        bmaj_all[i], bmin_all[i], bpa_all[i] = bmaj_img, bmin_img, bpa_img

        if (i==0):
            img_all = np.zeros([2 * width_pix, 2 * width_pix, 1])
            temp = img_cont.reshape((2 * width_pix, 2 * width_pix, 1))
            img_all = temp
            
        else:
            temp = img_cont.reshape((2 * width_pix, 2 * width_pix, 1))
            img_all = np.append(img_all, temp, axis=2)

    stacked = np.sum(img_all, 2) / len(t)

    return stacked, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all


# ========================== Code ==========================

#### LOAD IN TABLES FROM PAPER SUPPLEMENTAL MATERIAL
TS = Table.read('../input/t1_mrf.txt', format='ascii.cds')
TD = Table.read('../input/t2_mrf.txt', format='ascii.cds')
TG = Table.read('../input/t3_mrf.txt', format='ascii.cds')
T = join(TS, TD, join_type='inner')
T = join(T, TG, join_type='inner')
    
### INDEX DUST NON-DETECTIONS
ind_dust_nd = T['FCont'] / T['e_FCont'] < 3.0

### INDEX dust detections, 13CO non-detections, C18O non-detections
ind_gas_nd = ((T['FCont'] / T['e_FCont'] > 3.0)  & (~T['l_F13CO'].mask) & (~T['l_F18CO'].mask) )

### INDEX dust detections, 13CO detections, 18CO non-detections
ind_C18O_nd = ((T['FCont'] / T['e_FCont'] > 3.0)  & (T['l_F13CO'].mask) & (~T['l_F18CO'].mask) )

### STACK IMAGES
lines = ['13CO', 'C18O']
for n, nval in enumerate(lines):
    for i, val in enumerate([T[ind_dust_nd], T[ind_gas_nd], T[ind_C18O_nd]]):
        stacked, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all = stackme(val, lines[n])
        write_fits(stacked, lines[n], i, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all)

plt.close('all')
