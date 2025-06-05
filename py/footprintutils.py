from astropy.io import fits
import numpy as np

def fake_header(crval1, crval2, crpix=512, crpix2=2044,crpix1=2044, cdelt1=0.11, cdelt2=0.11,
                crota2=0.0,naxis1=4088,naxis2=4088):

    # crota2 - degree
    # cdelt1 - arcsec
    # cdelt2 - arcsec

    hdu = fits.PrimaryHDU()
    hdu.header

    # http://stsdas.stsci.edu/documents/SUG/UG_21.html

    theta = crota2*np.pi/180. # radians
    cdelt1 /= 3600. # deg
    cdelt2 /= 3600. # deg

    #cd1_1 = -cdelt1*np.cos(theta)
    #cd1_2 = -cdelt2*np.sin(theta)
    #cd2_1 = -cdelt1*np.sin(theta)
    #cd2_2 = cdelt2*np.cos(theta)
    
    R = np.array([
        [-1*np.cos(theta), -1*np.sin(theta)],
        [-1*np.sin(theta), np.cos(theta)],
    ])


    cd1_1 = cdelt1*R[0,0]
    cd1_2 = cdelt2*R[0,1]
    cd2_1 = cdelt1*R[1,0]
    cd2_2 = cdelt2*R[1,1]
                    
    hdu.header.set('NAXIS',2) # pixels
    hdu.header.set('NAXIS1',naxis1) # pixels
    hdu.header.set('NAXIS2',naxis2) # pixels

    hdu.header.set('WCSAXES',2) # pixels
    hdu.header.set('CTYPE1','RA---TAN')
    hdu.header.set('CTYPE2','DEC--TAN')
    hdu.header.set('CRVAL1',crval1) # deg
    hdu.header.set('CRVAL2',crval2) # deg
    hdu.header.set('CRPIX1',crpix1) # pixels
    hdu.header.set('CRPIX2',crpix2) # pixels
    #hdu.header.set('CDELT1',cdelt1)
    #hdu.header.set('CDELT2',cdelt2)
    #hdu.header.set('CROTA2',crota2)
    hdu.header.set("CD1_1",cd1_1)
    hdu.header.set("CD1_2",cd1_2)
    hdu.header.set("CD2_1",cd2_1)
    hdu.header.set("CD2_2",cd2_2)

    return hdu.header
