import os
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import importlib
import argparse

#on NERSC; should change this so that it is already sourced by environment script
sys.path.append('/global/common/software/m4943/grizli0/observing-program/py/')
sys.path.append('/global/common/software/m4943/grizli0/grism_sim/py/')
os.environ["WEBBPSF_PATH"]="/global/cfs/cdirs/m4943/grismsim/webbpsf-data"
out_root = '/global/cfs/cdirs/m4943/footprint/'
os.environ['github_dir']='/global/common/software/m4943/grizli0/'

#similarly, should change to make this an environment variable
code_data_dir = '/global/common/software/m4943/grizli0/observing-program/data/'

#import roman_coords_transform as ctrans
#rctrans = ctrans.RomanCoordsTransform(file_path=code_data_dir)

import pysiaf
from pysiaf.utils.rotations import attitude
rsiaf = pysiaf.Siaf('Roman')
wfi_cen = rsiaf['WFI_CEN']

import footprintutils as fp
import optical_model
optmod = optical_model.RomanOpticalModel()



#should move this to code library
def test_foot(xpix,ypix,min_pix=0,max_pix=4088,det=1,min_lam_4foot=1.6,max_lam_4foot=1.93):
    #xpix,ypix is the pixel center for the object in the direct image
    #min_pix,max_pix represent the detector bounds in pixels
    #det is the SCA detector number
    #min_lam_4foot is the minimum wavelength required to be considered
    #max_lam_4foot is the maximum wavelength required to be considered
    test = optmod._get_beam_trace(xpix,ypix,det,width=1)
    min_wv_ind = int((min_lam_4foot-.9)/0.001)#trace has a minimum wavelength of 0.9 microns and spacing of 0.001
    if min_wv_ind < 0:
        min_wv_ind = 0
    max_wv_ind = int((max_lam_4foot-.9)/0.001)+1
    lt = len(test['trace_pix_x'][0])
    if max_wv_ind > lt:
        max_wv_ind = lt
    if np.min(test['trace_pix_x'][0][min_wv_ind:max_wv_ind])>=min_pix and np.max(test['trace_pix_x'][0][min_wv_ind:max_wv_ind]) < max_pix and np.min(test['trace_pix_y'][0][min_wv_ind:max_wv_ind])>=min_pix and np.max(test['trace_pix_y'][0][min_wv_ind:max_wv_ind]) < max_pix:
        return 1
    else:
        return 0

#just make a square grid of ra,dec points; should be ok for small area; large area should use properly produced randoms (or similar)
def mkgrid(ra,dec,sz,res):    
    ral = np.arange(ra-sz/2,ra+sz/2,sz/res)
    decl = np.arange(dec-sz/2,dec+sz/2,sz/res)
    if len(ral) != len(decl):
        print('error, mismatched ra/dec lists!')
        return
    ral_tot  = []
    decl_tot = []
    for i in range(0,len(ral)):
        for j in range(0,len(ral)):
            ral_tot.append(ral[i])
            decl_tot.append(decl[j])
    return ral_tot,decl_tot

def cutran(ramin,ramax,decmin,decmax):
    sel = ral > ramin
    sel &= ral < ramax
    sel &= decl > decmin
    sel &= decl < decmax
    return ral[sel],decl[sel]

def get_pixl(coords,detfoot,detnum,PA):
    #coords are ra,dec points to test in astropy format
    #detfoot is the object containing the footprint of the SCA detectors, centers get used below
    #detnum is the SCA detector number
    #PA is the rotation relative to lines of constant ra; note that detfoot is rotated 60 deg. w.r.t. this
    ra = detfoot[0][int(detnum)]['ra_cen']
    dec = detfoot[0][int(detnum)]['dec_cen']
    h = fp.fake_header(ra, dec,crota2=PA)
    w = WCS(h)
    pixels = w.world_to_pixel(coords)
    return pixels

def get_pixl_siaf(ra,dec,att_in,detnum):
    rap = f'WFI{detnum :02}_FULL'
    wfi = rsiaf[rap]
    wfi.set_attitude_matrix(att_in)
    pixels = wfi.sky_to_sci(ra,dec)
    cen_ra,cen_dec = wfi.idl_to_sky(0, 0)
    ddec = abs(dec-cen_dec)
    sel = ddec > 0.1
    pixels[0][sel] *= -999
    return pixels

def plot_dets_rsiaf(att_in,ax):
    #roman_apertures = [f'WFI{i + 1:02}_FULL' ]
      
    #for rap in roman_apertures:
    for i in range(18):
        rap = f'WFI{i + 1:02}_FULL'
        wfi = rsiaf[rap]
        wfi.set_attitude_matrix(att_in)

        wfi_ra, wfi_dec = wfi.idl_to_sky(0, 0)
        #plt.plot(wfi_ra,wfi_dec,'k,')
        ax.text(wfi_ra,wfi_dec,str(i+1))
        corners_x = np.array([1,4088,4088,1,1])
        corners_y = np.array([1,1,4088,4088,1])
        corners_ra,corners_dec = wfi.sci_to_sky(corners_x, corners_y)
        ax.plot(corners_ra,corners_dec,'k')


parser = argparse.ArgumentParser()
parser.add_argument("--ra0", help="ra center",default=2)
parser.add_argument("--randens", help="density of random points /deg2 to use",default=2500)
parser.add_argument("--dec0", help="dec center",default=0)
parser.add_argument("--dithra", help="dither step in the ra direction",default=0.025)
parser.add_argument("--dithdec", help="dither step in the dec direction",default=0.05)
parser.add_argument("--pa_step_dec", help="step in the dec direction to take when flipping ~180deg",default=0.087)
parser.add_argument("--pa_step_ra", help="step in the ra direction with each pa angle",default=0)


args = parser.parse_args()

ra0 = args.ra0
dec0 = args.dec0

decm = dec0-0.6
decx = dec0+0.6
ram = ra0-0.4
rax = ra0+0.4
randens = args.randens
fullsky_area = 360.*360/np.pi
nran = randens*fullsky_area
acosl = np.random.rand(nran)*2-1 #distribute randomly in arccos
ral = np.random.rand(nran)*360-180
decl = np.arccos(acosl)

ral_tot,decl_tot = cutran(ram,rax,decm,decx)#mkgrid(0,0,1,100)
print(str(len(ral_tot))+' random points will be used')
coords = SkyCoord(ra=ral_tot*u.degree,dec=decl_tot*u.degree, frame='icrs')


dithstep = args.dithra,args.dithdec
ndith = 2
decpa = args.pa_step_dec
rapa = args.pa_step_ra
outdir = out_root+'ra'+str(ra0)+'dec'+str(dec0)+'_dithra'+str(args.dithra)+'dec'+str(args.dithdec)+'_para'+str(args.pa_step_ra)+str(args.pa_step_dec)+'/'
print('results will be written to '+outdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)
decoffl = [decpa,decpa,-decpa,-decpa]
raoffl = [rapa,rapa/3,-rapa/3,-rapa]
pa_off = 60 #this undoes the rotation applied by default to a pa=0 
pal = [-2.5+pa_off,2.5+pa_off,177.5+pa_off,182.5+pa_off]

fig = plt.figure()
ax = fig.add_subplot(111)
pa = pal[0]
decoff = decoffl[0]#decpa/2
raoff = raoffl[0]
#if pa > 100:
#    decoff = -decpa/2
#a = rctrans.wfi_sky_pointing(ra0+raoff, dec0+decoff, pa, ds9=False,ax=ax)#+dith*dithstep[0]
att = attitude(wfi_cen.V2Ref, wfi_cen.V3Ref, ra0+raoff, dec0+decoff, pa)
plot_dets_rsiaf(att,ax)
ax.set_xlim(ram, rax)
ax.set_ylim(decm,decx)
plt.grid()
plt.title('1 exposure')
plt.savefig(outdir+'oneexp.png')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
pa = pal[0]
for dith in range(0,ndith):
    decoff = decoffl[0]#decpa/2
    raoff = raoffl[0]
    #if pa > 100:
    #    decoff = -decpa/2
    #a = rctrans.wfi_sky_pointing(ra0+raoff+dith*dithstep[0], dec0+decoff+dith*dithstep[1], pa, ds9=False,ax=ax)#+dith*dithstep[0]
    att = attitude(wfi_cen.V2Ref, wfi_cen.V3Ref, ra0+raoff+dith*dithstep[0], dec0+decoff+dith*dithstep[1], pa)
    plot_dets_rsiaf(att,ax)
ax.set_xlim(ram, rax)
ax.set_ylim(decm,decx)
plt.title('2 dithers')
plt.savefig(outdir+'twodithers.png')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
for pa,decoff,raoff in zip(pal,decoffl,raoffl):
    #decoff = decpa/2
    #if pa > 100:
    #    decoff = -decpa/2
    #a = rctrans.wfi_sky_pointing(ra0+raoff, dec0+decoff, pa, ds9=False,ax=ax)#+dith*dithstep[0]
    att = attitude(wfi_cen.V2Ref, wfi_cen.V3Ref, ra0+raoff, dec0+decoff, pa)
    plot_dets_rsiaf(att,ax)

ax.set_xlim(ram, rax)
ax.set_ylim(decm,decx)
plt.title('4 rolls')
plt.grid()
plt.savefig(outdir+'fourrolls.png')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
for pa,decoff,raoff in zip(pal,decoffl,raoffl):
    for dith in range(0,ndith):
        #decoff = decpa/2
        #if pa > 100:
        #    decoff = -decpa/2
        #a = rctrans.wfi_sky_pointing(ra0+raoff+dith*dithstep[0], dec0+decoff+dith*dithstep[1], pa, ds9=False,ax=ax)#+dith*dithstep[0]
        att = attitude(wfi_cen.V2Ref, wfi_cen.V3Ref, ra0+raoff+dith*dithstep[0], dec0+decoff+dith*dithstep[1], pa)
        plot_dets_rsiaf(att,ax)

ax.set_xlim(ram, rax)
ax.set_ylim(decm,decx)
plt.title('all 8 exposures')
plt.savefig(outdir+'all8exp.png')
plt.show()

minwav = 1
maxwav = 1.9
nobs = np.zeros(len(ral_tot))
#minwav = (1+ztest)*halpha
#minwav = (1+ztest)*hbeta
dets = np.arange(1,19)
for PA,decoff,raoff in zip(pal,decoffl,raoffl):
    for dith in range(0,ndith):
        #decoff = decpa/2
        #if PA > 100:
        #    decoff = -decpa/2

        #dfoot = rctrans.wfi_sky_pointing(ra0+raoff+dith*dithstep[0], dec0+decoff+dith*dithstep[1], PA, ds9=False,do_plot=False)
        att = attitude(wfi_cen.V2Ref, wfi_cen.V3Ref, ra0+raoff+dith*dithstep[0], dec0+decoff+dith*dithstep[1], PA)
        for det in dets:
            #pixels = get_pixl(coords,dfoot,det,PA-pa_off)
            pixels = get_pixl_siaf(np.array(ral_tot),np.array(decl_tot),att,det)
            for i in range(0,len(pixels[0])):
                xpix = pixels[0][i]
                ypix = pixels[1][i]
                test = 0
                if xpix > -1000 and xpix < 5088 and ypix > -1000 and ypix < 5088:
                    test = test_foot(xpix,ypix,det=det,min_lam_4foot=minwav,max_lam_4foot=maxwav)
                nobs[i] += test

tout = Table()
tout['RA'] = ral_tot
tout['DEC'] = decl_tot
tout['NOBS'] = np.array(nobs,dtype=int)
tout.write(outdir+'nobs'+str(minwav)+str(maxwav)+'grid.ecsv',overwrite=True)

#make nobs figure
plt.clf()
cmap = plt.get_cmap('jet', 8)
sel = nobs > 0
plt.scatter(np.array(ral_tot)[sel],np.array(decl_tot)[sel],c=nobs[sel],s=.1,cmap=cmap,vmin=0.5,vmax=8.5)
plt.colorbar(ticks=np.arange(1,9 ),label='# of obs.')
#plt.colorbar(levels=[0,1,2,3,4])
plt.xlim(ram, rax)
plt.ylim(decm,decx)
plt.xlabel('ra (degrees)')
plt.ylabel('dec (degrees)')
plt.title(r'footprint 4 rolls/2 dithers')
plt.savefig(outdir+'nobs'+str(minwav)+str(maxwav)+'.png')
plt.show()
