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

import logging

# create logger
logname = 'Roman_coverage'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)



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

from time import time


def get_pixl_siaf(ra,dec,att_in,detnum):
    t0 = time()
    rap = f'WFI{detnum :02}_FULL'
    wfi = rsiaf[rap]
    wfi.set_attitude_matrix(att_in)
    cen_ra,cen_dec = wfi.idl_to_sky(0, 0)
    t1 = time()
    #print(str(t1-t0)+ 'setup')
    #ddec = abs(dec-cen_dec)
    #dra = abs(ra-cen_ra)
    ddec = dec-cen_dec
    dra = ra-cen_ra
    sel = ddec > -0.1#ddec < 0.1
    sel &= ddec < 0.1
    dfac = np.cos(np.min(dec)*np.pi/180)
    sel &= dra < 0.1/dfac
    sel &= dra > -0.1/dfac
    t2 = time()
    #print(str(t2-t1)+' masked array')
    #pixels = np.copy(in_array)
    #pixels[0][sel] = -999
    t3 = time()
    #print(str(t3-t2)+' initialize array')
    
    
    pixels_sel = wfi.sky_to_sci(ra[sel],dec[sel])
    #pixels[0][sel] = pixels_sel[0]
    #pixels[1][sel] = pixels_sel[1]
    #pixels[2] = sel
    t4 = time()
    #print(str(t4-t3)+' final result')
    return pixels_sel,sel#pixels

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
parser.add_argument("--wficen", help="if y, positions are detector center",default='y')
parser.add_argument("--wavmin", help="set minimum wavelength, if not None",default=None)
parser.add_argument("--wavmax", help="set maximum wavelength, if not None",default=None)


parser.add_argument("--randens", help="density of random points /deg2 to use",default=25,type=float)
parser.add_argument("--Nchunk", help="how many times to repeat",default=10,type=int)
parser.add_argument("--set", help="which set of chunks?",default=0,type=int)
parser.add_argument("--tiles", help="which set of tiles?",default='sd')



args = parser.parse_args()



decm = -60
decx = 10
ram = 0
rax = 360
randens = args.randens
fullsky_area = 360.*360/np.pi
cosmin = np.cos(np.pi/180*(decm+90))*-1
cosmax = np.cos(np.pi/180*(decx+90))*-1
fracdec = (cosmax-cosmin)/2
nran = int(randens*fullsky_area*fracdec)
print(cosmin,cosmax,nran,nran/fracdec,fracdec)

minwav = 1
maxwav = 1.9
wavstr = ''
if args.wavmin is not None:
    minwav = float(args.wavmin)
    wavstr += 'lam'+str(round(minwav,3))
if args.wavmax is not None:
    maxwav = float(args.wavmax)
    wavstr += str(round(maxwav,3))
wfistr =''
if args.wficen != 'y':
    wfistr = 'notwficen'    
outdir = out_root+'fullsurvey_'+args.tiles+'/ramin'+str(ram)+'decmin'+str(decm)+wavstr+wfistr+'/'
logger.info('results will be written to '+outdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)

if args.tiles == 'sd':
    tiles = np.loadtxt(os.environ['github_dir']+'observing-program/data/hlwas_tiling_241206.txt').transpose()
    gtiles = tiles[4] == 9
 
    racol = 2
    deccol = 1
    pacol = 3
    pad = 0.2
    if args.wficen != 'y':
        pad = 0.5
    print(len(tiles[0][gtiles]))
    tls = tiles[0][gtiles]
    
if args.tiles == 'socv0':
    tiles = Table.read(os.environ['github_dir']+'observing-program/data/tillingfromJavi_994_fixed_workaround.sim.ecsv')
    gtiles = tiles['BANDPASS'] == 'GRISM'
 
    racol = 'RA'
    deccol = 'DEC'
    pacol = 'PA'
    pad = 0.2
    if args.wficen != 'y':
        pad = 0.5
    tls = tiles[gtiles]
    
dets = np.arange(1,19)

ra_all = []
dec_all = []
cnts_all = []
indx_all = []

tottl = len(tls)

min_chunk = args.set*args.Nchunk
for chunk in range(int(min_chunk),args.Nchunk):
    rand_indx = []


    #selcos = cosl > cosmin
    #selcos &= cosl < cosmax
    cosl = np.random.rand(nran)*2*fracdec+cosmin#2-1 #distribute randomly in arccos
    
    ral = np.random.rand(nran)*360-180
    ral = ral
    decl = np.arccos(-1*cosl)*180/np.pi-90
    
    logger.info('made all sky randoms and cut to declination range')
    
    ral_tot,decl_tot = cutran(ram,rax,decm,decx)#mkgrid(0,0,1,100)
    logger.info(str(len(ral_tot))+' random points will be used')
    ran_indices = (np.arange(len(ral_tot))+ args.Nchunk*1e7).astype(int)
    in_array = np.ones((3,len(ral_tot)))*-9999
    coords = SkyCoord(ra=ral_tot*u.degree,dec=decl_tot*u.degree, frame='icrs')
        
    ral_tot = np.array(ral_tot)
    decl_tot = np.array(decl_tot)
    def get_idx_tl(tl):
        ra0 = tiles[racol][gtiles][tl]
        if ra0 > 180:
            ra0 -= 360
        dec0 = tiles[deccol][gtiles][tl]
        pa = tiles[pacol][gtiles][tl]
        if args.wficen == 'y':
            att = attitude(wfi_cen.V2Ref, wfi_cen.V3Ref, ra0, dec0, pa)
        else:
            att = attitude(0, 0, ra0, dec0, pa)
        idx = []
        for det in dets:
            #pixels = get_pixl(coords,dfoot,det,PA-pa_off)
            #pixels = get_pixl_siaf(np.array(ral_tot),np.array(decl_tot),att,det)
            pixel_sel,sel = get_pixl_siaf(ral_tot,decl_tot,att,det)
            selp = sel.astype(bool)#pixels[2].astype(bool)
            #print(np.sum(selp),len(selp))
            #for i in range(0,len(pixels[0][selp])):
            for i in range(0,len(pixel_sel[0])):
                #xpix = pixels[0][selp][i]
                #ypix = pixels[1][selp][i]
                xpix = pixel_sel[0][i]
                ypix = pixel_sel[1][i]
    
                test = 0
                if xpix > -1000 and xpix < 5088 and ypix > -1000 and ypix < 5088:
                    test = test_foot(xpix,ypix,det=det,min_lam_4foot=minwav,max_lam_4foot=maxwav)
                    if test == 1:
                        idx_det = ran_indices[selp][i]
                        idx.append(idx_det)
            #logger.info('completed detector '+str(det)+' on obs '+str(tl))
        #logger.info('completed '+str(tl)+' out of '+str(tottl))
        return idx
    
    par = 'y'
    if par == 'n':
        for tl in range(0,len(tls)):
            idx = get_idx_tl(tl)
            rand_indx.append(idx)
            print(str(tl)+' completed')
    
    if par == 'y':
        from concurrent.futures import ProcessPoolExecutor
        tl_idx = list(np.arange(len(tls)).astype(int))   
        with ProcessPoolExecutor() as executor:
            for idx in executor.map(get_idx_tl, tl_idx):
                rand_indx.append(idx)

    logger.info('completed chunk '+str(chunk))
    logger.info('length of list of random ids '+str(len(rand_indx)))
    rand_indx = np.concatenate(rand_indx)
    logger.info('length of concatenated array of random ids '+str(len(rand_indx)))
    rans,cnts = np.unique(rand_indx,return_counts=True)
    selobs = np.isin(ran_indices,rans)
    if np.array_equal(ran_indices[selobs],rans):
        logger.info('input/final ids are matched in order')
    else:
        sys.exit('ids are not matched')

    racut = ral_tot[selobs]
    deccut = decl_tot[selobs]
    ra_all.append(racut)
    dec_all.append(deccut)
    cnts_all.append(cnts)
    indx_all.append(rans)

    
ra_all = np.concatenate(ra_all)
dec_all = np.concatenate(dec_all)
cnts_all = np.concatenate(cnts_all)
indx_all = np.concatenate(indx_all)



tout = Table()
tout['RA'] = ra_all
tout['DEC'] = dec_all
tout['ID'] = indx_all
tout['NOBS'] = np.array(cnts_all,dtype=int)
tout.write(outdir+'nobs'+str(minwav)+str(maxwav)+'_chunkranset'+str(set)+'.ecsv',overwrite=True)

#make nobs figure
plt.clf()
cmap = plt.get_cmap('jet', 8)
#sel = nobs > 0
plt.scatter(np.array(ral_tot)[selobs],np.array(decl_tot)[selobs],c=cnts,s=.1,cmap=cmap,vmin=0.5,vmax=8.5)
plt.colorbar(ticks=np.arange(1,9 ),label='# of obs.')
#plt.colorbar(levels=[0,1,2,3,4])
plt.xlim(ram, rax)
plt.ylim(decm,decx)
plt.xlabel('ra (degrees)')
plt.ylabel('dec (degrees)')
#plt.title(r'footprint 4 rolls/2 dithers')
plt.savefig(outdir+'nobs'+str(minwav)+str(maxwav)+'.png')
plt.show()
