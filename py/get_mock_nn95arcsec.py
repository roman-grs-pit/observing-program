import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import fitsio
from matplotlib import pyplot as plt
import h5py
import numpy as np

indir = '/global/cscratch1/sd/ajross/RomanGRS/'

f = h5py.File(indir+'ATLAS_small_V2_0.hdf5')
data = f['data']
data = np.array(data)
print(len(data))
sel = data.transpose()[8]*.28685 > 1e-16
print(len(data[sel]))
data = data[sel]
for i in range(1,10):
    f = h5py.File(indir+'ATLAS_small_V2_'+str(i)+'.hdf5')
    dn = f['data']
    dn = np.array(dn)
    sel = dn.transpose()[8]*.28685 > 1e-16
    data = np.concatenate((data,dn[sel]))
print(len(data))

data = data.transpose()

cd = SkyCoord(ra=data[0]*u.degree, dec=data[1]*u.degree)

idxcn_all = []
idxcd_all = []
n_all = []
for i in range(0,10):
    f = h5py.File(indir+'ATLAS_small_V2_'+str(i)+'.hdf5')
    dn = f['data']
    dn = np.array(dn)
    sel = dn.transpose()[8]*.28685 > 1e-16
    dn = dn[sel]
    dn = dn.transpose()
    n_all.append(len(dn[0]))
    print(n_all)
    cn = SkyCoord(ra=dn[0]*u.degree, dec=dn[1]*u.degree)
    idxcn, idxcd, d2d, d3d = cn.search_around_sky(cd,93.0 * u.arcsec)
    idxcn_all.append(idxcn)
    idxcd_all.append(idxcd)
    
    print(i)

idxcatalog = np.concatenate(idxcn_all)

idx,nc = np.unique(idxcatalog,return_counts=True)

np.savetxt(indir+'num_neighbor_93arcsec.txt',nc)

del idxcatalog
del idx
del idxcn_all
del idxcd_all
del cd

#now add columns and sub-sampling

tab = Table(data)
names = ['RA','DEC','Z_COS','Z_OBS','VEL','STARMASS','SFR','HALOMASS','flux_Halpha','flux_OIII','Habs','CENTRAL']
for i in range(0,len(names)):
    tab.rename_column('col'+str(i), names[i])

n93 = np.loadtxt(indir+'num_neighbor_93arcsec.txt')
tab['NUM93arcsec'] = nc

subfrac = nc/500 #the distribution peaks at ~50, we sub-select 90% at 50 and alter linearly with nc

ran = np.random.random(len(subfrac))

sel1 = ran > subfrac

tab['sub10per'] = sel1

sel2 = ran > 2*subfrac

tab['sub20per'] = sel2

tab.write(indir+'ATLAS_small_concat_Ha16_wsub.fits',format='fits')


