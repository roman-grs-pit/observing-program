import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
#import imageio
import os 
from scipy.interpolate import griddata

#center_long,center_lat = 60,-50 #angular coordinates of center of the field
center_long,center_lat = 30,-41 #angular coordinates of center of the field
sec_size = 3.32 #sector size in degrees

#define the number of sectors in each strip
#nstrips = np.array([2,3,4,5,6,7,8,9,10,11,12,12,11,10,9,8,7,6,5,4,3,3]) #this was the fiducial one

nstrips = np.array([9,10,11,12,13,14,14,13,12,11,10,9,9,8])
print(sum(nstrips))
strips = np.arange(len(nstrips))
#print(strips)

lat_extent = sec_size*len(strips)

min_lat = center_lat-lat_extent/2.
max_lat = min_lat+lat_extent

lat_edges = np.linspace(min_lat,max_lat,len(strips)+1)
lat_mid = []
for ii in range(0,len(lat_edges)-1):
    lat_mid.append((lat_edges[ii]+lat_edges[ii+1])/2.)

#print(lat_mid)

strip_edges_low ={} #longitude coordinate on the low end of the latitude
strip_edges_high ={} #longitude coordinate on the high end of the latitude

def radec2thphi(ra,dec):
    #for healpix
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.


def get_long_minmax(lat,ns):
    latfac = np.cos(lat*np.pi/180.)
    long_extent = ns*sec_size/latfac
    min_long = center_long-long_extent/2.
    offset = 0
    #if lat > -40:
    #    min_long += 40+lat
    if hasattr(lat, "__len__"):
        sel = lat > -40
        min_long[sel] += 40+lat[sel]
    else:
        if lat > -40:
            min_long += 40+lat

    max_long = min_long+long_extent
    return min_long,max_long


def get_long_edges(lat,ns):
    latfac = np.cos(lat*np.pi/180.)
    long_extent = ns*sec_size/latfac
    min_long = center_long-long_extent/2.
    
    if lat > -40:
        min_long += 40+lat

    max_long = min_long+long_extent
    long_edges = np.linspace(min_long,max_long,ns+1)
    return long_edges


for strip,ns,lat in zip(strips,nstrips,lat_mid):
    #print(strip)
    lat_low = lat-sec_size/2.
    long_edges = get_long_edges(lat_low,ns)
    strip_edges_low[strip] = long_edges
    
    lat_high = lat_low + sec_size
    long_edges = get_long_edges(lat_high,ns)
    strip_edges_high[strip] = long_edges

def get_corners(sn,nn):
    #sn is the strip number
    #nn is sector number in the strip
    ra0 = strip_edges_low[sn][nn]
    dec0 = lat_edges[sn]
    ra1 = strip_edges_low[sn][nn+1]
    dec1 = lat_edges[sn]
    ra2 = strip_edges_high[sn][nn+1]
    dec2 = lat_edges[sn+1]
    ra3 = strip_edges_high[sn][nn]
    dec3 = lat_edges[sn+1]
    ral = [ra0,ra1,ra2,ra3,ra0]
    decl = [dec0,dec1,dec2,dec3,dec0]
    return ral,decl

def get_centers(sn,nn):
    #sn is the strip number
    #nn is sector number in the strip
    ra0 = strip_edges_low[sn][nn]
    dec0 = lat_edges[sn]
    ra1 = strip_edges_low[sn][nn+1]
    dec1 = lat_edges[sn]
    ra2 = strip_edges_high[sn][nn+1]
    dec2 = lat_edges[sn+1]
    ra3 = strip_edges_high[sn][nn]
    dec3 = lat_edges[sn+1]
    ral = [(ra0+ra1+ra2+ra3)/4.]
    decl = [(dec0+dec1+dec2+dec3)/4.]
    return ral,decl

def write_table_4back(year=2027,day=10,microns=1.4,daystep=10,ido_view='0'):
    '''
    this produces a table to input to https://irsa.ipac.caltech.edu/applications/BackgroundModel/
    '''
    rall = []
    decll = []
    ntot = 0
    for sn in strips:
        for nn in range(0,nstrips[sn]):
            ral,decl = get_centers(sn,nn)
            rall.append(ral)
            decll.append(decl)
            #plt.plot(ral,decl,'r-',lw=3)
            ntot +=1
    rat = np.concatenate(rall)
    dect = np.concatenate(decll)
    fo = open('background/radec4back'+str(year)+'daystep'+str(daystep)+str(microns)+ido_view+'.txt','w')
    fo.write('| ra       | dec      | wavelength |  year  |  day   | obsloc | ido_view |\n')
    fo.write('| double   | double   | double     |  char  |  char  |  char  |  char    |\n')
    fo.write('| deg      | deg      | microns    |        |        |        |          |\n')


    while day < 365:
        for i in range(0,len(rat)):
            fo.write(str(rat[i])+' '+str(dect[i])+' '+str(microns)+' '+str(year)+' '+str(day)+' 0 '+ido_view+'\n')
        day += 10
    print(ntot)
    fo.close()
    return True  

def get_ebv():
    import fitsio
    import healpy as hp
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ebv = fitsio.read('pixweight-1-dark.fits')['EBV']
    ebvl = []
    ravl = []
    decvl = []
    ranvl = []
    decnvl = []

    nmiss = 0
    for sn in strips:
        for nn in range(0,nstrips[sn]):
            ral,decl = get_centers(sn,nn)
            th,phi = radec2thphi(np.array(ral),np.array(decl))
            px = hp.ang2pix(256,th,phi,nest=True)
            for i in range(0,len(px)):
                val = ebv[px[i]]
                if val != -1:
                    ebvl.append(val)
                    ravl.append(ral[i])
                    decvl.append(decl[i])
                else:
                    nmiss += 1
                    ranvl.append(ral[i])
                    decnvl.append(decl[i])
    print('number of missing pixels '+str(nmiss))
    print('number of included pixels '+str(len(ebvl)))
    print(min(ebvl),max(ebvl),np.median(ebvl),np.mean(ebvl))
    plt.scatter(ravl,decvl,c=ebvl)
    plt.plot(ranvl,decnvl,'k.')
    plt.title('centers in equatorial; colors are ebv values')
    plt.colorbar()
    plt.show()
    c = SkyCoord(np.array(ravl)* u.deg,np.array(decvl)* u.deg,frame='icrs')
    ec = c.transform_to('barycentricmeanecliptic')
    plt.scatter(ec.lon,ec.lat,c=ebvl)
    plt.title('centers in eclipitic; colors are ebv values')
    plt.colorbar()
    plt.show()
    
        

def plot_foot():
    ntot = 0
    for sn in strips:
        for nn in range(0,nstrips[sn]):
            ral,decl = get_corners(sn,nn)
            plt.plot(ral,decl,'r-',lw=3)
            ntot +=1
    print(ntot)  
    plt.xlabel('longitude coordinate (degrees)') 
    plt.ylabel('latitude coordinate (degrees)')    
    plt.title('Possible configuration for 155 Roman sectors (1700 sq. deg.)')  
    plt.savefig('footprint'+str(center_long)+str(center_lat)+'.png')
    plt.show()

def plot_foot_cen():
    ntot = 0
    for sn in strips:
        for nn in range(0,nstrips[sn]):
            ral,decl = get_centers(sn,nn)
            plt.plot(ral,decl,'rx',lw=3)
            ntot +=1
    print(ntot)  
    plt.xlabel('longitude coordinate (degrees)') 
    plt.ylabel('latitude coordinate (degrees)')    
    plt.title('Possible configuration for centers of 155 Roman sectors (1700 sq. deg.)')  
    plt.show()


def plot_foot_proj(proj='moll',coord='eq'):
    
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ntot = 0
    m = Basemap(projection=proj, lat_0=0, lon_0=0)    

    lat_line_range = [-90,90]
    lat_lines = 6
    lat_line_count = (lat_line_range[1]-lat_line_range[0])/lat_lines

    merid_range = [-180,180]
    merid_lines = 8
    merid_count = (merid_range[1]-merid_range[0])/merid_lines


    ct = 1001
    lats = np.arange(lat_line_range[0],lat_line_range[1],lat_line_count)
    m.drawparallels(lats,labels=[True,False,False,False],zorder=ct)
#     x0,y0 = m(-180,0)
#     for lat in lats:
#         x,y = m(-180,lat)
#         print(x,y)
#         xy = (x-0.5*x0,y)
#         print(xy)
#         plt.annotate(str(lat),xy=xy,xycoords='data',color='k',zorder=ct)
#         ct += 1

    ct += 1
    mers = np.arange(merid_range[0],merid_range[1],merid_count)
    m.drawmeridians(mers,zorder=ct)
    ct += 1
    for mer in mers:
        if mer < 0:
            mer = 360+mer
        plt.annotate(str(mer),xy=m(mer,-10),xycoords='data',color='k',zorder=ct)

    minlat = 90
    maxlat = -90
    for sn in strips:
        for nn in range(0,nstrips[sn]):
            ral,decl = get_corners(sn,nn)
            if coord != 'eq':
                c = SkyCoord(ral* u.deg,decl* u.deg,frame='icrs')
                ec = c.transform_to(coord)
                #print(ec)
                if coord != 'galactic':
                    xo,yo = m(np.array(ec.lon* u.deg),np.array(ec.lat* u.deg))
                    mn = np.min(np.array(ec.lat*u.deg))
                    if mn < minlat:
                        minlat = mn
                    mx = np.max(np.array(ec.lat*u.deg))
                    if mx > maxlat:
                        maxlat = mx
                    #print(np.min(ec.lat),np.max(ec.lat))
                else:
                    xo,yo = m(np.array(ec.l),np.array(ec.b))
            else:
                xo,yo = m(ral,decl)
                c = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='barycentricmeanecliptic')
                ec = c.transform_to('icrs')
                x,y = m(np.array(ec.ra),np.array(ec.dec))
                plt.plot(x,y,'r.',markersize=.5)
                gc = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='galactic')
                ecg = gc.transform_to('icrs')
                x,y = m(np.array(ecg.ra),np.array(ecg.dec))
                plt.plot(x,y,'b.',markersize=.5)


            plt.plot(xo,yo,'r-',lw=1)
            ntot +=1
    print(ntot,minlat,maxlat)  
    #plt.xlabel('longitude coordinate (degrees)') 
    #plt.ylabel('latitude coordinate (degrees)')    
    plt.title('Possible configuration for 155 Roman sectors (1700 sq. deg.)') 
    plt.gca().invert_xaxis() 
    plt.savefig('footprint_'+coord+proj+str(center_long)+str(center_lat)+'.png')
    plt.show()

def plot_back_proj(day=10,year=2027,microns=1.4,proj='moll',gs=10,vm=0.3,vx=0.7,daystep=10,back='total'):
    #d = np.loadtxt('back_'+str(year)+str(day)+'_'+str(microns)+'_out.txt').transpose()
    d = np.loadtxt('background/background_'+str(year)+'all_daystep'+str(daystep)+'_'+str(microns)+'.txt').transpose()
    sel = d[-5] == day
    if len(d[0][sel]) == 0:
        print('ERROR! no days found')
        return False
    ntot = 0
    m = Basemap(projection=proj, lat_0=0, lon_0=0)    
    #add ecliptic
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='barycentricmeanecliptic')
    ec = c.transform_to('icrs')
    x,y = m(np.array(ec.ra),np.array(ec.dec))
    plt.plot(x,y,'r.',markersize=.5)
    gc = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='galactic')
    ecg = gc.transform_to('icrs')
    x,y = m(np.array(ecg.ra),np.array(ecg.dec))
    plt.plot(x,y,'b.',markersize=.5)

    #plt.xlim(max(x), min(x))
    lat_line_range = [-90,90]
    lat_lines = 6
    lat_line_count = (lat_line_range[1]-lat_line_range[0])/lat_lines

    merid_range = [-180,180]
    merid_lines = 8
    merid_count = (merid_range[1]-merid_range[0])/merid_lines


    ct = 1001
    lats = np.arange(lat_line_range[0],lat_line_range[1],lat_line_count)
    m.drawparallels(lats,labels=[True,False,False,False],zorder=ct)
#     x0,y0 = m(-180,0)
#     for lat in lats:
#         x,y = m(-180,lat)
#         print(x,y)
#         xy = (x-0.5*x0,y)
#         print(xy)
#         plt.annotate(str(lat),xy=xy,xycoords='data',color='k',zorder=ct)
#         ct += 1

    ct += 1
    mers = np.arange(merid_range[0],merid_range[1],merid_count)
    m.drawmeridians(mers,zorder=ct)
    ct += 1
    for mer in mers:
        if mer < 0:
            mer = 360+mer
        plt.annotate(str(mer),xy=m(mer,-10),xycoords='data',color='k',zorder=ct)

    x,y = m(d[0][sel],d[1][sel])
    if back == 'total':
        backcol = -8
        titl = 'Total backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/background_'+str(year)+'_'+str(day)+'_'+str(microns)+'_new.png'
    if back == 'zodi':
        backcol = -12
        titl = 'Zodi backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/zodi_background_'+str(year)+'_'+str(day)+'_'+str(microns)+'_new.png'
    if back == 'ism':
        backcol = -11
        titl = 'ism backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/ism_background_'+str(year)+'_'+str(day)+'_'+str(microns)+'_new.png'
    if back == 'stars':
        backcol = -10
        titl = 'stars backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/stars_background_'+str(year)+'_'+str(day)+'_'+str(microns)+'_new.png'
   
    plt.scatter(x,y,c=d[backcol][sel],s=gs,marker='s',vmin=vm,vmax=vx)
    #plt.hexbin(x,y,d[-8],gridsize=gs)
    plt.colorbar(location='bottom')
    #plt.xlabel('longitude coordinate (degrees)') 
    #plt.ylabel('latitude coordinate (degrees)')    
    plt.title(titl)  
    plt.gca().invert_xaxis()
    plt.tight_layout()
    
    plt.savefig(plot_fn)
    #plt.show()
    plt.clf()
    return True

def plot_relnum(day=10,year=2027,microns=1.4,proj='moll',gs=20,vm=0.5,vx=1,daystep=10,back='total',newfoot='_newfoot',mdz ='_medianzodi',cmap='plasma',latmax=90, lat0=0, lon0=0,latmin=-90):
    #d = np.loadtxt('back_'+str(year)+str(day)+'_'+str(microns)+'_out.txt').transpose()
    d = np.loadtxt('data/background/background_'+str(year)+'all_daystep'+str(daystep)+'_'+str(microns)+newfoot+mdz+'.txt').transpose()
    if mdz == '':
        sel = d[-5] == day
    else:
        sel = np.zeros(len(d[-5]),dtype='bool')
        sel[:155] = 1
    if len(d[0][sel]) == 0:
        print('ERROR! no days found')
        return False
    ntot = 0
    #fig, ax = plt.subplots()
    
    m = Basemap(projection=proj, lat_0=lat0, lon_0=lon0,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=latmin,urcrnrlat=latmax)    
    #add ecliptic
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    c = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='barycentricmeanecliptic')
    ec = c.transform_to('icrs')
    x,y = m(np.array(ec.ra),np.array(ec.dec))
    fig = plt.figure()#figsize=(2.1, 1.6))
    plt.plot(x,y,'r.',markersize=.5)
    gc = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='galactic')
    ecg = gc.transform_to('icrs')
    x,y = m(np.array(ecg.ra),np.array(ecg.dec))
    plt.plot(x,y,'b.',markersize=.5)
    _,ymax = m(0,latmax)
    
    _,ymin = m(0,latmin)
    plt.ylim(ymin,ymax)
    #plt.xlim(max(x), min(x))
    
    lat_line_range = [latmin,latmax]
    latspace = 30
    lat_lines = int((latmax-latmin)/latspace)
    lat_line_count = (lat_line_range[1]-lat_line_range[0])/lat_lines

    merid_range = [-180,180]
    merid_lines = 8
    merid_count = (merid_range[1]-merid_range[0])/merid_lines


    ct = 1001
    lats = np.arange(latmin,latmax,latspace)#np.arange(lat_line_range[0],lat_line_range[1],lat_line_count)
    m.drawparallels(lats,labels=[True,False,False,False],zorder=ct,size=16)
#     x0,y0 = m(-180,0)
#     for lat in lats:
#         x,y = m(-180,lat)
#         print(x,y)
#         xy = (x-0.5*x0,y)
#         print(xy)
#         plt.annotate(str(lat),xy=xy,xycoords='data',color='k',zorder=ct)
#         ct += 1

    ct += 1
    mers = np.arange(merid_range[0],merid_range[1],merid_count)
    m.drawmeridians(mers,zorder=ct)
    ct += 1
    for mer in mers:
        if mer < 0:
            mer = 360+mer
        plt.annotate(str(mer),xy=m(mer,-10),xycoords='data',color='k',zorder=ct,size=16)

    x,y = m(d[0][sel],d[1][sel])
    titl = 'Variation in expected number density, based on '+str(back)+' background at '+str(microns)+' microns'
    if mdz == '':
        titl += '; '+str(year)+' day '+str(day)

    if back == 'total':
        backcol = -8
        plot_fn = 'data/background/plots/background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_new.png'
    if back == 'zodi':
        backcol = -12
        #titl = 'Zodi backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/zodi_background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_new.png'
    if back == 'ism':
        backcol = -11
        #titl = 'ism backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/ism_background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_new.png'
    if back == 'stars':
        backcol = -10
        #titl = 'stars backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/stars_background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_new.png'
   
    minback = np.min(d[-8][sel]) #minimum overall
    minsp = np.min(d[backcol][sel]) #minimum for quantity of interest
    varsp = d[backcol][sel]-minsp+minback
    mxnum = np.max(1./varsp)
    relnum = np.sqrt(1./varsp/mxnum)
    print(min(relnum),max(relnum))
    #x, y,z = np.meshgrid(d[0][sel],d[1][sel],relnum)
    #nx = len(np.unique(x))
    #ny = len(np.unique(y))
    #z = np.zeros((nx,ny))
    #for ii in range(0,len(x)):
    #    xind = int(.999*nx*(x[ii]-min(x))/max(x))
    #    yind = int(.999*ny*(y[ii]-min(y))/max(y))
    #    z[xind][yind] = relnum[ii]
    im1 = plt.scatter(x,y,c=relnum,s=gs,marker='s',vmin=vm,vmax=vx,cmap=cmap)
    #z = np.reshape(relnum,x.shape)
    #xg,yg = np.meshgrid(np.unique(x),np.unique(y))
    #im1 = m.pcolormesh(xg,yg,z.transpose(),cmap=cmap)
    #plt.hexbin(x,y,d[-8],gridsize=gs)
    #plt.colorbar(location='bottom')
    cbar = m.colorbar(im1,location='bottom')
    cbar.ax.tick_params(labelsize=18)
    
    #plt.xlabel('longitude coordinate (degrees)') 
    #plt.ylabel('latitude coordinate (degrees)')    
    #plt.title(titl,wrap=True)  
    
    
    plt.xlabel(titl,wrap=True,labelpad=40,size=18)#,location='top')
    #ax.xaxis.set_label_position('top') 
    plt.gca().invert_xaxis()
    
    plt.tight_layout()
    if mdz == '':
        plt.text(0.1, 0.95, 'Day '+str(day)+', '+str(year),fontsize=16, horizontalalignment='center',verticalalignment='center', transform=plt.gca().transAxes)
    plt.savefig(plot_fn)
    #plt.show()
    plt.clf()
    return True

def plot_relnum_half(day=10,year=2027,microns=1.4,proj='moll',gs=20,vm=0.5,vx=1,daystep=10,back='total',newfoot='_newfoot',mdz ='_medianzodi',cmap='Greys',latmax=90, lat0=0, lon0=0,latmin=-90,height=2,ratio=3):
    #d = np.loadtxt('back_'+str(year)+str(day)+'_'+str(microns)+'_out.txt').transpose()
    plt.clf()
    d = np.loadtxt('data/background/background_'+str(year)+'all_daystep'+str(daystep)+'_'+str(microns)+newfoot+mdz+'.txt').transpose()
    if mdz == '':
        sel = d[-5] == day
    else:
        sel = np.zeros(len(d[-5]),dtype='bool')
        sel[:155] = 1
    if len(d[0][sel]) == 0:
        print('ERROR! no days found')
        return False
    ntot = 0
    #fig, ax = plt.subplots()
    
    m = Basemap(projection=proj, lat_0=lat0, lon_0=lon0,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=latmin,urcrnrlat=latmax)    
    #add ecliptic
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    c = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='barycentricmeanecliptic')
    ec = c.transform_to('icrs')
    x,y = m(np.array(ec.ra),np.array(ec.dec))
    fig = plt.figure(figsize=(ratio*height, height))
    plt.plot(x,y,'r.',markersize=.5)
    gc = SkyCoord(np.arange(-180,180)* u.deg,np.zeros(360)* u.deg,frame='galactic')
    ecg = gc.transform_to('icrs')
    x,y = m(np.array(ecg.ra),np.array(ecg.dec))
    plt.plot(x,y,'b.',markersize=.5)
    _,ymax = m(0,latmax)
    
    _,ymin = m(0,latmin)
    
    xmin,_ = m(180,0)
    plt.ylim(ymin,ymax)
    tx,ty = xmin,ymin-0.1*(ymax-ymin)
    titl = r'Relative expected $n_{\rm gal}$, based on '+str(back)+' background '#at '+str(microns)+' microns'
    if mdz == '':
        #titl += '; '+str(year)+' day '+str(day)
        titl += '; day '+str(day)

    plt.text(tx,ty,titl,size=16,wrap=True)#,labelpad=40,size=18)

    plt.axis('off')
    #plt.xlim(max(x), min(x))
    lats4lon = np.arange(latmin,latmax,.1)
    merl = [-180,180]
    for mer in merl:
        mers = np.ones(len(lats4lon))*mer
        x,y = m(mers,lats4lon)
        plt.plot(x,y,'k-')
    lonspace = 45
    merl = np.arange(-180+lonspace,180,lonspace)
    for mer in merl:
        mers = np.ones(len(lats4lon))*mer
        x,y = m(mers,lats4lon)
        plt.plot(x,y,'k:')
        plt.annotate(str(mer),xy=m(mer,-10),xycoords='data',color='k',size=16)#,zorder=ct,size=16)
    latspace = 30
    latl = np.arange(latmin,latmax,latspace)
    lons4lat = np.arange(-180,180,.1)
    for lat in latl:
        lats = np.ones(len(lons4lat))*lat
        x,y = m(lons4lat,lats)
        plt.plot(x,y,'k:')
        x,y = m(-185,lat)
        xy = (x-0.5,y)
        if lat != -90:
            plt.annotate(str(lat),xy=xy,xycoords='data',color='k',size=16)#,zorder=ct)
#         ct += 1


    x,y = m(d[0][sel],d[1][sel])

    if back == 'total':
        backcol = -8
        plot_fn = 'data/background/plots/background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_half.png'
    if back == 'zodi':
        backcol = -12
        #titl = 'Zodi backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/zodi_background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_half.png'
    if back == 'ism':
        backcol = -11
        #titl = 'ism backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/ism_background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_half.png'
    if back == 'stars':
        backcol = -10
        #titl = 'stars backgroud (MJy/sr) at '+str(microns)+' microns; '+str(year)+' day '+str(day)
        plot_fn = 'data/background/plots/stars_background_'+str(year)+'_'+str(day)+'_'+str(microns)+newfoot+mdz+'_'+cmap+'_half.png'
   
    minback = np.min(d[-8][sel]) #minimum overall
    minsp = np.min(d[backcol][sel]) #minimum for quantity of interest
    varsp = d[backcol][sel]-minsp+minback
    mxnum = np.max(1./varsp)
    relnum = 1./varsp/mxnum
    print(min(relnum),max(relnum))
    #x, y,z = np.meshgrid(d[0][sel],d[1][sel],relnum)
    #nx = len(np.unique(x))
    #ny = len(np.unique(y))
    #z = np.zeros((nx,ny))
    #for ii in range(0,len(x)):
    #    xind = int(.999*nx*(x[ii]-min(x))/max(x))
    #    yind = int(.999*ny*(y[ii]-min(y))/max(y))
    #    z[xind][yind] = relnum[ii]
    im1 = plt.scatter(x,y,c=relnum,s=gs,marker='s',vmin=vm,vmax=vx,cmap=cmap)
    #z = np.reshape(relnum,x.shape)
    #xg,yg = np.meshgrid(np.unique(x),np.unique(y))
    #im1 = m.pcolormesh(xg,yg,z.transpose(),cmap=cmap)
    #plt.hexbin(x,y,d[-8],gridsize=gs)
    #plt.colorbar(location='bottom')
    cbar = m.colorbar(im1,location='bottom',pad=.25)
    cbar.ax.tick_params(labelsize=16)
    
    #plt.xlabel('longitude coordinate (degrees)') 
    #plt.ylabel('latitude coordinate (degrees)')    
    #plt.title(titl,wrap=True)  
    
    
    #plt.xlabel(titl,wrap=True,labelpad=40,size=18)#,location='top')
    #ax.xaxis.set_label_position('top') 
    plt.gca().invert_xaxis()
    
    plt.tight_layout()

    #plt.show()
    #return True

    
    plt.savefig(plot_fn)
    #plt.show()
    plt.clf()
    return True


def mkbackgif(back='total'):
    frame_length = 0.5 # seconds between frames
    end_pause = 4 # seconds to stay on last frame
    # loop through files, join them to image array, and write to GIF called 'wind_turbine_dist.gif'
    images = []
    pre_fn = ''
    sorted_files = [pre_fn+'background_2027_1_1.4.png']
    vm = 0.3
    vx = 0.7
    if back != 'total':
        pre_fn = back+'_'
        sorted_files = []
        #if back == 'zodi':
        vm = 0
        vx = 0.5
    day = 10
    while day < 370:
        fn = pre_fn+'background_2027_'+str(day)+'_1.4.png'
        suc = plot_back_proj(day,back=back,vm=vm,vx=vx)
        if os.path.isfile(fn):
        #try:
            #plot_back_proj(day)
            sorted_files.append(fn)
        #except:
        else:
            print('no '+fn)   
        day += 10     
    #,'background_2027_10_1.4.png','background_2027_20_1.4.png','background_2027_90_1.4.png','background_2027_180_1.4.png']
    for ii in range(0,len(sorted_files)):       
        file_path = sorted_files[ii]
        if ii==len(sorted_files)-1:
            for jj in range(0,int(end_pause/frame_length)):
                images.append(imageio.imread(file_path))
        else:
            images.append(imageio.imread(file_path))
    # the duration is the time spent on each image (1/duration is frame rate)
    gif_name = pre_fn+'background.gif'
    imageio.mimsave(gif_name, images,'GIF',duration=frame_length)


def get_strip(lat):
    return ((lat-min_lat)//sec_size)

def get_nstrip(lat):
    st = np.array(get_strip(lat).astype(int))
    #print(st)
    nstrip = nstrips[st]
    return nstrip
    
def in_foot(lat,long):
    #return the points that are in the footprint
    #lat should be the latitude coordinate in degrees (-90,90)
    #long should be the longitude coordinate in degrees (-180,180)
    sel = lat > min_lat
    sel &= lat < max_lat
    sn = get_nstrip(lat[sel])
    minl,maxl = get_long_minmax(lat[sel],sn)
    sell = long[sel] > minl
    sell &= long[sel] < maxl
    return lat[sel][sell],long[sel][sell]

if __name__ == '__main__':
    #do a test on N randoms and plot kept points in red
    import sys
    nran = int(sys.argv[1])
    acosl = np.random.rand(nran)*2-1 #distribute randomly in arccos
    ral = np.random.rand(nran)*360-180
    decl = np.arccos(acosl)
    decl = 180/np.pi*(decl-np.pi/2.)
    lat,long = in_foot(decl,ral)
    area = len(long)/len(ral)*360*360/np.pi
    print('based on the number of points kept, the implied area is '+str(round(area,1))+' square degrees')
    plt.plot(ral,decl,'k,')
    plt.plot(long,lat,'r,')
    plt.xlabel('longitude coordinate (degrees)')
    plt.ylabel('latitude coordinate (degrees)')
    plt.show()
    