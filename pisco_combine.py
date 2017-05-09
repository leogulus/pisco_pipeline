import numpy as np
import matplotlib
import os

from astropy.io import fits
import subprocess
import cosmics
import shlex
import sys
#from pisco_lib import *
#edited 4/4/17

"""
pisco_combine: run pisco pipeline to reduce the raw data to clean data with correct WCS
The pipeline is a combination of LA Cosmics, Astrometry, Sextractor, SCAMP and SWARP.

ARGUMENTS:
1. raw directory (e.g., 'ut170103/')
2. fieldname for object (e.g., 'Field027')

Examples: python pisco_pipeline/pisco_combine.py data/ Field026
"""
def filter_name(index):
    """
    filter_name: turn index [1,8] into letter band (g,r,i,z) for PISCO quadrant data
    INPUT:
    - index: number
    OUTPUT:
    - a pair of letter for corresponding band and dome band
    """
    if index==1 or index==2:
        filter_name='g'
        dome_name='g'
    elif index==3 or index==4:
        filter_name='r'
        dome_name='r'
    else:
        dome_name='iz'
        if index==5 or index==6:
            filter_name='i'
        elif index==7 or index==8:
            filter_name='z'
    return [filter_name,dome_name]

def list_file_name(dir,name,end=0):
    """
    list_file_name: list all filename which started with 'name' and end with
    'end' in 'dir' directory
    INPUT:
    - dir: directory to search in
    - name: begining of the file name
    - end: ending of the file name
    OUTPUT:
    - list of all filename in that directory
    """
    names=[]
    for file in os.listdir(dir):
        if file.startswith(name):
            if end==0:
                names.append(os.path.join(dir, file))
            else:
                if file.endswith(end):
                    names.append(os.path.join(dir, file))
    if len(names)==0:
        print 'Cannot find the files'
    return names

def open_files(names,index,bias=np.array([])):
    """
    open_files: use to open multiple bias or domeflat files at once and take the mean to
    to get the average bias/domeflat file for image reduction
    bias: take the mean of all bias files
    domeflat: subtracted by average 'bias' (also calcualted from open_files) before take the mean
    INPUT:
    - name: starting name of the bias/domeflat files (output from 'list_file_name')
    - index: extension of the fits file to read in (8 extension of PISCO - two each for different bands)
    - (optional) bias: average 2D bias image (required to calculate domeflat correctly)
    OUTPUT:
    - 2D array of average bias/domeflat images
    """
    ch_bs=[]
    for name in names:
        hdulist = fits.open(name)
        ch_b=hdulist[index].data
        if len(bias)==0:
            ch_bs.append(ch_b) #for bias to combine as a mean
        else:
            ch_bs.append(ch_b-bias) #for domeflat-bias before combine into a mean
    return np.mean(np.array(ch_bs),axis=0)

def plot_one_chip(ax,data,vmin,vmax):
    """
    plot_one_chip: plot 2D array 'data' on 'ax' with Normalize scale 'vmin' and 'vmax'
    INPUT:
    - ax: ax from fig,ax=plt.subplots(...)
    - data: 2D array data to be plot
    - vmin: normalize scale for the minimum
    - vmax: normalize scale for the maximum
    OUTPUT:
    - plot of the data at a specified normalize scale
    """
    norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    c_m = matplotlib.cm.Greys
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    ax.imshow(data, cmap=c_m, norm=norm)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

def save_fits(index,dir,outdir,fieldname,final_image,name):
    """
    save_fits: save the fits file from 2D array 'final_image' with a known header from the raw PISCO data
    'fieldname' (changed the size to accompany the attachment of two amplifiers) with the output 'name'
    in 'reduced/' directory
    INPUT:
    - index: specific the band (g, r, i, z) that we want to save on.
    - dir: input directory for raw PISCO data
    - fieldname: starting of the name for the raw PISCO data (e.g., 'Field027_B_73')
    - final_image: 2D array of image that we want to save to the fits file
    - name: output name of the fits file in 'reduced/' directory
    OUTPUT:
    - fits file in 'reduced/' directory
    """
    ch1_name=list_file_name(dir,fieldname)
    hdulist = fits.open(ch1_name[0])

    hdu0 = hdulist[0]
    hdu0.header['NAXIS']=2
    hdulist[index].header['NAXIS1']='1546'
    hdulist[index].header['DATASEC'] = '[1:1546,1:3092]'
    hdulist[index].header['TRIMSEC'] = '[1:1546,1:3092]'
    hdulist[index].header['ORIGSEC'] = '[1:1546,1:3092]'
    hdulist[index].header['CCDSEC'] = '[1:1546,3093:6184]'
    hdulist[index].header['DETSEC'] = '[1:1546,3093:6184]'

    hdu1 = fits.ImageHDU(final_image*1000,name='filter '+filter_name(index)[0],header=hdulist[index].header)
    hdu_l = fits.HDUList(hdus=[hdu0,hdu1])
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)
    outname=os.path.join(outdir,name)
    print 'saving the fits file '+outname
    hdu_l.writeto(outname, overwrite=True)

    data, header = fits.getdata(outname, header=True)
    fits.writeto(outname, data, header, clobber=True)

def reduce_data(dir,index,fieldname):
    """
    reduce_data: combine raw PISCO data with bias and domeflat to create 2D array of output image
    using function list_file_name, open_files
    INPUT:
    - dir: directory for the raw PISCO data
    - index: index for the band of the image that we want to reduce
    - fieldname: the begining of the file name (e.g., 'Field027_B_73')
    - (extra) cut: -27 is the number of pixel needed to be cut out for the gap in the image
    OUTPUT:
    - ch1: 2D array of raw input image
    - bias: 2D array for the bias image
    - domeflat: 2D array for the domeflat image
    - img: 2D array of the output image after subtraction of bias and normalization with domeflat
    """
    cut=-27

    ch1_name=list_file_name(dir,fieldname)
    print 'working on %s with the index=%i'%(ch1_name[0],index)
    hdulist = fits.open(ch1_name[0])
    ch1=hdulist[index].data

    bias_names=list_file_name(dir,'Bias_')
    domeflat_names=list_file_name(dir,"domeflat"+filter_name(index)[1])

    bias=open_files(bias_names,index)
    domeflat=open_files(domeflat_names,index,bias=bias)
    domeflat[domeflat==0]=1e-4
    img = (ch1-bias)/domeflat;
    ch1,bias,domeflat,img=ch1[:,:cut],bias[:,:cut],domeflat[:,:cut],img[:,:cut]

    if index%2==0:
        return np.fliplr(ch1), np.fliplr(bias), np.fliplr(domeflat), np.fliplr(img)
    else:
        return ch1, bias, domeflat, img

def cosmic_reduce(dir,field,band):
    """
    cosmic_reduce: read the FITS file and use L.A. Cosmic (http://www.astro.yale.edu/dokkum/lacosmic/)
    to remove cosmic rays in the images
    INPUT:
    - dir: directory input of the combine images ('reduced/')
    - field: beginning of the file name (e.g., 'Field027_A_72')
    - band: {'g','r','i','z'} band
    PARAMETERS for LA Cosmic:
    - gain and readnoise are the property from the telescope (PISCO: gain 4 ADU/e, readnoise 3 e -Brian[3/27/17])
    - satlevel: identify saturated level for bright stars
    - sigclip, sigfrac, objlim
    OUTPUT:
    - nField..._g.fits: not clean data (original with a mask cut)
    - cField..._g.fits: clean version, removed cosmic ray
    - mField..._g.fits: masked file to remove cosmic ray
    """
    array, header = cosmics.fromfits(os.path.join(dir,field+'_'+band+'.fits'))
    #cutting the circular aperature of the image out to only have good pixels in the center
    if band=='g':
        array_c=array[10:-10,350:2500]
    elif band=='r':
        array_c=array[10:-10,350:2500]
    elif band=='z':
        array_c=array[10:-10,650:2750]
    elif band=='i':
        array_c=array[10:-10,650:2750]
    c = cosmics.cosmicsimage(array_c, gain=4.0, readnoise=3.0, sigclip = 2.5, sigfrac = 0.5,\
                             objlim = 5.0, satlevel=3000.0, verbose=False)
    c.run(maxiter = 5)

    # if not os.path.exists(os.path.join(dir,'cosmics')):
    #     os.makedirs(os.path.join(dir,'cosmics'))
    cosmics.tofits(os.path.join(dir,'cosmics','n'+field+'_'+band+'.fits'), array_c, header)
    cosmics.tofits(os.path.join(dir,'cosmics','c'+field+'_'+band+'.fits'), c.cleanarray, header)
    cosmics.tofits(os.path.join(dir,'cosmics','m'+field+'_'+band+'.fits'), c.mask, header)


def astrometry_solve(cosmicdir,field,outdir):
    """
    astrometry_solve: apply astrometry algorithm to find celestial Coordinate (WCS) for the image
    REQUIRE:
    - appropriate index files in '/usr/local/astrometry/data/' for Astrometry to have enough patch of
    the sky to search for the position
    INPUT:
    - cosmicdir: input directory for cosmic-ray subtracted fits file ('reduced/cosmics/')
    - field: begining of the file name after cosmic-ray subtraction for each band and each exposure
     (e.g. 'cField027_B_73_z')
    - outdir: output directory for these outputs ('wcs/')
    OUTPUT:
    - wcs/.wcs file: for the header with appropriate coordinate
    - new_fits/..._new.fits: updated fits file with new wcs information in 'new_fits' directory
    """
    # if not os.path.exists(outdir):
    #     os.makedirs(os.path.join(outdir))
    if not os.path.isfile(os.path.join(outdir,field+'.wcs')):
        cmd='solve-field %s --downsample 2 --overwrite --scale-unit arcsecperpix --scale-low 0.08 --scale-high 0.3 --dir %s' \
        % (os.path.join(cosmicdir,field+'.fits'),outdir)
        print cmd
        sub=subprocess.check_call(shlex.split(cmd))
        if sub==0:
            print 'finish solve-field and updating fits headers'
        else:
            print 'solve-field does not work.'

    orig=fits.open(os.path.join(cosmicdir,field+'.fits'))
    wcs_file=fits.open(os.path.join(outdir,field+'.wcs'))
    header=wcs_file[0].header
    wcsaxes_index=np.where(np.array(header.keys())=='WCSAXES')[0][0]
    for i in range(wcsaxes_index,len(header)):
        orig[0].header[header.keys()[i]]=header.values()[i]
    orig.writeto(os.path.join('new_fits',field+'_new.fits'),clobber=True)

def sextracting(field):
    """
    sextracting: run Sextractor to find all the point sources in .ldac.fits format (suitable for SCAMP input)
    INPUT:
    - config.sex: sextractor config file
    - field: begining of the file name for each band and each exposure (e.g. 'cField027_B_73_z')
    OUTPUT:
    - new_fits/..._new.ldac.fits: source catalogs of all the point source from Sextractor
    """
    cmd='sex %s -c pisco_pipeline/config.sex -CATALOG_NAME %s' % \
    (os.path.join('new_fits',field+'_new.fits'),os.path.join('new_fits',field+'_new.ldac.fits'))
    print cmd
    sub=subprocess.check_call(shlex.split(cmd))

def scamp(fieldname):
    """
    scamp: run SCAMP to align the coordinate better after Astrometry with distortions.
    (need to run all exposure and all filters at once to get the best performance)
    INPUT:
    - config.scamp: SCAMP config file
    - fieldname: begining of the file name (e.g., 'cField027')
    OUTPUT:
    - new_fits/...ldac.head: SCAMP output which includes new celestial coordinates for fixing WCS
    """
    cmd='scamp %s -c pisco_pipeline/config.scamp' % ' '.join(list_file_name('new_fits',fieldname,end='_new.ldac.fits'))
    print cmd
    sub=subprocess.check_call(shlex.split(cmd))

def swarp(fieldname):
    """
    swarp: run SWARP to combine multiple exposure into a better image with SCAMP output to help correct the location
    INPUT:
    - config.swarp: SWARP config file
    - fieldname: begining of the file name (e.g., 'cField027')
    OUTPUT:
    - final/coadd_'fieldname'_'g'.fits: final image for each 'band' with corrected WCS
    """
    bands=['g','r','i','z']
    print 'Swarping...'
    for band in bands:
        coadd_files=list_file_name('new_fits',fieldname,end=band+'_new.fits')
        cmd='swarp %s %s -c pisco_pipeline/config.swarp -IMAGEOUT_NAME %s' %\
        (coadd_files[0],coadd_files[1],os.path.join('final','coadd_'+fieldname+'_'+band+'.fits'))
        print cmd
        sub=subprocess.check_call(shlex.split(cmd))
##--------

if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    ## Pipeline to run PISCO reduction data
    dir=str(sys.argv[1])
    fieldname=str(sys.argv[2])
    outdir='wcs'
    reducedir='reduced'
    cosmicdir=os.path.join(reducedir,'cosmics')

    if not os.path.exists(outdir):
        os.makedirs(os.path.join(outdir))
    if not os.path.exists(reducedir):
        os.makedirs(reducedir)
    if not os.path.exists(cosmicdir):
        os.makedirs(cosmicdir)
    if not os.path.exists('new_fits'):
        os.makedirs(os.path.join('new_fits'))
    if not os.path.exists('final'):
        os.makedirs('final')

    fields=[name.split('/')[1].split('.')[0] for name in list_file_name(dir,fieldname)]
    for field in fields:
        for index in [1,3,5,7]:
            ch1, bias1, domeflat1, img1=reduce_data(dir,index,field)
            ch2, bias2, domeflat2, img2=reduce_data(dir,index+1,field)
            final_image=np.concatenate((img1,img2),axis=1)
            save_fits(index,dir,reducedir,field,final_image,"%s_%s.fits"%(field,filter_name(index)[0]))

    #Cosmic ray reduction using L.A. Cosmic
    bands=['g','r','i','z']
    for field in fields:
        for band in bands:
            if not os.path.isfile(os.path.join(reducedir,'cosmics','c'+field+'_'+band+'.fits')):
                print 'working on the cosmic '+'c'+field+'_'+band
                cosmic_reduce(reducedir,field,band)
            else:
                print 'already did this band '+band
                #cosmic_reduce(reducedir,field,band)

    cfieldname='c'+fieldname
    print 'number of files in %s is %i' % (cosmicdir,len(list_file_name(cosmicdir,cfieldname)))

    for field_long in list_file_name(cosmicdir,cfieldname):
        field=field_long.split('/')[2].split('.')[0]
        #Astrometry to get a rough estimate on the World Coordinate System (WCS) for each images
        astrometry_solve(cosmicdir,field,outdir)
        #Sextracting
        sextracting(field)

    #SCAMP
    scamp(cfieldname)
    #SWARP
    swarp(cfieldname)
