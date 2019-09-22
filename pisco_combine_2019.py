import os, re, subprocess, shlex, sys, yaml

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
import cosmics

#from pisco_lib import *
# edited 5/9/17

"""
pisco_combine: run pisco pipeline to reduce the raw data to clean data with correct WCS
The pipeline is a combination of LA Cosmics, Astrometry, Sextractor, SCAMP and SWARP.

ARGUMENTS:
1. raw directory (e.g., 'ut170103/')
2. fieldname for object (e.g., 'Field027')

Examples: python pisco_pipeline/pisco_combine.py data/ Field026
python pisco_pipeline/pisco_combine_2019.py ut170619/ Field292 'twilight'
python pisco_pipeline/pisco_combine_2019.py ut171209/ CHIPS0152-5028 'twilight'
python pisco_pipeline/pisco_combine_2019.py ut171209/ CHIPS1011-0505 'twilight'
"""

def list_file_name(dir, name, end=''):
    names=[]
    myReg=re.compile('^%s((?:[^/]*/)*)(.*)%s'%(name,end))
    for text in os.listdir(dir):
        if myReg.search(text) != None:
            names.append(os.path.join(dir, myReg.search(text).group()))
    return names

def open_files(names, index, bias=np.array([]), twilight=False):
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
    ch_bs = []
    for name in names:
        hdulist = fits.open(name)
        ch_b = hdulist[index].data
        if len(bias) == 0:
            ch_bs.append(ch_b)  # for bias to combine as a mean
        else:
            ch_bs.append(ch_b - bias)  # for domeflat-bias before combine into a mean
    if twilight == True:
        print 'working on twlight flat'
        return np.median(np.array(ch_bs), axis=0)
    else:
        return np.mean(np.array(ch_bs), axis=0)

def reduce_data(dir, index, fieldname, flat='domeflat'):
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
    cut=-32
    print cut

    ch1_name = list_file_name(dir, fieldname)
    print 'working on %s with the index=%i' % (ch1_name[0], index)
    hdulist = fits.open(ch1_name[0])
    ch1 = hdulist[index].data

    bias_names = list_file_name(dir, 'Bias_')
    domeflat_names = list_file_name(dir, "twiflat_")
    
    bias = open_files(bias_names, index)
    domeflat = open_files(domeflat_names, index, bias=bias, twilight=True)

    domeflat[domeflat == 0] = 1e-4

    img = (ch1 - bias) / (domeflat / np.median(domeflat))
    print ch1.shape, bias.shape, domeflat.shape, np.median(ch1), np.median(bias), np.median(domeflat),\
        np.median(domeflat / np.median(domeflat)), np.median(img)

    ch1, bias, domeflat, img = ch1[:, :cut], bias[:, :cut], domeflat[:, :cut], img[:, :cut]
    if index % 2 == 0:
        return np.fliplr(ch1), np.fliplr(bias), np.fliplr(domeflat), np.fliplr(img)
    else:
        return ch1, bias, domeflat, img

def save_fits(index, dir, outdir, fieldname, final_image, name):
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
    ch1_name = list_file_name(dir, fieldname)
    hdulist = fits.open(ch1_name[0])

    hdu0 = hdulist[0]
    hdu0.header['NAXIS'] = 2
    hdulist[index].header['NAXIS1'] = '1546'
    hdulist[index].header['DATASEC'] = '[1:3094,1:6147]'
    hdulist[index].header['TRIMSEC'] = '[1:3094,1:6147]'
    hdulist[index].header['ORIGSEC'] = '[1:3094,1:6147]'
    hdulist[index].header['CCDSEC'] = '[1:3094,6147:6184]'
    hdulist[index].header['DETSEC'] = '[1:3094,6147:6184]'

    hdu1 = fits.ImageHDU(final_image, name='filter ' +
                         filter_name(index)[0], header=hdulist[index].header) #*1000 (edit: add 1000)
    hdu_l = fits.HDUList(hdus=[hdu0, hdu1])
    outname = os.path.join(outdir, name)
    print 'saving the fits file ' + outname
    hdu_l.writeto(outname, overwrite=True)

    data, header = fits.getdata(outname, header=True)
    fits.writeto(outname, data, header, overwrite=True)
    
def filter_name(index):
    """
    filter_name: turn index [1,8] into letter band (g,r,i,z) for PISCO quadrant data
    INPUT:
    - index: number
    OUTPUT:
    - a pair of letter for corresponding band and dome band
    """
    if index == 1 or index == 2:
        filter_name = 'g'
        dome_name = 'g'
    elif index == 3 or index == 4:
        filter_name = 'r'
        dome_name = 'r'
    else:
        dome_name = 'iz'
        if index == 5 or index == 6:
            filter_name = 'i'
        elif index == 7 or index == 8:
            filter_name = 'z'
    return [filter_name, dome_name]

def cosmic_reduce(dir, field, band):
    if not os.path.isfile(os.path.join(dir, 'cosmics', 'c' + field + '_' + band + '.fits')):
        print 'working on the cosmic ' + 'c' + field + '_' + band
        array, header = cosmics.fromfits(
            os.path.join(dir, field + '_' + band + '.fits'))
        print os.path.join(dir, field + '_' + band + '.fits')

        satfield = '_field2'

        with open("pisco_pipeline/params.yaml", 'r') as stream:
            try:
                param=yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        if band == 'g':
            array_c = array[35:-20,600:-1050]  # [20:-20,350:2550]
            satlv = param['satur_level%s_%s'%(satfield,band)] #2000.0
        elif band == 'r':
            array_c = array[50:-20,750:-850]  # [20:-20,350:2550]
            satlv = param['satur_level%s_%s'%(satfield,band)]# 1250.0
        elif band == 'i':
            array_c = array[20:-70,1200:-520]  # [20:-20,650:2800]
            satlv = param['satur_level%s_%s'%(satfield,band)] #600.0
        elif band == 'z':
            array_c = array[10:-20,960:-670]  # [20:-20,650:2800]
            satlv = param['satur_level%s_%s' % (satfield, band)]  # 1500.0
        c = cosmics.cosmicsimage(array_c, gain=0.25, readnoise=3.0, sigclip=6.0, sigfrac=0.4, 
                                 objlim=5.0, satlevel=satlv, verbose=False)  # sigclip=8.0, sigfrac=0.4

        #sigclip 4.5, sigfrac=0.5, objlim=2, niter=4
        #IMAC: sigclip: 6, niter=4, objlim=5.0, sigfrac=0.3, (gain/readnoise from PISCO)
        c.run(maxiter=4)
        cosmics.tofits(os.path.join(dir, 'cosmics', 'c' + field + '_' + band + '.fits'), c.cleanarray, header)
#         cosmics.tofits(os.path.join(dir, 'cosmics', 'm' + field + '_' + band + '.fits'), c.mask, header)
#         cosmics.tofits(os.path.join(dir, 'cosmics', 'n' + field + '_' + band + '.fits'), array_c, header)
    else:
        print 'already did the cosmic with this band ' + band

def astrometry_solve(cosmicdir, field, outdir):
    myReg=re.compile(r'(CHIPS).+?(?=\_)')
    fieldc=myReg.search(field).group()

    myReg2=re.compile(r'(cCHIPS)\w+.\w+')
    field_all=myReg2.search(field_location).group()
    
    if not os.path.isfile(os.path.join(outdir, field_all+ '.wcs')):
        total=pd.read_csv('/Users/taweewat/Documents/xray_project/ned-result/total_776_new_pan.csv',index_col=0)


        if fieldc=='CHIPS1423-3125':
            ra,dec=215.88332999999997,-31.423749999999998
        elif fieldc=='CHIPS2007-4434':
            ra,dec=301.98419,-44.58055
        elif fieldc=='CHIPS1605-3115':
            ra,dec=241.4466667,-31.2577778
        else:
            ra,dec=total.loc[total['chips']==fieldc,'ra_01'].values[0], total.loc[total['chips']==fieldc,'dec_01'].values[0]
        
        print ra, dec, field, fieldc
        cmd = 'solve-field %s --downsample 2 --overwrite --scale-unit arcsecperpix --scale-low 0.08 --scale-high 0.3 --dir %s --ra %s --dec %s --radius 2' \
            % (os.path.join(cosmicdir, field_all+ '.fits'), outdir, str(ra), str(dec))
        print cmd

        sub = subprocess.check_call(shlex.split(cmd))
        if sub == 0:
            print 'finish solve-field and updating fits headers'
        else:
            print 'solve-field does not work.'
    else:
        print 'already have ' + field_all + '.wcs'

    orig = fits.open(os.path.join(cosmicdir, field_all+ '.fits'))
    wcs_file = fits.open(os.path.join(outdir, field_all+ '.wcs'))
    header = wcs_file[0].header
    wcsaxes_index = np.where(np.array(header.keys()) == 'WCSAXES')[0][0]
    for i in range(wcsaxes_index, len(header)):
        orig[0].header[header.keys()[i]] = header.values()[i]
    orig.writeto(os.path.join('new_fits', field_all + '_new.fits'), overwrite=True)

def sextracting(field, band, dir):
    with open("../pisco_code/pisco_pipeline/params.yaml", 'r') as stream:
        try:
            param=yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)

    satfield = '_field2'
    
    myReg3=re.compile(r'(CHIPS)[^\_]*\_[^\_]*')
    seeing = float(fits.open(list_file_name(dir,myReg3.search(field).group())[0])[0].header['FWHM1'])
    print seeing
    
    pxscale=0.12
    cmd = 'sex %s -c pisco_pipeline/config.sex -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -CHECKIMAGE_NAME %s,%s -PIXEL_SCALE %s' % \
        (os.path.join('new_fits', field + '_new.fits'),
         os.path.join('new_fits', field + '_new.ldac.fits'), str(seeing), str(param['satur_level%s_%s'%(satfield,band)]),'check_%s.fits'%(band),'segment_%s.fits'%(band),str(pxscale))
    print cmd
    subprocess.check_call(shlex.split(cmd))

    cmd = 'sex %s -c pisco_pipeline/config.sex -CATALOG_NAME %s -CATALOG_TYPE ASCII -SEEING_FWHM %s -SATUR_LEVEL %s -CHECKIMAGE_NAME %s,%s -PIXEL_SCALE %s' % \
        (os.path.join('new_fits', field + '_new.fits'),
         os.path.join('new_fits', 'tmp_%s.cat' % band), str(seeing), str(param['satur_level%s_%s'%(satfield, band)]), 'check_%s.fits'%(band), 'segment_%s.fits'%(band),str(pxscale))
    print cmd
    subprocess.check_call(shlex.split(cmd))

    name = ['NUMBER', 'EXT_NUMBER', 'XWIN_WORLD', 'YWIN_WORLD', 'MAG_AUTO', 'MAGERR_AUTO', 'MAG_APER', 'MAGERR_APER', 'XWIN_IMAGE',
            'YWIN_IMAGE', 'ERRAWIN_IMAGE', 'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE', 'FLUX_AUTO', 'FLUXERR_AUTO', 'FLAGS',
            'FLUX RADIUS', 'CLASS_STAR', 'ALPHA_J2000', 'DELTA_J2000']
    df0 = pd.read_csv(os.path.join('new_fits', 'tmp_%s.cat' %
                                   band), delim_whitespace=True, names=name)
    hdu = fits.open(os.path.join('new_fits', field + '_new.ldac.fits'))

    print 'number of total stars (objects) found', df0.shape[0]
    # df0 = df0[(df0['FLUX_AUTO'] < float(param['flux_auto_%s'%band])).values]
    # df0=df0[(df0['FLAGS']==0)] #CHIPS1309-0406
    df0=df0[(df0['FLAGS']<5)]
    # df0=df0[(df0['FLAGS']<5).values & (df0['CLASS_STAR']>0.95).values] #edit: Field292, CHIPS2223-3455
    # df0=df0[(df0['FLAGS']<5).values & (df0['CLASS_STAR']>0.90).values] #edit: SDSS501
    # df0=df0[(df0['FLAGS']<2)&(df0['CLASS_STAR']>0.9)]

    print 'number of stars (CLASS_STAR>0.8 & FLAGS<5) using in Sextractor', len(np.array(df0.index))
    hdu[2].data = hdu[2].data[np.array(df0.index)]
    hdu.writeto(os.path.join('new_fits', field + '_new.ldac.fits'), overwrite=True)

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
    cmd = 'scamp %s -c pisco_pipeline/config.scamp' % ' '.join(
        list_file_name('new_fits', fieldname, end='_new.ldac.fits'))  # % band))
    print cmd
    subprocess.check_call(shlex.split(cmd))

def swarp(fieldname):
    """
    swarp: run SWARP to combine multiple exposure into a better image with SCAMP output to help correct the location
    INPUT:
    - config.swarp: SWARP config file
    - fieldname: begining of the file name (e.g., 'cField027')
    OUTPUT:
    - final/coadd_'fieldname'_'g'.fits: final image for each 'band' with corrected WCS
    """
    bands = ['g', 'r', 'i', 'z']
    print 'Swarping...'
    for band in bands:
        coadd_files = list_file_name(
            'new_fits', fieldname, end=band + '_new.fits')

        cmd = 'swarp %s -c pisco_pipeline/config.swarp -IMAGEOUT_NAME %s' %\
            (' '.join(coadd_files), os.path.join(
                'final', 'coadd_' + fieldname + '_' + band + '.fits'))
        print cmd
        subprocess.check_call(shlex.split(cmd))

def save_rgb_image(field):
    cmd = "ds9 -zscale -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom to fit -saveimage final/img%s.eps -exit" % \
        (field, field, field, field)  # -exit
    print cmd

##--------

if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    dir = str(sys.argv[1])
    fieldname = str(sys.argv[2])
    outdir = 'wcs'
    reducedir = 'reduced'
    cosmicdir = os.path.join(reducedir, 'cosmics')

    if len(sys.argv) > 3:
        flattype = str(sys.argv[3])
    else:
        flattype = 'domeflat'

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

    fields = [name.split('/')[-1].split('.')[0] for name in list_file_name(dir, fieldname)]
    print 'All fields:', fields

    for field in fields:
        for index in [1, 3, 5, 7]: #
            if not os.path.isfile(os.path.join(reducedir, "%s_%s.fits" % (field, filter_name(index)[0]))):
                ch1, bias1, domeflat1, img1 = reduce_data(dir, index, field, flat=flattype)
                ch2, bias2, domeflat2, img2 = reduce_data(dir, index + 1, field, flat=flattype)
                final_image = np.concatenate((img1, img2 * np.median(img1[:,750:]) / np.median(img2[:,:-100])), axis=1)
                save_fits(index, dir, reducedir, field, final_image, "%s_%s.fits" % (field, filter_name(index)[0]))
            else:
                print 'already have', "%s_%s.fits" % (field, filter_name(index)[0])
    
    bands = ['g', 'r', 'i', 'z']
    for field in fields:
        for band in bands:
            cosmic_reduce(reducedir, field, band)

    cfieldname = 'c' + fieldname
    print 'number of files in %s is %i' % (cosmicdir, len(list_file_name(cosmicdir, cfieldname)))

    for field_location in list_file_name(cosmicdir, cfieldname):
        myReg2=re.compile(r'(cCHIPS)\w+.\w+')
        field_all=myReg2.search(field_location).group()
        band = field_all[-1]

        astrometry_solve(cosmicdir, field_location, reducedir)
        sextracting(field_all, band, dir)

    scamp(cfieldname)
    swarp(cfieldname)
    save_rgb_image(fieldname)

    