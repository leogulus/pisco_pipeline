import numpy as np
import matplotlib
import os
import re
import pandas as pd

from astropy.io import fits
import subprocess
import cosmics
import shlex
import sys

import yaml

import matplotlib.pyplot as plt
#from pisco_lib import *
# edited 5/9/17

"""
pisco_combine: run pisco pipeline to reduce the raw data to clean data with correct WCS
The pipeline is a combination of LA Cosmics, Astrometry, Sextractor, SCAMP and SWARP.

ARGUMENTS:
1. raw directory (e.g., 'ut170103/')
2. fieldname for object (e.g., 'Field027')

Examples: python pisco_pipeline/pisco_combine.py data/ Field026
python pisco_pipeline/pisco_combine.py ut170619/ Field022 'twilight'
python pisco_pipeline/pisco_combine.py ut171209/ CHIPS0152-5028 'twilight'
"""

# add twilight option to use twilight flats
# add look up seeing (FWHM) from the header
# scamp_v2 for doing scamp for g-r and i-z simultaniously, instead of doing all 4 bands simultaniously

#edit 2/13/17: added yaml parameter to control all the sextractor parameters

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


def list_file_name(dir, name, end=0):
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
    names = []
    for file in os.listdir(dir):
        if file.startswith(name):
            if end == 0:
                names.append(os.path.join(dir, file))
            else:
                if file.endswith(end):
                    names.append(os.path.join(dir, file))
    if len(names) == 0:
        print 'Cannot find the files'
    return names


def list_file_name_seeing(dir, name, end=0, startdir=0):
    names = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith(name):
                if end == 0:
                    if startdir == 0:
                        names.append(os.path.join(root, file))
                    else:
                        if root.split('/')[-1][:2] == startdir:
                            names.append(os.path.join(root, file))
                else:
                    if file.endswith(end):
                        if startdir == 0:
                            names.append(os.path.join(root, file))
                        else:
                            if root.split('/')[-1][:2] == startdir:
                                names.append(os.path.join(root, file))
    if len(names) == 0:
        print 'Cannot find the files'
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
            # for domeflat-bias before combine into a mean
            ch_bs.append(ch_b - bias)
    if twilight == True:
        print 'working on twlight flat'
        return np.median(np.array(ch_bs), axis=0)
    else:
        return np.mean(np.array(ch_bs), axis=0)


def plot_one_chip(ax, data, vmin, vmax):
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
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    c_m = matplotlib.cm.Greys
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    ax.imshow(data, cmap=c_m, norm=norm)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)


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
    hdulist[index].header['DATASEC'] = '[1:1546,1:3092]'
    hdulist[index].header['TRIMSEC'] = '[1:1546,1:3092]'
    hdulist[index].header['ORIGSEC'] = '[1:1546,1:3092]'
    hdulist[index].header['CCDSEC'] = '[1:1546,3093:6184]'
    hdulist[index].header['DETSEC'] = '[1:1546,3093:6184]'

    hdu1 = fits.ImageHDU(final_image * 1000, name='filter ' +
                         filter_name(index)[0], header=hdulist[index].header)
    hdu_l = fits.HDUList(hdus=[hdu0, hdu1])
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)
    outname = os.path.join(outdir, name)
    print 'saving the fits file ' + outname
    hdu_l.writeto(outname, overwrite=True)

    data, header = fits.getdata(outname, header=True)
    fits.writeto(outname, data, header, overwrite=True)


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
    print fieldname[0:5]
    if (fieldname[0:5] == 'Field') or (fieldname[0:3] == 'PKS'):
        cut = -27
    elif fieldname[0:5] == 'CHIPS':
        cut = -32
    else:
        print 'there is no fieldname for that' + fieldname
    print cut
    # if index == 3:
    #     edge=20
    # elif index == 6:
    #     edge=30
    # else:
    #     edge=10

    ch1_name = list_file_name(dir, fieldname)
    print 'working on %s with the index=%i' % (ch1_name[0], index)
    hdulist = fits.open(ch1_name[0])
    ch1 = hdulist[index].data

    bias_names = list_file_name(dir, 'Bias_')
    # print bias_names
    if flat == 'domeflat':
        domeflat_names = list_file_name(
            dir, "domeflat" + filter_name(index)[1])

    if flat == 'twilight':
        if (fieldname[0:5] == 'Field') or (fieldname[0:3] == 'PKS'):
            domeflat_names = list_file_name(dir, "twiflat_")
        elif fieldname[0:5] == 'CHIPS':
            domeflat_names = list_file_name(dir, "twiflats_")
        else:
            print 'there is no fieldname for twilight flat' + fieldname
        # domeflat_names = list_file_name(dir, "twiflat_")

    bias = open_files(bias_names, index)
    if flat == 'domeflat':
        domeflat = open_files(domeflat_names, index, bias=bias, twilight=False)
    if flat == 'twilight':
        domeflat = open_files(domeflat_names, index, bias=bias, twilight=True)

    domeflat[domeflat == 0] = 1e-4
    # domeflat = np.clip(domeflat, 1e-4, 1e10)

    # if index in [1,2,3,4]:
    #     mean=np.median(domeflat[350:2550, 10:-10])
    # elif index in [5,6,7,8]:
    #     mean=np.median(domeflat[650:2800, 10:-10])
    # domeflat=domeflat/mean

    print ch1.shape, bias.shape, domeflat.shape

    img = (ch1 - bias) / domeflat
    ch1, bias, domeflat, img = ch1[:, :cut], bias[:,
                                                  :cut], domeflat[:, :cut], img[:, :cut]

    if index % 2 == 0:
        return np.fliplr(ch1), np.fliplr(bias), np.fliplr(domeflat), np.fliplr(img)
    else:
        return ch1, bias, domeflat, img


def cosmic_reduce(dir, field, band):
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
    if not os.path.isfile(os.path.join(dir, 'cosmics', 'c' + field + '_' + band + '.fits')):
        print 'working on the cosmic ' + 'c' + field + '_' + band
        array, header = cosmics.fromfits(
            os.path.join(dir, field + '_' + band + '.fits'))
        print os.path.join(dir, field + '_' + band + '.fits')
        # cutting the circular aperature of the image out to only have good pixels
        # in the center

        with open("pisco_pipeline/params.yaml", 'r') as stream:
            try:
                param=yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        param['satur_level_%s'%band]

        if band == 'g':
            array_c = array[20:-20, 350:2550]  # [20:-20,350:2550]
            satlv = param['satur_level_%s'%band] #2000.0
        elif band == 'r':
            array_c = array[30:-20, 400:2550]  # [20:-20,350:2550]
            satlv = param['satur_level_%s'%band]# 1250.0
        elif band == 'i':
            array_c = array[20:-40, 650:2800]  # [20:-20,650:2800]
            satlv = param['satur_level_%s'%band] #600.0
        elif band == 'z':
            array_c = array[20:-20, 650:2750]  # [20:-20,650:2800]
            satlv = param['satur_level_%s'%band] #1500.0
        c = cosmics.cosmicsimage(array_c, gain=4.0, readnoise=3.0, sigclip=2.5, sigfrac=0.5,
                                 objlim=5.0, satlevel=satlv, verbose=False)
        c.run(maxiter=5)
        cosmics.tofits(os.path.join(dir, 'cosmics', 'c' + field +
                                    '_' + band + '.fits'), c.cleanarray, header)
    else:
        print 'already did the cosmic with this band ' + band

    # cosmics.tofits(os.path.join(dir, 'cosmics', 'n' + field +
    #                             '_' + band + '.fits'), array_c, header)
    # cosmics.tofits(os.path.join(dir, 'cosmics', 'm' + field +
    #                             '_' + band + '.fits'), c.mask, header)

def ra_dec(name):
    ra = float(name[6:8]) * 15 + float(name[8:10]) / 4.
    if name[10] == '-':
        dec = float(name[10:13]) - float(name[13:15]) / 60.
    else:
        dec = float(name[10:13]) + float(name[13:15]) / 60.
    return ra, dec

def ra_dec_field(field):
    de = pd.read_csv(
        '/Users/taweewat/Documents/red_sequence/field_chips_all_obj.csv')
    ra = de[de.name == field[1:9]].RA0.values[0]
    dec = de[de.name == field[1:9]].DEC0.values[0]
    return ra, dec


def astrometry_solve(cosmicdir, field, outdir):
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
    if not os.path.isfile(os.path.join(outdir, field + '.wcs')):
        print field, field[0:6]

        if field[0:6] == 'cCHIPS':
            ra, dec = ra_dec(field)
            # cmd = 'solve-field %s --downsample 2 --overwrite --scale-unit arcsecperpix --scale-low 0.08 --scale-high 0.3 --dir %s --ra %s --dec %s --radius 2' \
            #     % (os.path.join(cosmicdir, field + '.fits'), outdir, str(ra), str(dec))
        if field[0:4] == 'cPKS':
            ra=209.0225
            dec=-34.3530556
        else:
            ra, dec = ra_dec_field(field)
            # cmd = 'solve-field %s --downsample 2 --overwrite --scale-unit arcsecperpix --scale-low 0.08 --scale-high 0.3 --dir %s' \
            #     % (os.path.join(cosmicdir, field + '.fits'), outdir)
        cmd = 'solve-field %s --downsample 2 --overwrite --scale-unit arcsecperpix --scale-low 0.08 --scale-high 0.3 --dir %s --ra %s --dec %s --radius 2' \
            % (os.path.join(cosmicdir, field + '.fits'), outdir, str(ra), str(dec))
        print cmd
        sub = subprocess.check_call(shlex.split(cmd))
        if sub == 0:
            print 'finish solve-field and updating fits headers'
        else:
            print 'solve-field does not work.'
    else:
        print 'already have ' + field + '.wcs'

    orig = fits.open(os.path.join(cosmicdir, field + '.fits'))
    wcs_file = fits.open(os.path.join(outdir, field + '.wcs'))
    header = wcs_file[0].header
    wcsaxes_index = np.where(np.array(header.keys()) == 'WCSAXES')[0][0]
    for i in range(wcsaxes_index, len(header)):
        orig[0].header[header.keys()[i]] = header.values()[i]
    orig.writeto(os.path.join('new_fits', field + '_new.fits'), overwrite=True)


def sextracting(field, band):
    """
    sextracting: run Sextractor to find all the point sources in .ldac.fits format (suitable for SCAMP input)
    INPUT:
    - config.sex: sextractor config file
    - field: begining of the file name for each band and each exposure (e.g. 'cSDSS123_B_64_r')
    OUTPUT:
    - new_fits/..._new.ldac.fits: source catalogs of all the point source from Sextractor
    """

    with open("pisco_pipeline/params.yaml", 'r') as stream:
        try:
            param=yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    fieldname = field.split('_')[0][1:] + '_' + field.split('_')[1]
    seeing = float(fits.open(list_file_name_seeing(
        '/Users/taweewat/Documents/pisco_code/', fieldname, startdir='ut')[0])[0].header['FWHM1'])

    cmd = 'sex %s -c pisco_pipeline/config.sex -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -CHECKIMAGE_NAME %s,%s -PIXEL_SCALE 0.2' % \
        (os.path.join('new_fits', field + '_new.fits'),
         os.path.join('new_fits', field + '_new.ldac.fits'), str(seeing), str(param['satur_level_%s'%band]),'check_%s.fits'%(band),'segment_%s.fits'%(band))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    cmd = 'sex %s -c pisco_pipeline/config.sex -CATALOG_NAME %s -CATALOG_TYPE ASCII -SEEING_FWHM %s -SATUR_LEVEL %s -CHECKIMAGE_NAME %s,%s -PIXEL_SCALE 0.2' % \
        (os.path.join('new_fits', field + '_new.fits'),
         os.path.join('new_fits', 'tmp_%s.cat' % band), str(seeing), str(param['satur_level_%s'%band]),'check_%s.fits'%(band),'segment_%s.fits'%(band))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    name = ['NUMBER', 'EXT_NUMBER', 'XWIN_WORLD', 'YWIN_WORLD', 'MAG_AUTO', 'MAGERR_AUTO', 'MAG_APER', 'MAGERR_APER', 'XWIN_IMAGE',
            'YWIN_IMAGE', 'ERRAWIN_IMAGE', 'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE', 'FLUX_AUTO', 'FLUXERR_AUTO', 'FLAGS',
            'FLUX RADIUS', 'CLASS_STAR', 'ALPHA_J2000', 'DELTA_J2000']
    df0 = pd.read_csv(os.path.join('new_fits', 'tmp_%s.cat' %
                                   band), delim_whitespace=True, names=name)
    hdu = fits.open(os.path.join('new_fits', field + '_new.ldac.fits'))

    print 'number of total stars (objects) found', df0.shape[0]
    df0 = df0[(df0['FLUX_AUTO'] < float(param['flux_auto_%s'%band])).values]
    # if band == 'g':
    #     df0 = df0[(df0['FLUX_AUTO'] < 7e6).values]  # 7e6 3e5
    # if band == 'r':
    #     df0 = df0[(df0['FLUX_AUTO'] < 8e6).values]  # 8e6 3e5
    # if band == 'i':
    #     df0 = df0[(df0['FLUX_AUTO'] < 6e5).values]  # 6e5 3e5
    # if band == 'z':
    #     df0 = df0[(df0['FLUX_AUTO'] < 7e6).values]  # 7e6 3e5

    # if band=='i':
    #     df0=df0.sort_values(by='FLUX_AUTO').iloc[:int(df0.shape[0]/1.09)]
    # if band=='r':
    #     df0=df0.sort_values(by='FLUX_AUTO').iloc[:int(df0.shape[0]/1.08)]
    # if band=='g':
    #     df0=df0.sort_values(by='FLUX_AUTO').iloc[:int(df0.shape[0]/1.09)]
    # if band=='z':
    #     df0=df0.sort_values(by='FLUX_AUTO').iloc[:int(df0.shape[0]/1.09)]

    df0=df0[(df0['FLAGS']<5).values]  #(df0['CLASS_STAR']>0.8).values &

    print 'number of stars (CLASS_STAR>0.8 & FLAGS<5) using in Sextractor', len(np.array(df0.index))
    hdu[2].data = hdu[2].data[np.array(df0.index)]
    hdu.writeto(os.path.join('new_fits', field +
                             '_new.ldac.fits'), overwrite=True)

    # if band == 'i':
    # ploting PNG for viewing purpose
    # img = fits.open('segment_%s.fits' % band)[0].data
    # cmap = plt.cm.spring
    # cmap.set_bad(color='black')
    # img0 = np.ma.masked_where(img < 0.05, img)
    # plt.imshow(img0, origin='lower', cmap=cmap, interpolation='none')
    # plt.savefig('plt_bf_%s_%s.png' % (band, str(np.random.randint(100))))
    #
    # img[~np.in1d(img, df0['NUMBER']).reshape(img.shape)] = 0
    # img = np.ma.masked_where(img < 0.05, img)
    # plt.imshow(img, origin='lower', cmap=cmap, interpolation='none')
    # plt.savefig('plt_%s_%s.png' % (band, str(np.random.randint(100))))


def scamp(fieldname):  # , band):
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
    sub = subprocess.check_call(shlex.split(cmd))


def scamp_v2(fieldname):  # , band):
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
        list_file_name('new_fits', fieldname, end='g_new.ldac.fits') + list_file_name('new_fits', fieldname, end='r_new.ldac.fits'))  # % band))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    cmd = 'scamp %s -c pisco_pipeline/config.scamp' % ' '.join(
        list_file_name('new_fits', fieldname, end='i_new.ldac.fits') + list_file_name('new_fits', fieldname, end='z_new.ldac.fits'))  # % band))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))


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
        sub = subprocess.check_call(shlex.split(cmd))


def save_rgb_image(field):
    cmd = "ds9 -zscale -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom to fit -saveimage final/img%s.eps -exit" % \
        (field, field, field, field)  # -exit
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))
    print 'finished saving final/img%s.eps' % field


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))


def run_bash_code(cmd):
    print cmd
    try:
        sub = subprocess.check_output(shlex.split(cmd))
    except subprocess.CalledProcessError as e:
        print e

# --------


if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    # Pipeline to run PISCO reduction data
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

    fields = [name.split('/')[-1].split('.')[0]
              for name in list_file_name(dir, fieldname)]

    # Combine two amplifiers and bias and flat fielding
    for field in fields:
        for index in [1, 3, 5, 7]:
            ch1, bias1, domeflat1, img1 = reduce_data(
                dir, index, field, flat=flattype)
            ch2, bias2, domeflat2, img2 = reduce_data(
                dir, index + 1, field, flat=flattype)
            final_image = np.concatenate((img1, img2), axis=1)
            save_fits(index, dir, reducedir, field, final_image,
                      "%s_%s.fits" % (field, filter_name(index)[0]))

    # Cosmic ray reduction using L.A. Cosmic
    bands = ['g', 'r', 'i', 'z']
    for field in fields:
        for band in bands:
            # if not os.path.isfile(os.path.join(reducedir, 'cosmics', 'c' + field + '_' + band + '.fits')):
            #     print 'working on the cosmic ' + 'c' + field + '_' + band
            cosmic_reduce(reducedir, field, band)
            # else:
            #     print 'already did the cosmic with this band ' + band
            #     # cosmic_reduce(reducedir,field,band)

    cfieldname = 'c' + fieldname
    print 'number of files in %s is %i' % (cosmicdir, len(list_file_name(cosmicdir, cfieldname)))

    os.system('rm plt*')

    for field_long in list_file_name(cosmicdir, cfieldname):
        field = field_long.split('/')[2].split('.')[0]
        print 'Field', field
        # Astrometry to get a rough estimate on the World Coordinate System
        # (WCS) for each images
        astrometry_solve(cosmicdir, field, outdir)
        # Sextracting
        band = field.split('_')[3]
        print band
        sextracting(field, band[0])

    if fieldname == 'Field173':
        cmd = 'rm new_fits/cField173_A_100_i_new.fits'
        run_bash_code(cmd)

    # # SCAMP
    # # for band in bands:
    # scamp_v2(cfieldname)
    scamp(cfieldname)
    # # SWARP
    swarp(cfieldname)
    # # save eps file for RGB image
    save_rgb_image(fieldname)

    purge('reduced', "%s_.*\.fits" % fieldname)
    purge('new_fits', "c%s_.*\.fits" % fieldname)
    purge('wcs', "c%s_.*\.new" % fieldname)
