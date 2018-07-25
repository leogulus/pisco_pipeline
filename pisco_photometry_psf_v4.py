import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from photutils import aperture_photometry
from photutils import SkyCircularAperture
from photutils import SkyCircularAnnulus

import FITS_tools
import yaml

import os
import re
import subprocess
import shlex
import sys
# edited 8/23/17
# ------
#v3: instead of using photoutil, we reproject the frame and use dual mode of sextractor to get the magnitudes
# also add seeing (FWHM) look up for sextractor

def list_file_name_seeing(dir, name, end=0, startdir=0):
    names=[]
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith(name):
                if end == 0:
                    if startdir == 0:
                        names.append(os.path.join(root, file))
                    else:
                        if root.split('/')[-1][:2]==startdir:
                            names.append(os.path.join(root, file))
                else:
                    if file.endswith(end):
                        if startdir == 0:
                            names.append(os.path.join(root, file))
                        else:
                            if root.split('/')[-1][:2]==startdir:
                                names.append(os.path.join(root, file))
    if len(names) == 0:
        print 'Cannot find the files'
    return names

def aperature_proj(field,band):
    with open("pisco_pipeline/params.yaml", 'r') as stream:
        try:
            param=yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # see=float(fits.open(list_file_name_seeing('/Users/taweewat/Documents/pisco_code/',field,startdir='ut')[0])[0].header['FWHM1'])
    # seeing = float(fits.open(list_file_name_seeing(
    #     '/Users/taweewat/Documents/pisco_code/', field, startdir='ut')[0])[0].header['FWHM1'])
    # def see_px(see):
    #     return (0.24-0.1)/(1.4-0.6)*(see)
    df_see=pd.read_csv('/Users/taweewat/Documents/red_sequence/total_chips_field_seeing.csv',index_col=0)
    if field[0:5]=='CHIPS':
        seeing = df_see[df_see.chips==field]['seeing_q25'].values[0] #np.min(df_see[df_see.chips==field][['seeing_025','seeing_gra_025']].values)
        print seeing
    elif (field[0:5]=='Field')|(field[0:3]=='PKS'):
        seeing = df_see[df_see.name==field]['seeing_q25'].values[0] #np.min(df_see[df_see.name==field][['seeing_025','seeing_gra_025']].values)
        print seeing

    slrdir = 'slr_output'
    to_be_projected = 'final/coadd_c%s_%s.fits'%(field,band)
    reference_fits  = 'final/coadd_c%s_i.fits'%field
    im1,im2, header = FITS_tools.match_fits(to_be_projected,reference_fits,return_header=True)
    outname = 'final/proj_coadd_c%s_%s.fits'%(field,band)
    print 'projecting from %s band to i band the fits file '%band + outname
    fits.writeto(outname, im1, header, overwrite=True)

    cmd='sex final/coadd_c%s_i.fits,final/proj_coadd_c%s_%s.fits -c pisco_pipeline/config_slr.sex -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PIXEL_SCALE %s -CHECKIMAGE_NAME %s' % \
    (field,field,band,"%s/mag_i%s.fits"%(slrdir,band),str(seeing),str(param['satur_level_%s'%band]),str(see_px(seeing)),"%s/check_%s.fits"%(slrdir,band))
    print cmd
    sub=subprocess.check_call(shlex.split(cmd))

    table=Table.read(slrdir+'/mag_i%s.fits'%band)

    for name in table.colnames[1:]:
        table.rename_column(name, name + '_%s' % band)
    return table

def slr_running_psf(field, infile="None", mode="psf", bigmacs="pisco_pipeline/big-macs-calibrate-master"):
    """
    slr_running: running SLR script from github.com/patkel/big-macs-calibrate to get a calibrated magnitude
    INPUT:
    - field: object of interset e.g., 'Field026'
    - bigmacs: the location for "big-macs-calibrate" directoty
    OUTPUT:
    - a new table with added columns with name MAG_g,...,MAGERR_g,...
    """
    slrdir = 'slr_output'
    # infile = slrdir+'/star_psf_%s_%s.fits' % (mode,field)
    pyfile = os.path.join(bigmacs, 'fit_locus.py')
    cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 15 -l -r ALPHA_J2000_i -d DELTA_J2000_i -j --plot=PLOTS_%s_%s" \
        % (pyfile, infile, os.path.join(bigmacs, "coadd_mag_sex_%s.columns"%mode), mode, field)
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

def update_color(fname, table, mode='psf'):
    """
    update_color: using the output from SLR, update to the correct magnitude
    INPUT:
    - fname: input file from SLR output (...offsets.list)
    - table: the table that we want to update the value (from column magg,etc to MAG_g,etc)
    OUTPUT:
    - a new table with added columns with name MAG_g,...,MAGERR_g,...
    """
    with open(fname) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    band = [x.split(' ')[0][-1] for x in content[5:-1]]
    corr = [float(x.split(' ')[1]) for x in content[5:-1]]
    ecorr = [float(x.split(' ')[3]) for x in content[5:-1]]
    print 'bands = ', band

    if mode=='psf':
        MODE1='PSF'
    elif mode=='model':
        MODE1='MODEL'
    elif mode=='auto':
        MODE1='AUTO'

    # table['MAG1_%s_'%MODE1 + band[0]] = table['MAG_%s_'%MODE1 + band[0]] + corr[0]
    # table['MAG1_%s_'%MODE1 + band[1]] = table['MAG_%s_'%MODE1 + band[1]] + corr[1]
    # table['MAG1_%s_'%MODE1 + band[2]] = table['MAG_%s_'%MODE1 + band[2]] + corr[2]
    # table['MAG1_%s_'%MODE1 + band[3]] = table['MAG_%s_'%MODE1 + band[3]] + corr[3]
    # table['MAGERR1_%s_'%MODE1 + band[0]] = table['MAGERR_%s_'%MODE1 + band[0]] + ecorr[0]
    # table['MAGERR1_%s_'%MODE1 + band[1]] = table['MAGERR_%s_'%MODE1 + band[1]] + ecorr[1]
    # table['MAGERR1_%s_'%MODE1 + band[2]] = table['MAGERR_%s_'%MODE1 + band[2]] + ecorr[2]
    # table['MAGERR1_%s_'%MODE1 + band[3]] = table['MAGERR_%s_'%MODE1 + band[3]] + ecorr[3]
    table['MAG_' + band[0]] = table['MAG_%s_'%MODE1 + band[0]] + corr[0]
    table['MAG_' + band[1]] = table['MAG_%s_'%MODE1 + band[1]] + corr[1]
    table['MAG_' + band[2]] = table['MAG_%s_'%MODE1 + band[2]] + corr[2]
    table['MAG_' + band[3]] = table['MAG_%s_'%MODE1 + band[3]] + corr[3]
    table['MAGERR_' + band[0]] = (table['MAGERR_%s_'%MODE1 + band[0]]**2 + ecorr[0]**2)**0.5
    table['MAGERR_' + band[1]] = (table['MAGERR_%s_'%MODE1 + band[1]]**2 + ecorr[1]**2)**0.5
    table['MAGERR_' + band[2]] = (table['MAGERR_%s_'%MODE1 + band[2]]**2 + ecorr[2]**2)**0.5
    table['MAGERR_' + band[3]] = (table['MAGERR_%s_'%MODE1 + band[3]]**2 + ecorr[3]**2)**0.5
    return table

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))

"""
pisco_photometry: run pisco output data from pisco_combine to correct for the photometry of each object
and determine which objects are stars/galaxies.
The pipeline is a combination of SLR algorithm (cite: https://github.com/patkel/big-macs-calibrate)
and Photutils for photometry aperatures

ARGUMENTS:
1. fieldname for object (e.g., 'Field027')

EXAMPLES:
python pisco_pipeline/pisco_photometry_v3.py SDSS123
"""

if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    slrdir = 'slr_output'
    if not os.path.exists(slrdir):
        os.makedirs(slrdir)

    field = str(sys.argv[1])

    total3=Table.from_pandas(pd.read_csv("/Users/taweewat/Documents/pisco_code/slr_output/star_psf_total_%s.csv"%field))
    total3=total3[['NUMBER','ALPHA_J2000_i','DELTA_J2000_i',\
                  'MAG_AUTO_i','MAGERR_AUTO_i','MAG_AUTO_g','MAGERR_AUTO_g','MAG_AUTO_r','MAGERR_AUTO_r','MAG_AUTO_z','MAGERR_AUTO_z',\
                  'CLASS_STAR_i','CLASS_STAR_g','CLASS_STAR_r','CLASS_STAR_z',\
                  'FLAGS_g','FLAGS_r','FLAGS_i','FLAGS_z',\
                  'MAG_PSF_g','MAG_PSF_r','MAG_PSF_i','MAG_PSF_z','MAGERR_PSF_g','MAGERR_PSF_r','MAGERR_PSF_i','MAGERR_PSF_z',\
                  'MAG_MODEL_g','MAG_MODEL_r','MAG_MODEL_i','MAG_MODEL_z','MAGERR_MODEL_g','MAGERR_MODEL_r','MAGERR_MODEL_i','MAGERR_MODEL_z',\
                  'SPREAD_MODEL_g','SPREAD_MODEL_r','SPREAD_MODEL_i','SPREAD_MODEL_z']]

    print 'number of stars =', len(total3)
    total3.write(slrdir+'/star_psf_psf_%s_%i.fits' % (field,0), overwrite=True)

    # dff=dff[(dff.FLAGS_i<5)&(dff.FLAGS_g<5)&(dff.FLAGS_r<5)&(dff.FLAGS_z<5)]  #remove bad data
    # dff=dff[(dff.MAG_PSF_i<0)&(dff.MAG_MODEL_i<0)&(dff.MAG_AUTO_i<0)&(dff.MAG_PSF_r<0)&(dff.MAG_MODEL_r<0)&(dff.MAG_AUTO_r<0)&(dff.MAG_PSF_g<0)&(dff.MAG_MODEL_g<0)&(dff.MAG_AUTO_g<0)]

    #only select stars for SLR
    # dff1=dff[np.abs(dff.SPREAD_MODEL_i)<0.005]
    # dff=dff[dff.CLASS_STAR_i > 0.9]

    # filename='/Users/taweewat/Documents/pisco_code/slr_output/all_psf_%s.fits'%field
    # filename='/Users/taweewat/Documents/pisco_code/slr_output/allstar_psf_%s.fits'%field
    # print filename
    # dff=Table(fits.open(filename)[1].data).to_pandas()
    #

    # Table.from_pandas(dff).write(slrdir+'/star_psf_psf_%s_%i.fits' % (field,0), overwrite=True)
    # Table.from_pandas(dff).write(slrdir+'/star_psf_model_%s.fits' % field, overwrite=True)
    # Table.from_pandas(dff).write(slrdir+'/star_psf_auto_%s.fits' % field, overwrite=True)

    slr_running_psf(field, infile=slrdir+'/star_psf_psf_%s_%i.fits' % (field,0), mode='psf')

    total_gal=Table.from_pandas(pd.read_csv("/Users/taweewat/Documents/pisco_code/slr_output/galaxy_psf_total_%s.csv"%field))
    ntotal_gal = update_color(slrdir+'/star_psf_psf_%s_%i.fits.offsets.list'%(field,0), total_gal, mode='model')
    ntotal_gal.write(os.path.join(slrdir, 'galaxy_psf_ntotal_%s.csv' % field), overwrite=True)

    # ntotal = update_color(slrdir+'/star_psf_psf_%s_%i.fits.offsets.list'%(field,0), Table.from_pandas(dff), 'psf')
    # ntotal.write(os.path.join(slrdir, 'ntotal_psf_%s.csv' % field), overwrite=True)

    # slr_running_psf(field, 'model')
    # ntotal = update_color(slrdir+'/star_psf_model_%s.fits.offsets.list'%field, Table.from_pandas(dff), 'model')
    # ntotal.write(os.path.join(slrdir, 'ntotal_model_%s.csv' % field), overwrite=True)
    #
    # slr_running_psf(field, 'auto')
    # ntotal = update_color(slrdir+'/star_psf_auto_%s.fits.offsets.list'%field, Table.from_pandas(dff), 'auto')
    # ntotal.write(os.path.join(slrdir, 'ntotal_auto_%s.csv' % field), overwrite=True)

    # purge('final', "proj_coadd_c%s_.*\.fits" % field)



    # cmd = 'mv star_%s.fits slr_output/' % field
    # try:
    #     sub = subprocess.check_call(shlex.split(cmd))
    # except (ValueError, RuntimeError, TypeError, NameError):
    #     pass
    #
    # cmd = 'mv star_%s.fits.offsets.list slr_output/' % field
    # try:
    #     sub = subprocess.check_call(shlex.split(cmd))
    # except (ValueError, RuntimeError, TypeError, NameError):
    #     pass
