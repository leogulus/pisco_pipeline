import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
import subprocess
import shlex
import sys
# edited 8/23/17
# ------
#v3: instead of using photoutil, we reproject the frame and use dual mode of sextractor to get the magnitudes
# also add seeing (FWHM) look up for sextractor

#v4: MAG_PSF and MAG_MODEL new photometry calibration

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
    seeing = float(fits.open(list_file_name_seeing(
        '/Users/taweewat/Documents/pisco_code/', field, startdir='ut')[0])[0].header['FWHM1'])

    slrdir = 'slr_output'
    to_be_projected = 'final/coadd_c%s_%s.fits'%(field,band)
    reference_fits  = 'final/coadd_c%s_i.fits'%field
    im1,im2, header = FITS_tools.match_fits(to_be_projected,reference_fits,return_header=True)
    outname = 'final/proj_coadd_c%s_%s.fits'%(field,band)
    print 'projecting from %s band to i band the fits file '%band + outname
    fits.writeto(outname, im1, header, overwrite=True)

    # cmd='sex final/coadd_c%s_i.fits,final/proj_coadd_c%s_%s.fits -c pisco_pipeline/config_slr.sex -CATALOG_NAME %s/mag_i%s.fits -SEEING_FWHM %s -SATUR_LEVEL %s -CHECKIMAGE_NAME %s/check_%s.fits' % \
    # (field,field,band,slrdir,band,str(seeing),str(param['satur_level_%s'%band]),slrdir,band)
    # print cmd
    # sub=subprocess.check_call(shlex.split(cmd))

    cmd='sex final/coadd_c%s_i.fits,final/proj_coadd_c%s_%s.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s'%\
    (field,field,band,'sex_psf.param','psfex_output/psf_%s_%s.fits'%(field,band),str(seeing),str(param['satur_level_%s'%band]))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    cmd='psfex %s -c pisco_pipeline/pisco.psfex' % ('psfex_output/psf_%s_%s.fits'%(field,band))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    cmd='sex final/coadd_c%s_i.fits -c pisco_pipeline/config.sex -PSF_NAME %s -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PIXEL_SCALE 0.2 -CATALOG_TYPE FITS_1.0 -CHECKIMAGE_NAME %s'%\
    (field,'psfex_output/psf_%s_%s.psf'%(field,band),'sex_after_psf.param','%s/a_psf_%s_%s.fits'%(slrdir,field,band),str(seeing),str(param['satur_level_%s'%band]),'check_%s_%s.fits'%(field,band))
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    table=Table.read('%s/a_psf_%s_%s.fits'%(slrdir,field,band))
    for name in table.colnames[1:]:
        table.rename_column(name, name + '_%s' % band)
    return table

def slr_running(field, mode, bigmacs="pisco_pipeline/big-macs-calibrate-master"):
    """
    slr_running: running SLR script from github.com/patkel/big-macs-calibrate to get a calibrated magnitude
    INPUT:
    - field: object of interset e.g., 'Field026'
    - bigmacs: the location for "big-macs-calibrate" directoty
    OUTPUT:
    - a new table with added columns with name MAG_g,...,MAGERR_g,...
    """
    slrdir = 'slr_output'
    infile = slrdir+'/star_%s.fits' % field
    pyfile = os.path.join(bigmacs, 'fit_locus.py')
    cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 5 -l -r ALPHA_J2000_i -d DELTA_J2000_i -j --plot=PLOTS_%s" \
        % (pyfile, infile, os.path.join(bigmacs, "coadd_mag_sex_%s.columns"%mode), field)
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

def update_color(fname, table, mode):
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
        table['MAG_cPSF_' + band[0]] = table['MAG_PSF_' + band[0]] + corr[0]
        table['MAG_cPSF_' + band[1]] = table['MAG_PSF_' + band[1]] + corr[1]
        table['MAG_cPSF_' + band[2]] = table['MAG_PSF_' + band[2]] + corr[2]
        table['MAG_cPSF_' + band[3]] = table['MAG_PSF_' + band[3]] + corr[3]
    elif mode=='model':
        table['MAG_cMODEL_' + band[0]] = table['MAG_MODEL_' + band[0]] + corr[0]
        table['MAG_cMODEL_' + band[1]] = table['MAG_MODEL_' + band[1]] + corr[1]
        table['MAG_cMODEL_' + band[2]] = table['MAG_MODEL_' + band[2]] + corr[2]
        table['MAG_cMODEL_' + band[3]] = table['MAG_MODEL_' + band[3]] + corr[3]
    elif mode=='auto':
        table['MAG_cAUTO_' + band[0]] = table['MAG_AUTO_' + band[0]] + corr[0]
        table['MAG_cAUTO_' + band[1]] = table['MAG_AUTO_' + band[1]] + corr[1]
        table['MAG_cAUTO_' + band[2]] = table['MAG_AUTO_' + band[2]] + corr[2]
        table['MAG_cAUTO_' + band[3]] = table['MAG_AUTO_' + band[3]] + corr[3]
    return table


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
    mag_ig=aperature_proj(field,'g')
    mag_ii=aperature_proj(field,'i')
    mag_ir=aperature_proj(field,'r')
    mag_iz=aperature_proj(field,'z')

    total=join(join(join(mag_ii,mag_ig,keys='NUMBER'), mag_ir,keys='NUMBER'),mag_iz,keys='NUMBER')
    # print total[1]
    total2=total[['ALPHA_J2000_i','DELTA_J2000_i',\
                  'MAG_AUTO_i','MAGERR_AUTO_i','MAG_AUTO_g','MAGERR_AUTO_g','MAG_AUTO_r','MAGERR_AUTO_r','MAG_AUTO_z','MAGERR_AUTO_z',\
                  'CLASS_STAR_i','CLASS_STAR_g','CLASS_STAR_r','CLASS_STAR_z',\
                  'FLAGS_g','FLAGS_r','FLAGS_i','FLAGS_z',\
                  'MAG_PSF_g','MAG_PSF_r','MAG_PSF_i','MAG_PSF_z','MAGERR_PSF_g','MAGERR_PSF_r','MAGERR_PSF_i','MAGERR_PSF_z',\
                  'MAG_MODEL_g','MAG_MODEL_r','MAG_MODEL_i','MAG_MODEL_z','MAGERR_MODEL_g','MAGERR_MODEL_r','MAGERR_MODEL_i','MAGERR_MODEL_z',\
                  'SPREAD_MODEL_g','SPREAD_MODEL_r','SPREAD_MODEL_i','SPREAD_MODEL_z']]

    # total3=total2[total2['CLASS_STAR_i'] > 0.90]
    total2.write(slrdir+'/all_psf_%s.fits' % field, overwrite=True)

    # slr_running(field,'auto')
    # ntotal = update_color(slrdir+'/star_%s.fits.offsets.list'%field, total2, 'auto')
    # ntotal.write(os.path.join(slrdir, 'n_auto_total_%s.csv' % field), overwrite=True)
    #
    # slr_running(field,'psf')
    # ntotal = update_color(slrdir+'/star_%s.fits.offsets.list'%field, total2, 'psf')
    # ntotal.write(os.path.join(slrdir, 'n_psf_total_%s.csv' % field), overwrite=True)
    #
    # slr_running(field,'model')
    # ntotal = update_color(slrdir+'/star_%s.fits.offsets.list'%field, total2, 'model')
    # ntotal.write(os.path.join(slrdir, 'n_model_total_%s.csv' % field), overwrite=True)



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


# ALPHA_J2000
# DELTA_J2000
