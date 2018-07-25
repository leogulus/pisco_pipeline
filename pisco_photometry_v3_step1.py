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
    myReg=re.compile(r'%s_A_\d{1,4}\.fits'%field)
    for root, dirs, files in os.walk('/Users/taweewat/Documents/pisco_code/'):
        for file in files:
            if myReg.search(file) != None:
                seeing=float(fits.open(root+'/'+myReg.search(file).group())[0].header['FWHM1'])

    # print seeing
    # seeing = float(fits.open(list_file_name_seeing(
    #     '/Users/taweewat/Documents/pisco_code/', field, startdir='ut')[0])[0].header['FWHM1'])
    # def see_px(see):
    #     return (0.24-0.1)/(1.4-0.6)*(see)

    slrdir = 'slr_output'
    to_be_projected = 'final/coadd_c%s_%s.fits'%(field,band)
    reference_fits  = 'final/coadd_c%s_i.fits'%field
    im1,im2, header = FITS_tools.match_fits(to_be_projected,reference_fits,return_header=True)
    outname = 'final/proj_coadd_c%s_%s.fits'%(field,band)
    print 'projecting from %s band to i band the fits file '%band + outname
    fits.writeto(outname, im1, header, overwrite=True)

    # cmd='sex final/coadd_c%s_i.fits,final/proj_coadd_c%s_%s.fits -c pisco_pipeline/config_slr.sex -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PIXEL_SCALE %s -CHECKIMAGE_NAME %s' % \
    # (field,field,band,"%s/mag_i%s.fits"%(slrdir,band),str(seeing),str(param['satur_level_%s'%band]),str(see_px(seeing)),"%s/check_%s.fits"%(slrdir,band))
    # print cmd
    # sub=subprocess.check_call(shlex.split(cmd))

    df_see=pd.read_csv('/Users/taweewat/Documents/red_sequence/total_chips_field_seeing.csv',index_col=0)
    if field[0:5]=='CHIPS':
        seeing = df_see[df_see.chips==field]['seeing_q25'].values[0] #np.min(df_see[df_see.chips==field][['seeing_025','seeing_gra_025']].values)
        print seeing
    elif (field[0:5]=='Field')|(field[0:3]=='PKS'):
        seeing = df_see[df_see.name==field]['seeing_q25'].values[0] #np.min(df_see[df_see.name==field][['seeing_025','seeing_gra_025']].values)
        print seeing

    # if field=='CHIPS1011-0505':
    #     seeing=0.95
    # if field=='Field179':
    #     seeing=1.12

    # if seeing <= 0.65:
    #     seeing=0.9
    # # elif seeing > 1.3:
    # #     seeing=1.34
    # elif seeing > 1.:
    #     seeing=seeing
    # else:
    #     seeing=1. #0.95 (1011), 1.0 (0005)

    # seeing=1.1
    # seeing=0.95
    minarea=1.7 #field159

    cmd='sex final/coadd_c%s_i.fits,final/proj_coadd_c%s_%s.fits -c pisco_pipeline/config_slr.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE 0.22 -DETECT_MINAREA %s -CHECKIMAGE_NAME checki.fits'%\
    (field,field,band,'sex_slr.param',"%s/mag_i%s.fits"%(slrdir,band),str(seeing),str(param['satur_level_%s'%band]),str(1.1/1.7*np.pi*(seeing/0.22)**2)); print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    table=Table.read(slrdir+'/mag_i%s.fits'%band)

    for name in table.colnames[1:]:
        table.rename_column(name, name + '_%s' % band)
    return table

def slr_running(field, bigmacs="pisco_pipeline/big-macs-calibrate-master"):
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
    # infile = slrdir+'/star_bleem_%s.fits' % field
    pyfile = os.path.join(bigmacs, 'fit_locus.py')
    cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 5 -l -r ALPHA_J2000_i -d DELTA_J2000_i -j --plot=PLOTS_%s" \
        % (pyfile, infile, os.path.join(bigmacs, "coadd_mag_sex.columns"), field)
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

def update_color(fname, table):
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

    table['MAG_' + band[0]] = table['MAG_AUTO_' + band[0]] + corr[0]
    table['MAG_' + band[1]] = table['MAG_AUTO_' + band[1]] + corr[1]
    table['MAG_' + band[2]] = table['MAG_AUTO_' + band[2]] + corr[2]
    table['MAG_' + band[3]] = table['MAG_AUTO_' + band[3]] + corr[3]
    table['MAGERR_' + band[0]] = table['MAGERR_AUTO_' + band[0]] + ecorr[0]
    table['MAGERR_' + band[1]] = table['MAGERR_AUTO_' + band[1]] + ecorr[1]
    table['MAGERR_' + band[2]] = table['MAGERR_AUTO_' + band[2]] + ecorr[2]
    table['MAGERR_' + band[3]] = table['MAGERR_AUTO_' + band[3]] + ecorr[3]
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
python pisco_pipeline/pisco_photometry_v3_step1.py SDSS123
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

    total.write(os.path.join(slrdir, 'total_%s.csv' % field), overwrite=True)

    # total2=total[['NUMBER','ALPHA_J2000_i','DELTA_J2000_i','MAG_AUTO_i','MAGERR_AUTO_i','MAG_AUTO_g','MAGERR_AUTO_g',\
    #               'MAG_AUTO_r','MAGERR_AUTO_r','MAG_AUTO_z','MAGERR_AUTO_z','CLASS_STAR_i','CLASS_STAR_g',\
    #               'CLASS_STAR_r','CLASS_STAR_z','FLAGS_i','FLAGS_g','FLAGS_r','FLAGS_z']]
    # total2=total2[(total2['FLAGS_g']<5)&(total2['FLAGS_r']<5)&(total2['FLAGS_i']<5)&(total2['FLAGS_z']<5)]
    # total3=total2[total2['CLASS_STAR_i'] > 0.9]
    # print len(total3)
    # total3.write(slrdir+'/star_%s.fits' % field, overwrite=True)
    #
    # slr_running(field)
    # # ntotal = update_color(slrdir+'/star_bleem_%s.fits.offsets.list'%field, total)
    # # ntotal.write(os.path.join(slrdir, 'ntotal_bleem_%s.csv' % field), overwrite=True)
    # ntotal = update_color(slrdir+'/star_%s.fits.offsets.list'%field, total)
    # ntotal.write(os.path.join(slrdir, 'ntotal_%s.csv' % field), overwrite=True)
    #
    # purge('final', "proj_coadd_c%s_.*\.fits" % field)

    print 'test'
