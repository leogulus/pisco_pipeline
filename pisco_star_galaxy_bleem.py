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

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))

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



"""
pisco_star_galaxy_bleem:

ARGUMENTS:
1. fieldname for object (e.g., 'Field027')

EXAMPLES:
python pisco_pipeline/pisco_star_galaxy_bleem.py CHIPS0137-1248
"""

if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    sg_dir = 'star_galaxy'
    if not os.path.exists(sg_dir):
        os.makedirs(sg_dir)

    field = str(sys.argv[1])

    with open("pisco_pipeline/params.yaml", 'r') as stream:
        try:
            param=yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    seeing = float(fits.open(list_file_name_seeing(
        '/Users/taweewat/Documents/pisco_code/', field, startdir='ut')[0])[0].header['FWHM1'])

    df_see=pd.read_csv('/Users/taweewat/Documents/red_sequence/total_chips_field_seeing.csv',index_col=0)
    if field[0:5]=='CHIPS':
        seeing = df_see[df_see.chips==field]['seeing_q25'].values[0] #np.min(df_see[df_see.chips==field][['seeing_025','seeing_gra_025']].values)
        print seeing
    elif (field[0:5]=='Field')|(field[0:3]=='PKS'):
        seeing = df_see[df_see.name==field]['seeing_q25'].values[0] #np.min(df_see[df_see.name==field][['seeing_025','seeing_gra_025']].values)
        print seeing

    # seeing=0.95
    minarea=1.7

    data, header = fits.getdata('final/coadd_c%s_i.fits'%field, header=True)
    data2=data**2
    fits.writeto('final/coadd_c%s_sq_i.fits'%field, data2, header=header, overwrite=True)

    cmd='sex final/coadd_c%s_i.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -CATALOG_TYPE FITS_1.0 -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE 0.22 -DETECT_MINAREA %s -CHECKIMAGE_NAME checki.fits,segmenti.fits'%\
    (field,'sex.param',sg_dir+'/%s_catalog.fits'%(field),str(seeing),str(param['satur_level_i']),str(1.1/1.7*np.pi*(seeing/0.22)**2))
    sub = subprocess.check_call(shlex.split(cmd)); print cmd

    cmd='sex final/coadd_c%s_i.fits,final/coadd_c%s_sq_i.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -CATALOG_TYPE FITS_1.0 -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE 0.22 -DETECT_MINAREA %s'%\
    (field,field,'sex.param',sg_dir+'/%s_sq_catalog.fits'%(field),str(seeing),str(param['satur_level_i']),str(1.1/1.7*np.pi*(seeing/0.22)**2))
    sub = subprocess.check_call(shlex.split(cmd)); print cmd
