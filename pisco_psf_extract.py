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

import FITS_tools
import yaml

import matplotlib.pyplot as plt

"""
pisco_psf_extract: including 3 steps
1. sextractor to get appropriate catalog for psfex input file
2. run psfex to fit the psf profile to get the best .psf file
3. reapply .psf file to sextractor to get MAG_PSF, MAG_MODEL, etc

Examples: python pisco_pipeline/pisco_psf_extract.py CHIPS0005-2758 i
"""

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


def find_seeing_fits(field):
    home='/Users/taweewat/Documents/pisco_code/'
    dirs=['ut170103/','ut170104/','ut170619/','ut170621/','ut170624/','ut171208/',\
    'ut171209/','ut171212/','ut190412/','ut190413/']
    myReg=re.compile(r'(%s_A).*'%field)
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                seeing=float(fits.open(dir+myReg.search(text).group())[0].header['FWHM1'])
    return seeing

##------


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

field = str(sys.argv[1])

band = str(sys.argv[2])

with open("pisco_pipeline/params.yaml", 'r') as stream:
    try:
        param=yaml.load(stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print(exc)

# seeing = float(fits.open(list_file_name_seeing(
#     '/Users/taweewat/Documents/pisco_code/', field, startdir='ut')[0])[0].header['FWHM1'])
seeing=find_seeing_fits(field)

# df_see=pd.read_csv('/Users/taweewat/Documents/red_sequence/total_chips_field_seeing.csv',index_col=0)
# if field[0:5]=='CHIPS':
#     seeing = df_see[df_see.chips==field]['seeing_q25'].values[0] #np.min(df_see[df_see.chips==field][['seeing_025','seeing_gra_025']].values)
#     print seeing
# elif (field[0:5]=='Field')|(field[0:3]=='PKS'):
#     seeing = df_see[df_see.name==field]['seeing_q25'].values[0] #np.min(df_see[df_see.name==field][['seeing_025','seeing_gra_025']].values)
#     print seeing
#

# slrdir = 'slr_output'
# to_be_projected = 'final/coadd_c%s_%s.fits'%(field,band)
# reference_fits  = 'final/coadd_c%s_i.fits'%field
# im1,im2, header = FITS_tools.match_fits(to_be_projected,reference_fits,return_header=True)
# outname = 'final/proj_coadd_c%s_%s.fits'%(field,band)
# print 'projecting from %s band to i band the fits file '%band + outname
# fits.writeto(outname, im1, header, overwrite=True)

pxscale=0.11
minarea=1.7

# cmd='sex final/coadd_c%s_%s.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PIXEL_SCALE %s -PHOT_APERTURES 15'%\
# (field,band,'sex_psf.param','psf_%s_%s.fits'%(field,band),str(seeing),str(param['satur_level_%s'%band]),str(pxscale))


cmd='sex final/coadd_c%s_%s.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 23 -PIXEL_SCALE %s -DETECT_MINAREA %s'%\
(field,band,'sex_psf.param','psfex_output/psf_%s_%s.fits'%(field,band),str(seeing),str(param['satur_level_%s_psf'%band]),str(pxscale),str(1.1/minarea*np.pi*(seeing/pxscale)**2)) 


print cmd
sub = subprocess.check_call(shlex.split(cmd))

# cmd='psfex %s -c pisco_pipeline/pisco.psfex' % ('psf_%s_%s.fits'%(field,band))
# print cmd
# sub = subprocess.check_call(shlex.split(cmd))
#
# cmd='sex final/coadd_c%s_i.fits -c pisco_pipeline/config.sex -PSF_NAME %s -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -GAIN 0.5 -PIXEL_SCALE 0.2'%\
# (field,'psf_%s_%s.psf'%(field,band),'sex_after_psf.param','a_psf_%s_%s.fits'%(field,band),str(seeing),str(param['satur_level_%s'%band]))
# print cmd
# sub = subprocess.check_call(shlex.split(cmd))
