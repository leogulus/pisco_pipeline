import os
import subprocess
import shlex
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
from astropy.visualization import make_lupton_rgb

import aplpy


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


def make_images_jpeg(name, redshift):
    aplpy.make_rgb_cube(['/Users/taweewat/Documents/pisco_code/final/coadd_c%s_i.fits' % name, '/Users/taweewat/Documents/pisco_code/final/coadd_c%s_r.fits' % name,
                         '/Users/taweewat/Documents/pisco_code/final/coadd_c%s_g.fits' % name], 'test.fits', north=True)

    g = fits.open('test.fits')[0].data[2]
    g_nonan = g.ravel()[~np.isnan(g.ravel())]
    gmean = np.mean(g_nonan[np.abs(g_nonan) < 150])
    gstd = np.std(g_nonan[np.abs(g_nonan) < 150])
    g = (g - gmean) / gstd

    r = fits.open('test.fits')[0].data[1]
    r_nonan = r.ravel()[~np.isnan(r.ravel())]
    rmean = np.mean(r_nonan[np.abs(r_nonan) < 150])
    rstd = np.std(r_nonan[np.abs(r_nonan) < 150])
    r = (r - rmean) / rstd

    i = fits.open('test.fits')[0].data[0]
    i_nonan = i.ravel()[~np.isnan(i.ravel())]
    imean = np.mean(i_nonan[np.abs(i_nonan) < 150])
    istd = np.std(i_nonan[np.abs(i_nonan) < 150])
    i = (i - imean) / istd

    cmd = "rm test.fits"
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    rgb_default = make_lupton_rgb(i, r, g, Q=7, stretch=4)

    fig, ax = plt.subplots(1, figsize=(10, 20))
    ax.imshow(rgb_default, origin='lower')
    ax.plot([200, 328.57], [800, 800], color='white', lw=3)
    ax.annotate('30"', xy=(200, 800 - 58), xycoords='data',
                color='white', fontsize=15)

    if redshift == -1:
        ax.annotate('no z', xy=(200, 800 + 23),
                    xycoords='data', color='white', fontsize=14)
    else:
        kpc = 30 * (1 / cosmo.arcsec_per_kpc_proper(redshift)).value
        ax.annotate('z=%.2f: %.0f kpc' % (redshift, kpc), xy=(
            200, 800 + 23), xycoords='data', color='white', fontsize=14);

    ax.axis([0, 1600, 700, 1800])
    plt.tight_layout()
    plt.savefig('/Users/taweewat/Documents/pisco_code/Chips_images/%s_img.jpeg' %
                name, bbox_inches='tight')
    return 0


# name='CHIPS0022-5140'

# Problems
# CHIPS1011-0505: the i band does not have a correct wcs
# CHIPS0609-0247: wrong alignment for one of band (i)
# CHIPS1036-3513: error: can't allocate region
# CHIPS0552-2336: one of the band is misalingned
# CHIPS2133-2815: the i band does not have a correct wcs
# CHIPS0724-0715: Super dense field, not enough to get a good wcs

if __name__ == "__main__":

    mode = 'field'

    if mode == 'chips':
        allnames = list_file_name(
            '/Users/taweewat/Documents/pisco_code/final/', 'coadd_cCHIPS', end='i.fits')
        base = pd.read_csv(
            '/Users/taweewat/Documents/xray_project/red_sequence/chips_all_obj.csv', index_col=0)
    elif mode == 'field':
        allnames = list_file_name(
            '/Users/taweewat/Documents/pisco_code/final/', 'coadd_cField', end='i.fits')
        all_names = list(set([i.split('/')[-1][7:-7] for i in allnames]))
        base = pd.read_csv('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs.csv')

    print 'the total number of objects:', len(all_names)
    for i, name in enumerate(all_names[:]):

        if mode == 'chips':
            redshift = base[base.chips == name].redshift.values[0]
            RA = base[base.chips == name].ra.values[0]
            DEC = base[base.chips == name].dec.values[0]
        elif mode == 'field':
            redshift = base[base.name == name].redshift.values[0]
            RA = base[base.name == name].ra.values[0]
            DEC = base[base.name == name].dec.values[0]

        print i, name
        # make image with ds9 with the center is crosshaired
        cmd = "ds9 -zscale -crosshair %f %f wcs fk5 -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom out -saveimage Chips_images/%s_ds9.eps -exit" % \
            (RA, DEC, name, name, name, name)
        print cmd
        sub = subprocess.check_call(shlex.split(cmd))

        # make jpeg image from the python script with the scale size
        if not os.path.isfile(os.path.join('Chips_images', '%s_img.jpeg'%name)):
            print 'working on the the image %s_img.jpeg'%name
            make_images_jpeg(name, redshift)
