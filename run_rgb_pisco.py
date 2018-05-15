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
    g = (g - gmean) / gstd#* 2
    gmin=(gmean-gstd)/gstd

    r = fits.open('test.fits')[0].data[1]
    r_nonan = r.ravel()[~np.isnan(r.ravel())]
    rmean = np.mean(r_nonan[np.abs(r_nonan) < 150])
    rstd = np.std(r_nonan[np.abs(r_nonan) < 150])
    r = (r - rmean) / rstd#* 2
    rmin=(rmean-rstd)/rstd

    i = fits.open('test.fits')[0].data[0]
    i_nonan = i.ravel()[~np.isnan(i.ravel())]
    imean = np.mean(i_nonan[np.abs(i_nonan) < 150])
    istd = np.std(i_nonan[np.abs(i_nonan) < 150])
    i = (i - imean) / istd * 1.3
    imin=(imean-istd)/istd

    cmd = "rm test.fits"; print cmd; sub = subprocess.check_call(shlex.split(cmd))

    # rgb_default = make_lupton_rgb(i, r, g, Q=7, stretch=4)
    rgb_default = make_lupton_rgb(i, r, g, Q=6, stretch=7, minimum=[imin,rmin,gmin])

    fig, ax = plt.subplots(1, figsize=(20, 30))
    ax.imshow(rgb_default, origin='lower')
    # ax.plot([200, 328.57], [820, 820], color='white', lw=3)
    # ax.annotate('30"', xy=(200, 800 - 28), xycoords='data',
    #             color='white', fontsize=15)
    #
    # if redshift == -1:
    #     ax.annotate('no z', xy=(200, 800 + 33),
    #                 xycoords='data', color='white', fontsize=14)
    # else:
    #     kpc = 30 * (1 / cosmo.arcsec_per_kpc_proper(redshift)).value
    #     ax.annotate('z=%.2f: %.0f kpc' % (redshift, kpc), xy=(
    #         200, 800 + 33), xycoords='data', color='white', fontsize=14);

    # ax.axis([100,1570,730,1700])
    ax.axis([0,1600,400,2100])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axvline(800, c='lightgreen',ls='--')
    plt.axhline(1250, c='lightgreen',ls='--')
    plt.tight_layout()
    plt.savefig('/Users/taweewat/Documents/pisco_code/Chips_images/large_%s_img_.jpeg' %
                name, bbox_inches='tight')#,dpi='figure')
    return 0


# name='CHIPS0022-5140'

# Problems
# CHIPS1011-0505: the i band does not have a correct wcs, the g band is extremely faint
# CHIPS2357-1125


# fixed 35 CHIPS2357-1125: still don't know the problem
# fixed 36 CHIPS0140-1533: still don't know the problem  #work!!
# fixed 59 CHIPS0152-5028: very bad alignment, run wcs again
# fixed 118 CHIPS0106-2358: stil bad alignment

# 43 CHIPS1011-0505: still don't know the problem (can't do the alignment properly), not working, just can't align it
# 102 CHIPS0132-1608: bad alignment, the star is way too bright to get any good measurement

# New: CHIPS0107-1310, CHIPS0050-5249, CHIPS0112-2919, CHIPS0116-1136, CHIPS0127-4016
# Not: CHIPS0132-1608,
if __name__ == "__main__":

    mode = 'chips'

    if mode == 'chips':
        allnames = list_file_name(
            '/Users/taweewat/Documents/pisco_code/final/', 'coadd_cCHIPS', end='i.fits')
        base = pd.read_csv(
            '/Users/taweewat/Documents/red_sequence/chips_all_obj.csv', index_col=0)
    elif mode == 'field':
        allnames = list_file_name(
            '/Users/taweewat/Documents/pisco_code/final/', 'coadd_cField', end='i.fits')
        base = pd.read_csv('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs.csv')

    all_names = list(set([i.split('/')[-1][7:-7] for i in allnames]))
    print 'the total number of objects:', len(all_names)

    # chips='CHIPS0229-5232'
    # print chips, all_names[0]
    # print np.where(np.array(all_names) == chips)[0][0]; ind=np.where(np.array(all_names) == chips)[0][0]

    # field='Field060'
    # print field, all_names[0]
    # print np.where(np.array(all_names) == field)[0][0]; ind=np.where(np.array(all_names) == field)[0][0]

    #DONE: 'Field018','Field020','Field021','Field022','Field025','Field026','Field027']
    # all_names=['CHIPS0609-0247']
               # 'Field028','Field029','Field030','Field033','Field034','Field036','Field037','Field038','Field039','Field040','Field042','Field044','Field045']

    for i, name in enumerate(np.append(all_names[:43],all_names[44:])):
    # for i, name in enumerate(all_names[ind:ind+1]):

        if mode == 'chips':
            redshift = base[base.chips == name].redshift.values[0]
            RA = base[base.chips == name].ra.values[0]
            DEC = base[base.chips == name].dec.values[0]
        elif mode == 'field':
            redshift = base[base.name == name].redshift.values[0]
            RA = base[base.name == name].ra.values[0]
            DEC = base[base.name == name].dec.values[0]

        print i, name
        ## make image with ds9 with the center is crosshaired
        # cmd = "ds9 -zscale -crosshair %f %f wcs fk5 -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom out -saveimage Chips_images/%s_ds9.eps -exit" % \
        #     (RA, DEC, name, name, name, name)
        # print cmd
        # sub = subprocess.check_call(shlex.split(cmd))

        ## make jpeg image from the python script with the scale size
        if not os.path.isfile(os.path.join('Chips_images', 'large_%s_img_.jpeg'%name)):
            print i, 'working on the the image large_%s_img_.jpeg'%name
            a=make_images_jpeg(name, redshift)
        # else:
            # make_images_jpeg(name, redshift)
