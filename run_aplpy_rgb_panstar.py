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
import re

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))

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


def make_images_aplpy(name, RA, DEC, redshift, mode='chips',  RA_W=0., RA_M=0., DEC_W=0., DEC_M=0.):

    data, header = fits.getdata('/Users/taweewat/Documents/xray_project/panstar/panstar_%s_i.fits'%name, header=True)
    fits.writeto('/Users/taweewat/Documents/xray_project/panstar/panstar_%s_i0.fits'%name, data, header=header, overwrite=True)
    data, header = fits.getdata('/Users/taweewat/Documents/xray_project/panstar/panstar_%s_r.fits'%name, header=True)
    fits.writeto('/Users/taweewat/Documents/xray_project/panstar/panstar_%s_r0.fits'%name, data, header=header, overwrite=True)
    data, header = fits.getdata('/Users/taweewat/Documents/xray_project/panstar/panstar_%s_g.fits'%name, header=True)
    fits.writeto('/Users/taweewat/Documents/xray_project/panstar/panstar_%s_g0.fits'%name, data, header=header, overwrite=True)

    aplpy.make_rgb_cube(['/Users/taweewat/Documents/xray_project/panstar/panstar_%s_i0.fits'%name,\
                     '/Users/taweewat/Documents/xray_project/panstar/panstar_%s_r0.fits'%name,\
                     '/Users/taweewat/Documents/xray_project/panstar/panstar_%s_g0.fits'%name],\
                    'test.fits', north=True)

    g=fits.open('test.fits')[0].data[2]
    g_nonan=g.ravel()[~np.isnan(g.ravel())]
    gmean=np.mean(g_nonan[np.abs(g_nonan)<8])
    gstd=np.std(g_nonan[np.abs(g_nonan)<8])
    if mode=='chips':
        gmin=(gmean-0.5*gstd); print gmin
    elif mode=='field':
        gmin=(gmean-gstd); print gmin
    gmax=(gmean+10*gstd); print gmax  #20

    r=fits.open('test.fits')[0].data[1]
    r_nonan=r.ravel()[~np.isnan(r.ravel())]
    rmean=np.mean(r_nonan[np.abs(r_nonan)<8])
    rstd=np.std(r_nonan[np.abs(r_nonan)<8])
    if mode=='chips':
        rmin=(rmean-0.5*rstd); print rmin
    elif mode=='field':
        rmin=(rmean-rstd); print rmin
    rmax=(rmean+10*rstd); print rmax #20

    i=fits.open('test.fits')[0].data[0]
    i_nonan=i.ravel()[~np.isnan(i.ravel())]
    imean=np.mean(i_nonan[np.abs(i_nonan)<8])
    istd=np.std(i_nonan[np.abs(i_nonan)<8])
    if mode=='chips':
        imin=(imean-0.5*istd); print imin
    elif mode=='field':
        imin=(imean-istd); print imin
    imax=(imean+10*istd); print imax #20, 12

    aplpy.make_rgb_image('test.fits','rgb_image_arcsinh.jpeg',\
                         stretch_r='log', stretch_g='log',stretch_b='log',\
                         vmin_r=imin, vmin_g=rmin, vmin_b=gmin,
                         vmax_r=imax, vmax_g=rmax, vmax_b=gmax)
                         # pmin_r=8, pmin_g=8, pmin_b=8,\
                         # pmax_r=99.6, pmax_g=99.75, pmax_b=99.75)

    #fix the dimension of the fits file so that it can be read by make_rgb_image
    data, header = fits.getdata("test.fits", header=True)
    header.pop('NAXIS3')
    header.set('NAXIS',2)
    data=data[0]
    fits.writeto('test_3.fits', data, header=header, overwrite=True)

    # img = aplpy.FITSFigure('/Users/taweewat/Documents/pisco_code/final/coadd_c%s_i.fits'%name)
    img = aplpy.FITSFigure('test_3.fits')
    img.show_rgb('rgb_image_arcsinh.jpeg')
    img.tick_labels.set_xformat('ddmm')
    img.tick_labels.set_yformat('ddmm')
    # img.add_grid()
    # img.grid.set_color('green')
    img.add_scalebar(4/60.)
    img.scalebar.set_label('4 arcmin')
    img.scalebar.set_color('white')
    img.recenter(RA, DEC, width=0.067, height=0.067)
    # img.show_markers(RA, DEC, marker='x', s=200, lw=0.5, layer='markers', edgecolor='white', facecolor='white')
    img.show_markers(RA_W, DEC_W, marker='o', s=200, lw=0.5, layer='markers3', edgecolor='white', facecolor='none')

    img.save('Chips_images/aplpy_panstar_%s_img.jpeg'%name, adjust_bbox=True, dpi=120)

    # purge('/Users/taweewat/Documents/xray_project/panstar/', 'panstar_%s_i0.fits'%name)
    # purge('/Users/taweewat/Documents/xray_project/panstar/', 'panstar_%s_r0.fits'%name)
    # purge('/Users/taweewat/Documents/xray_project/panstar/', 'panstar_%s_g0.fits'%name)


    # purge('/Users/taweewat/Documents/xray_project/panstar/', "panstar_%s_.*\0.fits" % name)

    # g = fits.open('test.fits')[0].data[2]
    # g_nonan = g.ravel()[~np.isnan(g.ravel())]
    # gmean = np.mean(g_nonan[np.abs(g_nonan) < 150])
    # gstd = np.std(g_nonan[np.abs(g_nonan) < 150])
    # g = (g - gmean) / gstd#* 2
    # gmin=(gmean-gstd)/gstd
    #
    # r = fits.open('test.fits')[0].data[1]
    # r_nonan = r.ravel()[~np.isnan(r.ravel())]
    # rmean = np.mean(r_nonan[np.abs(r_nonan) < 150])
    # rstd = np.std(r_nonan[np.abs(r_nonan) < 150])
    # r = (r - rmean) / rstd#* 2
    # rmin=(rmean-rstd)/rstd
    #
    # i = fits.open('test.fits')[0].data[0]
    # i_nonan = i.ravel()[~np.isnan(i.ravel())]
    # imean = np.mean(i_nonan[np.abs(i_nonan) < 150])
    # istd = np.std(i_nonan[np.abs(i_nonan) < 150])
    # i = (i - imean) / istd * 1.3
    # imin=(imean-istd)/istd
    #
    # cmd = "rm test.fits"; print cmd; sub = subprocess.check_call(shlex.split(cmd))
    #
    # # rgb_default = make_lupton_rgb(i, r, g, Q=7, stretch=4)
    # rgb_default = make_lupton_rgb(i, r, g, Q=6, stretch=7, minimum=[imin,rmin,gmin])
    #
    # fig, ax = plt.subplots(1, figsize=(20, 30))
    # ax.imshow(rgb_default, origin='lower')
    # ax.axis([0,1600,400,2100])
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    # plt.axvline(800, c='lightgreen',ls='--')
    # plt.axhline(1250, c='lightgreen',ls='--')
    # plt.tight_layout()
    # plt.savefig('/Users/taweewat/Documents/pisco_code/Chips_images/large_%s_img_.jpeg' %
    #             name, bbox_inches='tight')#,dpi='figure')
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


    total=pd.read_csv('/Users/taweewat/Documents/red_sequence/chips_all_obj.csv',index_col=0)
    total=total[total.panstar==True]
    all_names = list(total.chips.values)
    # all_names=['CHIPS0004-4736','CHIPS0005-2758','CHIPS0024-6820','CHIPS0302-2758','CHIPS0512-3257','CHIPS0957-7554','CHIPS0137-1248','CHIPS2210-5508','CHIPS2217-3034','CHIPS2303-6807','CHIPS2325-4800']
    # all_names=['Field101']
    # all_names=['Field025','Field053','Field058','Field062','Field063','Field066','Field086','Field094','Field101','Field105','Field159','Field208','Field217','Field292']

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

    # for i, name in enumerate(np.append(all_names[:43],all_names[44:])):
    all_names=['CHIPS1011-0505','CHIPS1911+4455']

    for i, name in enumerate(all_names[:]):
    # for i, name in enumerate(all_names[ind:ind+1]):

        # print i, name
        # if name=='CHIPS2249-2808':
        #     name='CHIPS2227-4333'
        # elif name=='CHIPS2246-2854':
        #     name='CHIPS2223-3455'

        if mode == 'chips':
            redshift = base[base.chips == name].redshift.values[0]
            RA = base[base.chips == name].ra.values[0]
            DEC = base[base.chips == name].dec.values[0]
            RA_W = base[base.chips == name].ra_w.values[0]
            DEC_W = base[base.chips == name].dec_w.values[0]

        elif mode == 'field':
            redshift = base[base.name == name].redshift.values[0]
            RA = base[base.name == name].ra.values[0]
            DEC = base[base.name == name].dec.values[0]
            RA0 = base[base.name == name].RA0.values[0]
            DEC0 = base[base.name == name].DEC0.values[0]
            RA_W = base[base.name == name].ra_w.values[0]
            DEC_W = base[base.name == name].dec_w.values[0]
            RA_M = base[base.name == name].ra_m.values[0]
            DEC_M = base[base.name == name].dec_m.values[0]

        # if name=='CHIPS2227-4333':
        #     name='CHIPS2249-2808'
        # elif name=='CHIPS2223-3455':
        #     name='CHIPS2246-2854'

        print i, name
        ## make image with ds9 with the center is crosshaired
        # cmd = "ds9 -zscale -crosshair %f %f wcs fk5 -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom out -saveimage Chips_images/%s_ds9.eps -exit" % \
        #     (RA, DEC, name, name, name, name)
        # print cmd
        # sub = subprocess.check_call(shlex.split(cmd))

        ## make jpeg image from the python script with the scale size
        if not os.path.isfile(os.path.join('Chips_images', 'aplpy2_panstar_%s_img.jpeg'%name)):
            print i, 'working on the the image aplpy2_panstar_%s_img.jpeg'%name
            if mode=='chips':
                a=make_images_aplpy(name, RA, DEC, redshift, mode, RA_W, 0, DEC_W, 0)
            elif mode=='field':
                a=make_images_aplpy(name, RA, DEC, redshift, mode, RA_W, RA_M, DEC_W, DEC_M)
        # else:
            # make_images_jpeg(name, redshift)
