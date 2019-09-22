import os, re
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

"""
python pisco_pipeline/run_aplpy_rgb_pisco.py
"""

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
    aplpy.make_rgb_cube(['/Users/taweewat/Documents/pisco_code/final/coadd_c%s_i.fits' % name,\
    '/Users/taweewat/Documents/pisco_code/final/coadd_c%s_r.fits' % name,
    '/Users/taweewat/Documents/pisco_code/final/coadd_c%s_g.fits' % name], 'test.fits', north=True)

    g=fits.open('test.fits')[0].data[2]
    g_nonan=g.ravel()[~np.isnan(g.ravel())]
    gmean=np.mean(g_nonan[np.abs(g_nonan)<150])
    gstd=np.std(g_nonan[np.abs(g_nonan)<150])
    # g=(g-gmean)/gstd
    if mode=='chips':
        gmin=(gmean-gstd); print gmin
    elif mode=='field':
        gmin=(gmean-gstd); print gmin
    elif (mode == 'pks') | (mode == 'sdss'):
        gmin=(gmean-gstd); print gmin
    # gmax=(gmean+10*gstd); print gmax
    gmax=(gmean+20*gstd); print gmax  #20

    r=fits.open('test.fits')[0].data[1]
    r_nonan=r.ravel()[~np.isnan(r.ravel())]
    rmean=np.mean(r_nonan[np.abs(r_nonan)<150])
    rstd=np.std(r_nonan[np.abs(r_nonan)<150])
    # r=(r-rmean)/rstd
    if mode=='chips':
        rmin=(rmean-rstd); print rmin
    elif mode=='field':
        rmin=(rmean-rstd); print rmin
    elif (mode=='pks') | (mode=='sdss'):
        rmin=(rmean-rstd); print rmin
    # rmax=(rmean+10*rstd); print rmax 
    rmax=(rmean+20*rstd); print rmax #20

    i=fits.open('test.fits')[0].data[0]
    i_nonan=i.ravel()[~np.isnan(i.ravel())]
    imean=np.mean(i_nonan[np.abs(i_nonan)<150])
    istd=np.std(i_nonan[np.abs(i_nonan)<150])
    # i=(i-imean)/istd*1.3
    if mode=='chips':
        imin=(imean-istd); print imin
    elif mode=='field':
        imin=(imean-istd); print imin
    elif (mode == 'pks') | (mode == 'sdss'):
        imin=(imean-istd); print imin
    # imax=(imean+6*istd); print imax
    imax=(imean+20*istd); print imax  # 12, 8
    

    aplpy.make_rgb_image('test.fits','rgb_image_arcsinh.png',\
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
    img.show_rgb('rgb_image_arcsinh.png')
    img.tick_labels.set_xformat('ddmm')
    img.tick_labels.set_yformat('ddmm')
    # img.add_grid()
    # img.grid.set_color('green')
    # img.add_scalebar(1/60.)
    img.add_scalebar(86.48 / 3600.)
    # img.scalebar.set_length(86.48/3600.)
    img.scalebar.set_label('300 kpc')
    img.scalebar.set_color('white')
    img.recenter(RA, DEC, width=0.067, height=0.067)
    # img.show_markers(RA, DEC, marker='x', s=150, lw=0.4, layer='markers', edgecolor='white', facecolor='white')
    # if mode=='field':
        # img.show_markers(RA_M, DEC_M, marker='P', s=90, lw=2, layer='markers2', edgecolor='red', facecolor='none')
    # img.show_markers(RA_W, DEC_W, marker='o', s=150, lw=0.4, layer='markers3', edgecolor='white', facecolor='none')
    img.set_theme('publication')
    img.axis_labels.hide()
    img.ticks.hide()
    img.tick_labels.hide()
    # filename = 'Chips_images/aplpy6_%s_img.jpeg' % name
    img.save(filename, adjust_bbox=True, dpi=120)

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


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))


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

    # mode = 'chips'
    # mode = 'field'

    # if mode == 'chips':
    #     allnames = list_file_name(
    #         '/Users/taweewat/Documents/pisco_code/final/', 'coadd_cCHIPS', end='i.fits')
    #     base = pd.read_csv(
    #         '/Users/taweewat/Documents/red_sequence/chips_all_obj.csv', index_col=0)
    # elif mode == 'field':
    #     allnames = list_file_name(
    #         '/Users/taweewat/Documents/pisco_code/final/', 'coadd_cField', end='i.fits')
    #     base = pd.read_csv('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs.csv')

    # all_names = list(set([i.split('/')[-1][7:-7] for i in allnames]))

    # chips='CHIPS0229-5232'
    # print chips, all_names[0]
    # print np.where(np.array(all_names) == chips)[0][0]; ind=np.where(np.array(all_names) == chips)[0][0]

    # field='Field060'
    # print field, all_names[0]
    # print np.where(np.array(all_names) == field)[0][0]; ind=np.where(np.array(all_names) == field)[0][0]

    #DONE: 'Field018','Field020','Field021','Field022','Field025','Field026','Field027']
    # all_names=['Field094','Field159'] #'CHIPS1011-0505','CHIPS0005-2758'
               # 'Field028','Field029','Field030','Field033','Field034','Field036','Field037','Field038','Field039',\
               # 'Field040','Field042','Field044','Field045']
    # all_names=['CHIPS2210-5508','CHIPS2217-3034']
    # all_names = ['Field227']
    # all_names = ['CHIPS0132-1608']

    dirs = ['ut170103/','ut170104/','ut170619/','ut170621/','ut170624/','ut171208/','ut171209/','ut171212/']
    dirs = ['ut190412/','ut190413']
    home='/Users/taweewat/Documents/pisco_code/'
    names=[]
    myReg=re.compile(r'(CHIPS\d{4}[+-]\d{4})|(Field\d{3})')
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                names.append(myReg.search(text).group())
    all_fields=list(set(names))
    # all_fields = ['CHIPS1950-4000']

    all_fields=['SDSS123', 'SDSS501', 'SDSS603']
    all_fields_cut = all_fields[:]

    for i, name in enumerate(all_fields_cut):
        print str(i) + '/' + str(len(all_fields_cut)), name

        if name in ['CHIPS1933-1511']:
            continue

        if name=='CHIPS2249-2808':
            name='CHIPS2227-4333'
        elif name=='CHIPS2246-2854':
            name='CHIPS2223-3455'

        if name[0:5]=='CHIPS':
            mode='chips'
            base = pd.read_csv('/Users/taweewat/Documents/red_sequence/chips_all_obj.csv', index_col=0)
        elif name[0:5]=='Field':
            mode='field'
            base = pd.read_csv('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs.csv')
        elif name[0:3]=='PKS':
            mode='pks'
        elif name[0:4] == 'SDSS':
            mode = 'sdss'

        #for new targets ut190412
        base=pd.read_csv('/Users/taweewat/Documents/xray_project/ned-result/total_776_new_pan.csv',index_col=0)

        if name=='CHIPS2227-4333':
            name='CHIPS2249-2808'
        elif name=='CHIPS2223-3455':
            name='CHIPS2246-2854'


        # if name=='CHIPS1423-3125':
        #     RA,DEC=215.88332999999997,-31.423749999999998
        # elif name=='CHIPS2007-4434':
        #     RA,DEC=301.98419,-44.58055
        # elif name=='CHIPS1605-3115':
        #     RA,DEC=241.4466667,-31.2577778
        # else: 
        #     RA = base[base.chips == name].ra.values[0]
        #     DEC = base[base.chips == name].dec.values[0] 

        # if mode == 'chips':
        #     redshift = base[base.chips == name].redshift.values[0]
        #     RA = base[base.chips == name].ra.values[0]
        #     DEC = base[base.chips == name].dec.values[0]
        #     RA_W = base[base.chips == name].ra_w.values[0]
        #     DEC_W = base[base.chips == name].dec_w.values[0]
        # elif mode == 'field':
        #     redshift = base[base.name == name].redshift.values[0]
        #     RA = base[base.name == name].ra.values[0]
        #     DEC = base[base.name == name].dec.values[0]
        #     RA0 = base[base.name == name].RA0.values[0]
        #     DEC0 = base[base.name == name].DEC0.values[0]
        #     RA_W = base[base.name == name].ra_w.values[0]
        #     DEC_W = base[base.name == name].dec_w.values[0]
        #     RA_M = base[base.name == name].ra_m.values[0]
        #     DEC_M = base[base.name == name].dec_m.values[0]
        # elif mode == 'pks':
        #     redshift=0.2230
        #     RA = 209.0225
        #     DEC = -34.3530556
        
        if mode == 'sdss':
            de = pd.read_csv(
                '/Users/taweewat/Documents/xray_project/ned-result/final_sdss_cut5.csv', index_col=0)
            RA = de[de.name == name].RA.values[0]
            DEC = de[de.name == name].DEC.values[0]
            redshift = de[de.name == name].redshift.values[0]


        print i, name

        ## make jpeg image from the python script with the scale size
        global filename
        filename = 'Chips_images/aplpy4_%s_img4.jpeg' % name
        purge('.',filename)

        if not os.path.isfile(filename):
            print i, 'working on the the image %s'%filename
            if mode=='chips':
                # a=make_images_aplpy(name, RA, DEC, redshift, mode, RA_W, 0, DEC_W, 0)
                a=make_images_aplpy(name, RA, DEC, -1, mode, 0, 0, 0, 0)

            elif mode=='field':
                a=make_images_aplpy(name, RA, DEC, redshift, mode, RA_W, RA_M, DEC_W, DEC_M)
            elif mode=='pks':
                a=make_images_aplpy(name, RA, DEC, redshift, mode, 0, 0, 0, 0)
            elif mode=='sdss':
                a = make_images_aplpy(name, RA, DEC, redshift, mode, 0, 0, 0, 0)
        # else:
            # make_images_jpeg(name, redshift)
