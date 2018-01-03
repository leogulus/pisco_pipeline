import numpy as np
import matplotlib
import os
import pandas as pd

from astropy.io import fits
import subprocess
import cosmics
import shlex
import sys

import matplotlib.pyplot as plt
#from pisco_lib import *
# edited 8/30/17

##----
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
    # cut = -27

    ch1_name = list_file_name(dir, fieldname)
    print 'working on %s with the index=%i' % (ch1_name[0], index)
    hdulist = fits.open(ch1_name[0])
    ch1 = hdulist[index].data

    bias_names = list_file_name(dir, 'Bias_')
    if flat == 'domeflat':
        domeflat_names = list_file_name(dir, "domeflat" + filter_name(index)[1])
    if flat == 'twilight':
        domeflat_names = list_file_name(dir, "twiflat_")

    bias = open_files(bias_names, index)
    if flat == 'domeflat':
        domeflat = open_files(domeflat_names, index, bias=bias, twilight=False)
    if flat == 'twilight':
        domeflat = open_files(domeflat_names, index, bias=bias, twilight=True)

    domeflat[domeflat == 0] = 1e-4

    img = (ch1 - bias) / domeflat
    ch1, bias, domeflat, img = ch1[:, :], bias[:, :], domeflat[:, :], img[:, :]

    # if index % 2 == 0:
    #     return np.fliplr(ch1), np.fliplr(bias), np.fliplr(domeflat), np.fliplr(img)
    # else:
    return ch1, bias, domeflat, img

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

    if len(sys.argv)>3:
        flattype = str(sys.argv[3])
    else:
        flattype='domeflat'

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

    Combine two amplifiers and bias and flat fielding
    for field in fields:
        for index in [1, 3, 5, 7]:
            ch1, bias1, domeflat1, img1 = reduce_data(dir, index, field, flat=flattype)
            ch2, bias2, domeflat2, img2 = reduce_data(dir, index + 1, field, flat=flattype)
            final_image = np.concatenate((img1, img2), axis=1)
            save_fits(index, dir, reducedir, field, final_image,
                      "%s_%s.fits" % (field, filter_name(index)[0]))
