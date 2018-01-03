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
