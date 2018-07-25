import sys, os, re, yaml, subprocess, shlex, FITS_tools
import pandas as pd
import numpy as np

from astropy.table import Table, join


def pisco_photometry_psf_v4_nostar(field, mode='psf', slr=True):
    def update_color(table, mode='psf'):
        """
        update_color: using the output from SLR, update to the correct magnitude
        INPUT:
        - fname: input file from SLR output (...offsets.list)
        - table: the table that we want to update the value (from column magg,etc to MAG_g,etc)
        OUTPUT:
        - a new table with added columns with name MAG_g,...,MAGERR_g,...
        """
        band=['g','r','i','z']
        print 'bands = ', band

        if mode=='psf':
            MODE1='PSF'
        elif mode=='model':
            MODE1='MODEL'
        elif mode=='auto':
            MODE1='AUTO'
        elif mode=='aper':
            MODE1='APER'
        elif mode=='hybrid':
            MODE1='HYBRID'

        for i in range(4):
            table['MAG_' + band[i]] = table['MAG_%s_'%MODE1 + band[i]]# + corr[0]
            table['MAGERR_' + band[i]] = (table['MAGERR_%s_'%MODE1 + band[i]])#**2 + ecorr[0]**2)**0.5
        return table
    slrdir = 'slr_output'
    total_gal=Table.from_pandas(pd.read_csv("/Users/taweewat/Documents/pisco_code/slr_output/galaxy_psf_total_%s.csv"%field))
    print 'mode=', mode
    ntotal_gal = update_color(total_gal, mode=mode)
    ntotal_gal.write(os.path.join(slrdir, 'galaxy_psf_ntotal_nostar_%s.csv' % field), overwrite=True)

##-----------
if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)
    field = str(sys.argv[1])
    mode = str(sys.argv[2])

    all_fields=['Field166','Field054','Field042','Field037','Field060','Field105','Field116','CHIPS0411-5149',\
                'CHIPS0040-2902','CHIPS0025-5427','CHIPS0005-2758','CHIPS0150-3338','CHIPS1034-2837','CHIPS0409-2839',\
                'CHIPS0022-5140','CHIPS0512-3257','CHIPS0536-3401','CHIPS0514-5046','CHIPS0050-5249','CHIPS2227-4333',\
                'CHIPS2228-3220','CHIPS2306-3439','CHIPS2307-4236','CHIPS0004-4736','CHIPS2303-6807','CHIPS2348-3831',\
                'CHIPS0449-4350','CHIPS0449-3910','CHIPS2254-3635','CHIPS2243-3034','CHIPS2251-3210','CHIPS2311-4718',\
                'CHIPS0335-4715','CHIPS0355-6645','CHIPS0024-6820','CHIPS0342-3703','CHIPS2325-4800','CHIPS1011-0505',\
                'CHIPS0222-4159','CHIPS0302-2758','CHIPS0304-3556','CHIPS0157-1043','CHIPS0253-5441','CHIPS0552-2336',\
                'CHIPS0219-3626','CHIPS0229-5232','CHIPS1141-1407','CHIPS1142-1422','CHIPS0325-4926','CHIPS2218-2947',\
                'CHIPS0152-5028','CHIPS0153-3143','CHIPS0316-2247','CHIPS0724-0715','CHIPS1205-2633','CHIPS0827-2026',\
                'CHIPS0847-0703','CHIPS0012-1628','CHIPS0140-1533','CHIPS0132-1608','CHIPS2340-2302','CHIPS0018-1840','CHIPS0118-1430']

    for field in all_fields:
        pisco_photometry_psf_v4_nostar(field, mode=mode)
