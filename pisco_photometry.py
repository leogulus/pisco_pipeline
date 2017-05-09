import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import aperture_photometry
from photutils import SkyCircularAperture
from photutils import SkyCircularAnnulus

import os
import subprocess
import shlex
import sys
##------
def write_ds9_region_sky(table,allcol,color,fname):
    dir="region/"
    if not os.path.exists(dir):
        os.makedirs(dir)
    text_file = open(os.path.join(dir,fname), "w")
    text_file.write("# Region file format: CIAO version 1.0\n")
    for i in range(0,len(table)):
        ra=table[allcol][i].ra.to_string(unit=u.hr, sep=':')
        dec=table[allcol][i].dec.to_string(unit=u.degree, sep=':')
        text_file.write("j2000; circle(%s,%s,0.1') # color=%s\n"  % (ra,dec,color))
    text_file.close()

def sex_band(field,band):
    outname=os.path.join('final','test%s_%s.fits'%(field,band))
    cmd="sex %s -c pisco_pipeline/config_slr.sex -CATALOG_NAME %s" % (os.path.join('final','coadd_c%s_%s.fits'%(field,band)),outname)
    print cmd
    sub=subprocess.check_call(shlex.split(cmd))
    return outname

def aperature_f(field,band,aper_rad=2.,annu_in=3.,annu_out=6.):

    hdu=fits.open("final/coadd_c%s_%s.fits" % (field,band))

    pos_filename=os.path.join('final','test%s_%s.fits'%(field,'i'))
    if not os.path.isfile(pos_filename):
        print 'sextracting the i band position for %s' % field
        pos_filename=sex_band(field,'i')

    t7 = Table.read(pos_filename)
    positions = SkyCoord(t7['XWIN_WORLD'], t7['YWIN_WORLD'], frame='icrs')
    apertures = SkyCircularAperture(positions, r=aper_rad*u.arcsec)
    annulus_apertures = SkyCircularAnnulus(positions, r_in=annu_in*u.arcsec, r_out=annu_out*u.arcsec)
    apers = [apertures, annulus_apertures]

    phot_table = aperture_photometry(hdu[0], apers)
    phot_table2 = phot_table[np.isfinite(phot_table['aperture_sum_0'])]
    phot_table3 = phot_table2[np.where(phot_table2['aperture_sum_0']>0)]

    snr=phot_table3['aperture_sum_0']/(np.sqrt(phot_table3['aperture_sum_0']))
    phot_table3['snr'] = snr

    for name in phot_table3.colnames[1:]:
        phot_table3.rename_column(name,name+'_%s' % band)
    return phot_table3

def create_star_fits_slr(field,g,r,i,z):
    pos_filename=os.path.join('final','test%s_%s.fits'%(field,'i'))
    if not os.path.isfile(pos_filename):
        print 'sextracting the i band position for %s' % field
        pos_filename=sex_band(field,'i')
    t7=Table.read(pos_filename)

    print 'joining g,r,i,z band together for %s...' % field
    total=join(join(join(i, g, keys='id'), r, keys='id'), z, keys='id')
    total["CLASS_STAR"]=t7[total['id']-1]["CLASS_STAR"]

    total['magg']=-2.5*np.array(map(np.log10,total['aperture_sum_0_g']))
    total['magr']=-2.5*np.array(map(np.log10,total['aperture_sum_0_r']))
    total['magi']=-2.5*np.array(map(np.log10,total['aperture_sum_0_i']))
    total['magz']=-2.5*np.array(map(np.log10,total['aperture_sum_0_z']))

    total['maggerr']=1.0857/total["snr_g"]
    total['magrerr']=1.0857/total["snr_r"]
    total['magierr']=1.0857/total["snr_i"]
    total['magzerr']=1.0857/total["snr_z"]

    total["XWIN_WORLD"]=total["celestial_center_i"].ra.degree
    total["YWIN_WORLD"]=total["celestial_center_i"].dec.degree

    star=total[total["CLASS_STAR"]>0.9]
    gal=total[total["CLASS_STAR"]<0.6]
    write_ds9_region_sky(star,"celestial_center_i",'red','star_%s.reg' % field)
    write_ds9_region_sky(gal,"celestial_center_i",'green','gal_%s.reg' % field)

    starr=star[['id','XWIN_WORLD','YWIN_WORLD','magg','magr','magi','magz','maggerr','magrerr','magierr','magzerr','CLASS_STAR']]
    starr.meta['aperture_photometry_args']=u"method='exact'"
    starr.write('star_%s.fits' % field, overwrite=True)
    print "writeout 'star_%s.fits' file for slr to run" % field
    return total[['id','XWIN_WORLD','YWIN_WORLD','magg','magr','magi','magz','maggerr','magrerr','magierr','magzerr','CLASS_STAR']]

def slr_running(field,bigmacs="pisco_pipeline/big-macs-calibrate-master"):
    infile='star_%s.fits' % field
    pyfile = os.path.join(bigmacs,'fit_locus.py')
    cmd="python %s --file %s --columns %s --extension 1 --bootstrap 5 -l -r XWIN_WORLD -d YWIN_WORLD -j --plot=PLOTS_%s" \
    % (pyfile,infile,os.path.join(bigmacs,"coadd_mag.columns"),field)
    print cmd
    sub=subprocess.check_call(shlex.split(cmd))
    if sub==0:
        print 'finish slr_file and return starr.fits.offsets.list'
    else:
        print 'slr stops working.'

def update_color(fname,table):
    with open(fname) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    band=[x.split(' ')[0][-1] for x in content[5:-1]]
    corr=[float(x.split(' ')[1]) for x in content[5:-1]]
    ecorr=[float(x.split(' ')[3]) for x in content[5:-1]]
    print 'bands = ',band

    table['MAG_'+band[0]]=table['mag'+band[0]]+corr[0]
    table['MAG_'+band[1]]=table['mag'+band[1]]+corr[1]
    table['MAG_'+band[2]]=table['mag'+band[2]]+corr[2]
    table['MAG_'+band[3]]=table['mag'+band[3]]+corr[3]
    table['MAGERR_'+band[0]]=table['mag'+band[0]+'err']+ecorr[0]
    table['MAGERR_'+band[1]]=table['mag'+band[1]+'err']+ecorr[1]
    table['MAGERR_'+band[2]]=table['mag'+band[2]+'err']+ecorr[2]
    table['MAGERR_'+band[3]]=table['mag'+band[3]+'err']+ecorr[3]
    return table

"""
pisco_photometry: run pisco output data from pisco_combine to correct for the photometry of each object and determine which objects are stars/galaxies.
The pipeline is a combination of SLR algorithm (cite: https://github.com/patkel/big-macs-calibrate) and Photutils for photometry aperatures

ARGUMENTS:
1. fieldname for object (e.g., 'Field027')

EXAMPLES:
python pisco_pipeline/pisco_photometry.py Field026
"""


if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    slrdir='slr_output'
    if not os.path.exists(slrdir):
        os.makedirs(slrdir)

    field=str(sys.argv[1])

    r=aperature_f(field,'r')
    i=aperature_f(field,'i')
    z=aperature_f(field,'z')
    g=aperature_f(field,'g')
    print 'len in g, r, i, z bands:', len(g), len(r), len(i), len(z)

    total=create_star_fits_slr(field,g,r,i,z)
    slr_running(field)
    ntotal=update_color('star_%s.fits.offsets.list'%field,total)
    ntotal.write(os.path.join(slrdir,'ntotal_%s.csv' % field), overwrite=True)

    cmd='mv star_%s.fits slr_output/' % field
    try:
        sub=subprocess.check_call(shlex.split(cmd))
    except (ValueError, RuntimeError, TypeError, NameError):
        pass

    cmd='mv star_%s.fits.offsets.list slr_output/' % field
    try:
        sub=subprocess.check_call(shlex.split(cmd))
    except (ValueError, RuntimeError, TypeError, NameError):
        pass
