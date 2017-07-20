import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from photutils import aperture_photometry
from photutils import SkyCircularAperture
from photutils import SkyCircularAnnulus

import os
import subprocess
import shlex
import sys
# edited 5/9/17
# ------


def write_ds9_region_sky(table, allcol, color, fname):
    """
    write_ds9_region_sky: write out a ds9 region file
    INPUT:
    - table: the table that we are interested in plotting
    - allcol: the column name for the celestial coordinates
    - color: the color for the region in ds9 (usually 'green', 'red', 'blue')
    - fname: the output filename (should be 'xxx.reg')
    OUTPUT:
    - create a fname reg file for ds9 to see on the celestial coordinates
    """
    dir = "region/"
    if not os.path.exists(dir):
        os.makedirs(dir)
    text_file = open(os.path.join(dir, fname), "w")
    text_file.write("# Region file format: CIAO version 1.0\n")
    for i in range(0, len(table)):
        ra = table[allcol][i].ra.to_string(unit=u.hr, sep=':')
        dec = table[allcol][i].dec.to_string(unit=u.degree, sep=':')
        text_file.write("j2000; circle(%s,%s,0.1') # color=%s\n" %
                        (ra, dec, color))
    text_file.close()


def sex_band(field, band):
    """
    sex_band: running sextractor in order to use its position for photometry (mainly using i band)
    INPUT:
    - field: our field of interest 'Field026'
    - band: our band of interest
        - in this software, we are only interested in 'i' band since it is the best band to get
        all the objects, but the function is flexible enough to be able to run for other bands
    OUTPUT:
    - create sextractor output file (e.g., posField026_i.fits) which have all the location of
    objects that we are interested in.
    """
    # seeing=Table.read('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/PISCO_Jan17_seeing.csv')
    # see=seeing[seeing['Field']==int(field[-3:])]['Seeing'][0]
    see=1.
    outname=os.path.join('final','pos%s_%s.fits'%(field,band))
    cmd="sex %s -c pisco_pipeline/config_slr.sex -CATALOG_NAME %s -SEEING_FWHM %s" % (os.path.join('final','coadd_c%s_%s.fits'%(field,band)),outname,str(see))
    print cmd
    sub=subprocess.check_call(shlex.split(cmd))
    return outname

def cut(x):
    if x>0:
        return x
    else:
        return 0.

def aperature_f(field, band, aper_rad=3., annu_in=5., annu_out=10.):
    """
    aperature_f: running "photutils: aperture_photometry" to get the number of count within
    a radius 'aper_rad', comparing with the background annulus with inner radius 'annu_in' and
    outer radius 'annu_out'
    note about the algorithm: since we already subtracted the background count so that on average,
    the background count is roughly zero, we calculate SNR = N*/sqtr(N*+npix(N_bg)) where N_bg~0.
    INPUT:
    - field: object of interset e.g., 'Field026'
    - band: str for the name of the band that we are interested in.
    - aper_rad: the circular radius (in arcsec) for the photometry
    - annu_in: the inner radius (in arcsec) of the annulus for the background counts
    - annu_out: the outer radius (in arcsec) of the annulus for the background counts
    OUTPUT:
    - a new table with the photometry for the number of count for a specific 'band'
    """
    hdu = fits.open("final/coadd_c%s_%s.fits" % (field, band))
    w = WCS("final/coadd_c%s_%s.fits" % (field, band))

    pos_filename = os.path.join('final', 'pos%s_%s.fits' % (field, 'i'))
    if not os.path.isfile(pos_filename):
        print 'sextracting the i band position for %s' % field
        pos_filename = sex_band(field, 'i')

    t6 = Table.read(pos_filename)
    t7 = t6[t6['CLASS_STAR']>0.9]   #only select stas for calibartion
    positions = SkyCoord(t7['ALPHA_J2000'], t7['DELTA_J2000'], frame='icrs')

    radii=np.arange(0.2,7,0.1)
    apertures = [SkyCircularAperture(positions, r=r*u.arcsec) for r in radii]
    bkg_apertures=[]
    for r in radii:
        bkg_apertures.append(SkyCircularAnnulus(positions, r_in=(r+2)*u.arcsec, r_out=(r+5)*u.arcsec))
    apertures.extend(bkg_apertures)

    phot_table = aperture_photometry(hdu[0], apertures)
    for i in range(len(radii)):
        bkg_per_px = phot_table['aperture_sum_%s'%str(i+len(radii))]/apertures[i+len(radii)].to_pixel(w).area()
        bkg_per_star = bkg_per_px*apertures[i].to_pixel(w).area()
        phot_table['final_%i'%i]=phot_table['aperture_sum_%s'%str(i)]-bkg_per_star

    table=phot_table[phot_table.colnames[-1-len(radii):-1]]
    dfd=table.to_pandas()
    dfd=dfd.applymap(cut)
    kronr=(dfd*radii).sum(axis=1)/dfd.sum(axis=1)
    kronr[np.isnan(kronr)]=1.  #remove all the zeros to be 1.

    for i in range(len(table)):
        RA=phot_table['celestial_center'].ra.value[i]*u.degree
        DEC=phot_table['celestial_center'].dec.value[i]*u.degree
        positions = SkyCoord(RA, DEC, frame='fk5')
        rad = 2.5*kronr[i]
        aperture=SkyCircularAperture(positions, r=rad*u.arcsec)
        bkg_aperture=SkyCircularAnnulus(positions, r_in=(rad+1)*u.arcsec, r_out=(rad+3)*u.arcsec)
        apertures=[aperture,bkg_aperture]
        phot_tabled = aperture_photometry(hdu[0], apertures)
        bkg_per_px = phot_tabled['aperture_sum_1']/apertures[1].to_pixel(w).area()
        bkg_per_star = bkg_per_px*apertures[0].to_pixel(w).area()
        phot_tabled['final_0']=cut(phot_tabled['aperture_sum_0']-bkg_per_star)
        phot_tabled['snr'] = phot_tabled['final_0'] / (np.sqrt(phot_tabled['final_0'])+1e-3)
        if i == 0:
            phot_tablei=phot_tabled
        else:
            phot_tablei.add_row(phot_tabled[0])

    phot_tablei["CLASS_STAR"] = t7["CLASS_STAR"]
    phot_tablei["FLAGS"] = t7["FLAGS"]
    phot_tablei['id']=np.arange(0,len(phot_tablei),1)
    # phot_table2=phot_tabled[phot_tabled['final_0']>0]

    for name in phot_tablei.colnames[1:]:
        phot_tablei.rename_column(name, name + '_%s' % band)
    return phot_tablei


def create_star_fits_slr(field, g, r, i, z):
    """
    create_star_fits_slr: convert the count to magnitude value and combine all four bands (g,r,i,z)
    together into one big table with magnitude and its uncertainty
    INPUT:
    - field: object of interset e.g., 'Field026'
    - g, r, i, z:
    OUTPUT:
    - a new table with celestial coordinates (ALPHA_J2000, DELTA_J2000) and magnitude from g, r, i, z band
    """
    pos_filename = os.path.join('final', 'pos%s_%s.fits' % (field, 'i'))
    if not os.path.isfile(pos_filename):
        print 'sextracting the i band position for %s' % field
        pos_filename = sex_band(field, 'i')
    t7 = Table.read(pos_filename)

    print 'joining g,r,i,z band together for %s...' % field
    total = join(join(join(i, g, keys='id'), r, keys='id'), z, keys='id')
    # total["CLASS_STAR"] = t7[total['id'] - 1]["CLASS_STAR"]
    # total["FLAGS"] = t7[total['id'] - 1]["FLAGS"]

    total=total[total['final_0_g']>0]; print 'after remove zero g', len(total)
    total=total[total['final_0_r']>0]; print 'after remove zero r', len(total)
    total=total[total['final_0_i']>0]; print 'after remove zero i', len(total)
    total=total[total['final_0_z']>0]; print 'after remove zero z', len(total)

    total['magg'] = -2.5 * np.array(map(np.log10, total['final_0_g']))
    total['magr'] = -2.5 * np.array(map(np.log10, total['final_0_r']))
    total['magi'] = -2.5 * np.array(map(np.log10, total['final_0_i']))
    total['magz'] = -2.5 * np.array(map(np.log10, total['final_0_z']))

    total['maggerr'] = 1.0857 / total["snr_g"]
    total['magrerr'] = 1.0857 / total["snr_r"]
    total['magierr'] = 1.0857 / total["snr_i"]
    total['magzerr'] = 1.0857 / total["snr_z"]

    total["ALPHA_J2000"] = np.array([i.ra.degree for i in total['celestial_center_i']])
    total["DELTA_J2000"] = np.array([i.dec.degree for i in total['celestial_center_i']])

    star = total[(total["CLASS_STAR_i"] > 0.9)]# & (total['magi']>-13.4)]
    #magnitude cut for saturated stars (pix~1400 for the whole aperature of radius 7.5 px => mag~-13.4)
    gal = total[total["CLASS_STAR_i"] < 0.75]

    # write_ds9_region_sky(star, "celestial_center_i",
    #                      'red', 'star_%s.reg' % field)
    # write_ds9_region_sky(gal, "celestial_center_i",
    #                      'green', 'gal_%s.reg' % field)

    starr = star[['id', 'ALPHA_J2000', 'DELTA_J2000', 'magg', 'magr', 'magi',
                  'magz', 'maggerr', 'magrerr', 'magierr', 'magzerr', 'CLASS_STAR_i','FLAGS_i']]
    # starr.meta['aperture_photometry_args'] = u"method='exact'"
    starr.meta[u'aperture_photometry_args'] = "nan"
    starr.write('star_%s.fits' % field, overwrite=True)
    print "writeout 'star_%s.fits' file for slr to run" % field
    return total[['id', 'ALPHA_J2000', 'DELTA_J2000', 'magg', 'magr', 'magi', 'magz',
                  'maggerr', 'magrerr', 'magierr', 'magzerr', 'CLASS_STAR_i','FLAGS_i']]


def slr_running(field, bigmacs="pisco_pipeline/big-macs-calibrate-master"):
    """
    slr_running: running SLR script from github.com/patkel/big-macs-calibrate to get a calibrated magnitude
    INPUT:
    - field: object of interset e.g., 'Field026'
    - bigmacs: the location for "big-macs-calibrate" directoty
    OUTPUT:
    - a new table with added columns with name MAG_g,...,MAGERR_g,...
    """
    infile = 'star_%s.fits' % field
    pyfile = os.path.join(bigmacs, 'fit_locus.py')
    cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 5 -l -r ALPHA_J2000 -d DELTA_J2000 -j --plot=PLOTS_%s" \
        % (pyfile, infile, os.path.join(bigmacs, "coadd_mag.columns"), field)
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))
    # if sub == 0:
    #     print 'finish slr_file and return starr.fits.offsets.list'
    # else:
    #     print 'slr stops working.'


def update_color(fname, table):
    """
    update_color: using the output from SLR, update to the correct magnitude
    INPUT:
    - fname: input file from SLR output (...offsets.list)
    - table: the table that we want to update the value (from column magg,etc to MAG_g,etc)
    OUTPUT:
    - a new table with added columns with name MAG_g,...,MAGERR_g,...
    """
    with open(fname) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    band = [x.split(' ')[0][-1] for x in content[5:-1]]
    corr = [float(x.split(' ')[1]) for x in content[5:-1]]
    ecorr = [float(x.split(' ')[3]) for x in content[5:-1]]
    print 'bands = ', band

    table['MAG_' + band[0]] = table['mag' + band[0]] + corr[0]
    table['MAG_' + band[1]] = table['mag' + band[1]] + corr[1]
    table['MAG_' + band[2]] = table['mag' + band[2]] + corr[2]
    table['MAG_' + band[3]] = table['mag' + band[3]] + corr[3]
    table['MAGERR_' + band[0]] = table['mag' + band[0] + 'err'] + ecorr[0]
    table['MAGERR_' + band[1]] = table['mag' + band[1] + 'err'] + ecorr[1]
    table['MAGERR_' + band[2]] = table['mag' + band[2] + 'err'] + ecorr[2]
    table['MAGERR_' + band[3]] = table['mag' + band[3] + 'err'] + ecorr[3]
    return table


"""
pisco_photometry: run pisco output data from pisco_combine to correct for the photometry of each object
and determine which objects are stars/galaxies.
The pipeline is a combination of SLR algorithm (cite: https://github.com/patkel/big-macs-calibrate)
and Photutils for photometry aperatures

V2: update the photometry to calculate Kron radius instead of using a fix radius for each objects

ARGUMENTS:
1. fieldname for object (e.g., 'Field027')

EXAMPLES:
python pisco_pipeline/pisco_photometry.py Field026
"""

if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    slrdir = 'slr_output'
    if not os.path.exists(slrdir):
        os.makedirs(slrdir)

    field = str(sys.argv[1])
    r = aperature_f(field, 'r')
    print 'finished r', len(r)
    i = aperature_f(field, 'i')
    print 'finished i', len(i)
    z = aperature_f(field, 'z')
    print 'finished z', len(z)
    g = aperature_f(field, 'g')
    print 'finished g', len(g)
    print 'len in g, r, i, z bands:', len(g), len(r), len(i), len(z)

    total = create_star_fits_slr(field, g, r, i, z)
    slr_running(field)

    ntotal = update_color('star_%s.fits.offsets.list' % field, total)
    ntotal.write(os.path.join(slrdir, 'ntotal_%s.csv' % field), overwrite=True)

    cmd = 'mv star_%s.fits slr_output/' % field
    try:
        sub = subprocess.check_call(shlex.split(cmd))
    except (ValueError, RuntimeError, TypeError, NameError):
        pass

    cmd = 'mv star_%s.fits.offsets.list slr_output/' % field
    try:
        sub = subprocess.check_call(shlex.split(cmd))
    except (ValueError, RuntimeError, TypeError, NameError):
        pass
