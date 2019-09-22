import sys, os, re, yaml, subprocess, shlex, FITS_tools
import pandas as pd
import numpy as np

import pickle

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import image
import matplotlib.cm as cm
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy import interpolate
from scipy.interpolate import interp1d
import scipy.stats

from astropy.io import fits
from astropy.table import Table, join
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=71, Om0=0.3, Tcmb0=2.725)

import extra_program as ex
from PIL import Image as Image_PIL
import ebvpy  #Galactic Reddening

"""
Example: 
python pisco_pipeline/panstarr_photometry_all.py PKS1353 psf allslr 2mass
python pisco_pipeline/panstarr_photometry_all.py PKS1353 psf allslr no2mass
python pisco_pipeline/panstarr_photometry_all.py PKS1353 psf noslr no2mass
python pisco_pipeline/panstarr_photometry_all.py PKS1353 model noslr no2mass
python pisco_pipeline/panstarr_photometry_all.py PKS1353 auto slr 2mass

field: name of the fields
mode: psf, auto, aper, hybrid, model
allslr:
  - allslr: run everything including photometry_v4, cut_frame, SLR
  - slr: run just SLR and update the color
  - noslr: don't run slr, just update the color with different modes
2mass
  - 2mass: run SLR with 2MASS to match
  - no2mass: run SLR without 2MASS
"""

###--------------------------------------------------------------------------###
def find_seeing(field,band):
    df_see=pd.read_csv('/Users/taweewat/Documents/red_sequence/total_chips_field_seeing.csv',index_col=0)
    if field[0:5]=='CHIPS':
        seeing = df_see[df_see.chips==field]['seeing_q25_%s'%band].values[0] #_%s'%band
        return seeing
    elif (field[0:5]=='Field')|(field[0:3]=='PKS')|(field[0:4]=='SDSS'):
        seeing = df_see[df_see.name==field]['seeing_q25_%s'%band].values[0] #_%s'%band
        return seeing
def find_seeing_fits(field):
    home='/Users/taweewat/Documents/pisco_code/'
    dirs=['ut170103/','ut170104/','ut170619/','ut170621/','ut170624/','ut171208/','ut171209/','ut171212/']
    myReg=re.compile(r'(%s_A).*'%field)
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                seeing=float(fits.open(dir+myReg.search(text).group())[0].header['FWHM1'])
    seeing=0.5
    return seeing
def read_param():
    with open("pisco_pipeline/params.yaml", 'r') as stream:
        try:
            param=yaml.load(stream)
            return param
        except yaml.YAMLError as exc:
            print(exc)
def read_param_izp(mode):
    if mode=='psf':
        mode_izp=''
    elif mode=='model':
        mode_izp='' #'_model'
    else:
        mode_izp=''
    # print "/Users/taweewat/Documents/pisco_code/pisco_pipeline/params_izeropoint%s.yaml" % mode_izp
    with open("/Users/taweewat/Documents/pisco_code/pisco_pipeline/params_izeropoint%s.yaml"%mode_izp, 'r') as stream:
        try:
            param=yaml.load(stream)
            return param
        except yaml.YAMLError as exc:
            print(exc)

def star_galaxy_bleem(field):
    sg_dir = 'star_galaxy'
    if not os.path.exists(sg_dir):
        os.makedirs(sg_dir)

    param=read_param()
    seeing=find_seeing(field,'i')
    # seeing=1.5
    # seeing=0.95
    minarea=1.7
    data, header = fits.getdata('final/coadd_c%s_i.fits'%field, header=True)
    data2=data**2
    pxscale=0.22
    fits.writeto('final/coadd_c%s_sq_i.fits'%field, data2, header=header, overwrite=True)
    cmd='sex final/coadd_c%s_i.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -CATALOG_TYPE FITS_1.0 -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE %s -DETECT_MINAREA %s -CHECKIMAGE_NAME checki.fits,segmenti.fits'%\
    (field,'sex.param',sg_dir+'/%s_catalog.fits'%(field),str(seeing),str(param['satur_level_i_psf']),str(pxscale),str(1.1/minarea*np.pi*(seeing/pxscale)**2)); print cmd
    subprocess.check_call(shlex.split(cmd))
    cmd='sex final/coadd_c%s_i.fits,final/coadd_c%s_sq_i.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -CATALOG_TYPE FITS_1.0 -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE %s -DETECT_MINAREA %s'%\
    (field,field,'sex.param',sg_dir+'/%s_sq_catalog.fits'%(field),str(seeing),str(param['satur_level_i_sq_psf']),str(pxscale),str(1.1/minarea*np.pi*(seeing/pxscale)**2)); print cmd
    subprocess.check_call(shlex.split(cmd))

def pisco_photometry_v4(field):
    def aperature_proj(field,band):
        param=read_param()
        # seeing=find_seeing(field,band)
        # seeing=1.5
        seeing=0.5
        seeing_class=1.8
        # saturation=9.0
        saturation=54000.

        data, header = fits.getdata('/Users/taweewat/Documents/red_sequence/panstar/coadd_panstar_{}_{}.fits'.format(field,band), header=True)
        data2=data*6000.
        fits.writeto('/Users/taweewat/Documents/red_sequence/panstar/coadd_scaled_panstar_{}_{}.fits'.format(field,band), data2, header=header, overwrite=True)

        slrdir = 'slr_output'
        # to_be_projected = '/Users/taweewat/Documents/red_sequence/panstar/coadd_panstar_{}_{}.fits'.format(field,band)
        to_be_projected = '/Users/taweewat/Documents/red_sequence/panstar/coadd_scaled_panstar_{}_{}.fits'.format(field,band)
        reference_fits  = '/Users/taweewat/Documents/red_sequence/panstar/coadd_panstar_{}_i.fits'.format(field)
        im1,im2, header = FITS_tools.match_fits(to_be_projected,reference_fits,return_header=True)
        # outname = 'final/proj_coadd_panstar_%s_%s.fits'%(field,band)
        outname = 'final/proj_coadd_panstar_%s_%s.fits'%(field,band)
        print 'projecting from %s band to i band the fits file '%band + outname
        fits.writeto(outname, im1, header, overwrite=True)

        minarea=1.7 #1.7
        pxscale=0.25
        cmd='sex /Users/taweewat/Documents/red_sequence/panstar/coadd_scaled_panstar_%s_%s.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 23 -PIXEL_SCALE %s -DETECT_MINAREA %s -CHECKIMAGE_NAME check_panstar_psf_%s.fits,segment_panstar_psf_%s.fits'%\
        (field,band,'sex_fwhm_psf.param','psfex_output/psf_%s_%s.fits'%(field,band),str(seeing_class),str(saturation),str(pxscale),str(1.1/minarea*np.pi*(seeing/pxscale)**2), band, band) 
        print cmd
        subprocess.check_call(shlex.split(cmd))

        Tf=Table(fits.open('psfex_output/psf_%s_%s.fits'%(field,band))[2].data)
        Tf=Tf[(Tf['FLUX_APER']>0)]

        df0 = pd.read_csv('/Users/taweewat/Documents/red_sequence/{}_star_list.csv'.format(field),index_col=0)
        x=np.array([286.0227455650082,285.9411038202907,286.0569138078614,285.9817436730952,286.00556207826133,286.01921620713756])
        real=np.array([287.7272544,287.8089409,287.6931021,287.7682687,287.7444376,287.7307988])
        p=np.poly1d(np.polyfit(x,real,1))
        # c0 = SkyCoord(ra=df0['raMean'].values*u.degree, dec=df0['decMean'].values*u.degree)
        # csex = SkyCoord(ra=p(np.array(Tf['ALPHA_J2000']))*u.degree, dec=np.array(Tf['DELTA_J2000'])*u.degree)
        # idxn, d2dn, d3dn=csex.match_to_catalog_sky(c0)

        # Tfcut=Tf[d2dn.to(u.arcsecond).value<2]
        Tfcut=Tf
        print "len of Tfcut after 2 arcsecond: {}".format(len(Tfcut))

        vignet_bad=[]
        for i in range(len(Tfcut)):
            vignet_bad.append(np.sum(Tfcut['VIGNET'][i].ravel()<-9e+29))
        Tfcut['VIGNET_bad']=vignet_bad

        # Tfcut=Tfcut[(Tfcut['FLAGS'] == 0) & (Tfcut['VIGNET_bad'] < 20)]# & (Tfcut['FLUX_APER'] > 300)].copy()
        Tfcut = Tfcut[(Tfcut['VIGNET_bad'] < 20)].copy() #(Tfcut['CLASS_STAR'] > 0.70) & (Tfcut['FLAGS'] < 4) & 
        Tfcut = Tfcut[(-2.5*np.log10(Tfcut['FLUX_APER'])<-16.)] #
        Tfcut_edge=Tfcut#[(Tfcut['XWIN_IMAGE']<np.max(Tfcut['XWIN_IMAGE'])-60)&(Tfcut['XWIN_IMAGE']>np.min(Tfcut['XWIN_IMAGE'])+60)&\
                #(Tfcut['YWIN_IMAGE']<np.max(Tfcut['YWIN_IMAGE'])-60)&(Tfcut['YWIN_IMAGE']>np.min(Tfcut['YWIN_IMAGE'])+60)].copy()
        Tfcut_more=Tfcut_edge[(np.abs(Tfcut_edge['FLUX_RADIUS']-np.mean(Tfcut_edge['FLUX_RADIUS']))<2*np.std(Tfcut_edge['FLUX_RADIUS']))]
        Tfcut_more2=Tfcut_more[(np.abs(Tfcut_more['ELONGATION']-np.mean(Tfcut_more['ELONGATION']))<2*np.std(Tfcut_more['ELONGATION']))].copy()
        print "length of Tf: all: {}, cut: {}, edges: {}, flux_radius: {}, elong: {}".format(len(Tf), len(Tfcut), len(Tfcut_edge), len(Tfcut_more), len(Tfcut_more2))
        hdu = fits.open('psfex_output/psf_%s_%s.fits'%(field,band))
        hdu[2].data = hdu[2].data[Tfcut_more2['NUMBER']-1]
        # hdu[2].data = hdu[2].data[Tfcut['NUMBER']-1]
        hdu.writeto('psfex_output/psf_%s_%s.fits'%(field,band), overwrite=True)

        cmd='psfex %s -c pisco_pipeline/panstarr.psfex' % ('psfex_output/psf_%s_%s.fits'%(field,band))
        print cmd
        subprocess.check_call(shlex.split(cmd))

        cmd='sex /Users/taweewat/Documents/red_sequence/panstar/coadd_scaled_panstar_%s_i.fits,final/proj_coadd_panstar_%s_%s.fits -c pisco_pipeline/config.sex -PSF_NAME %s -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PIXEL_SCALE %s -CATALOG_TYPE FITS_1.0 -PHOT_APERTURES 23 -DETECT_MINAREA %s -CHECKIMAGE_NAME check%s.fits,segment%s.fits'%\
            (field, field, band, 'psfex_output/psf_%s_%s.psf' % (field, band), 'sex_after_psf.param', '%s/a_psf_%s_%s.fits' % (slrdir, field, band),
             str(seeing_class), str(saturation), str(pxscale), str(1.1 / minarea * np.pi * (seeing / pxscale)**2), band, band)
        print cmd
        subprocess.check_call(shlex.split(cmd))

        table=Table.read('%s/a_psf_%s_%s.fits'%(slrdir,field,band))
        table['ALPHA_J2000']=p(np.array(table['ALPHA_J2000']))
        for name in table.colnames[:]:
            table.rename_column(name, name + '_%s' % band)
        return table

        # return Tf
    
    slrdir = 'slr_output'
    if not os.path.exists(slrdir):
        os.makedirs(slrdir)

    tableg=aperature_proj(field,'g')
    tablei=aperature_proj(field,'i')
    tabler=aperature_proj(field,'r')
    # tablez=aperature_proj(field,'z')

    print 'len of all table: {}, {}, {}'.format(len(tableg), len(tablei), len(tabler))

    ci=SkyCoord(ra=np.array(tablei['ALPHA_J2000_i'])*u.degree, dec=np.array(tablei['DELTA_J2000_i'])*u.degree)# print len(ci)
    cg=SkyCoord(ra=np.array(tableg['ALPHA_J2000_g'])*u.degree, dec=np.array(tableg['DELTA_J2000_g'])*u.degree)# print len(cg)
    cr=SkyCoord(ra=np.array(tabler['ALPHA_J2000_r'])*u.degree, dec=np.array(tabler['DELTA_J2000_r'])*u.degree)# print len(cr)
    # cz=SkyCoord(ra=np.array(tablez['ALPHA_J2000_z'])*u.degree, dec=np.array(tablez['DELTA_J2000_z'])*u.degree)# print len(cz)

    idxn, d2dn, d3dn=cg.match_to_catalog_sky(ci)
    Table_I=tablei[idxn][['NUMBER_i','XWIN_IMAGE_i','YWIN_IMAGE_i','ALPHA_J2000_i','DELTA_J2000_i','MAG_APER_i','MAGERR_APER_i','MAG_AUTO_i','MAGERR_AUTO_i','MAG_SPHEROID_i','MAGERR_SPHEROID_i',\
                  'CLASS_STAR_i','FLAGS_i','MAG_PSF_i','MAGERR_PSF_i','MAG_MODEL_i','MAGERR_MODEL_i','SPREAD_MODEL_i','SPREADERR_MODEL_i']]
    Table_I.rename_column('ALPHA_J2000_i','ALPHA_J2000')
    Table_I.rename_column('DELTA_J2000_i','DELTA_J2000')

    idxn, d2dn, d3dn=cg.match_to_catalog_sky(cr)
    Table_R=tabler[idxn][['NUMBER_r','ALPHA_J2000_r','DELTA_J2000_r','MAG_APER_r','MAGERR_APER_r','MAG_AUTO_r','MAGERR_AUTO_r','MAG_SPHEROID_r','MAGERR_SPHEROID_r',\
                  'CLASS_STAR_r','FLAGS_r','MAG_PSF_r','MAGERR_PSF_r','MAG_MODEL_r','MAGERR_MODEL_r','SPREAD_MODEL_r','SPREADERR_MODEL_r']]
    Table_R.rename_column('ALPHA_J2000_r','ALPHA_J2000')
    Table_R.rename_column('DELTA_J2000_r','DELTA_J2000')

    Table_G = tableg[['NUMBER_g', 'ALPHA_J2000_g', 'DELTA_J2000_g', 'MAG_APER_g', 'MAGERR_APER_g', 'MAG_AUTO_g', 'MAGERR_AUTO_g', 'MAG_SPHEROID_g', 'MAGERR_SPHEROID_g',
                  'CLASS_STAR_g','FLAGS_g','MAG_PSF_g','MAGERR_PSF_g','MAG_MODEL_g','MAGERR_MODEL_g','SPREAD_MODEL_g','SPREADERR_MODEL_g']]
    Table_G.rename_column('ALPHA_J2000_g','ALPHA_J2000')
    Table_G.rename_column('DELTA_J2000_g','DELTA_J2000')

    print 'len of all new table', len(Table_G), len(Table_I), len(Table_R)

    total=join(join(Table_I,Table_G,keys=['ALPHA_J2000','DELTA_J2000']),Table_R,keys=['ALPHA_J2000','DELTA_J2000'])

    total.write(os.path.join(slrdir, 'total0_psf_%s.csv' % field), overwrite=True)

    total2=total[['ALPHA_J2000','DELTA_J2000','NUMBER_i','NUMBER_r','NUMBER_g','XWIN_IMAGE_i','YWIN_IMAGE_i',\
                'MAG_APER_i','MAGERR_APER_i','MAG_APER_g','MAGERR_APER_g','MAG_APER_r','MAGERR_APER_r','MAG_AUTO_i',\
                'MAGERR_AUTO_i','MAG_AUTO_g','MAGERR_AUTO_g','MAG_AUTO_r','MAGERR_AUTO_r','MAG_SPHEROID_i',\
                'MAGERR_SPHEROID_i','MAG_SPHEROID_g','MAGERR_SPHEROID_g','MAG_SPHEROID_r','MAGERR_SPHEROID_r',\
                'CLASS_STAR_i','CLASS_STAR_g','CLASS_STAR_r','FLAGS_g','FLAGS_r','FLAGS_i','MAG_PSF_g',\
                'MAG_PSF_r','MAG_PSF_i','MAGERR_PSF_g','MAGERR_PSF_r','MAGERR_PSF_i','MAG_MODEL_g','MAG_MODEL_r',\
                'MAG_MODEL_i','MAGERR_MODEL_g','MAGERR_MODEL_r','MAGERR_MODEL_i','SPREAD_MODEL_g','SPREAD_MODEL_r',\
                'SPREAD_MODEL_i','SPREADERR_MODEL_g','SPREADERR_MODEL_r','SPREADERR_MODEL_i']]

    total2.write(os.path.join(slrdir, 'total_psf_%s.csv' % field), overwrite=True)

def pisco_cut_star(field,c_a,c_b,c_d,c_delta):
    seeing=find_seeing_fits(field)
    true_seeing=find_seeing(field,'i')

    df_i=Table(fits.open('/Users/taweewat/Documents/pisco_code/star_galaxy/%s_catalog.fits'%field)[1].data).to_pandas()
    df_isq=Table(fits.open('/Users/taweewat/Documents/pisco_code/star_galaxy/%s_sq_catalog.fits'%field)[1].data).to_pandas()

    #cut the object out so that it has the same number of object between the sq catalog list and the psf mag list.
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/total_psf_%s.csv"%field
    df0 = pd.read_csv(fname)
    df0['NUMBER'] = np.arange(0, len(df0), 1).tolist()
    cf_i=SkyCoord(ra=np.array(df_i['ALPHA_J2000'])*u.degree, dec=np.array(df_i['DELTA_J2000'])*u.degree)
    cf_isq=SkyCoord(ra=np.array(df_isq['ALPHA_J2000'])*u.degree, dec=np.array(df_isq['DELTA_J2000'])*u.degree)
    cf0=SkyCoord(ra=np.array(df0['ALPHA_J2000'])*u.degree, dec=np.array(df0['DELTA_J2000'])*u.degree)
    df0.rename(columns={'ALPHA_J2000': 'ALPHA_J2000_i'}, inplace=True)
    df0.rename(columns={'DELTA_J2000': 'DELTA_J2000_i'}, inplace=True)

    idxn, d2dn, d3dn=cf0.match_to_catalog_sky(cf_i)
    df_i_cut0=df_i.loc[idxn].copy()
    df_i_cut0['NUMBER']=np.arange(0,len(df0),1).tolist()
    df_i_cut=pd.merge(df_i_cut0,df0,on='NUMBER')

    idxn, d2dn, d3dn=cf0.match_to_catalog_sky(cf_isq)
    df_isq_cut0=df_isq.loc[idxn].copy()
    df_isq_cut0['NUMBER']=np.arange(0,len(df0),1).tolist()
    df_isq_cut=pd.merge(df_isq_cut0,df0,on='NUMBER')

    fig,ax=plt.subplots(2,3,figsize=(15,10))
    df_i0=df_i_cut[(df_i_cut.MAG_APER<0)&(df_isq_cut.MAG_APER<0)]
    df_isq0=df_isq_cut[(df_i_cut.MAG_APER<0)&(df_isq_cut.MAG_APER<0)]# print len(df_i), len(df_isq)

    # c_d=-7.5
    df_i2=df_i0[(df_i0.CLASS_STAR>c_a) & (df_i0.MAG_APER<c_d)]# & (df_i0.MAG_APER>c_c)]
    df_isq2=df_isq0[(df_i0.CLASS_STAR>c_a) & (df_i0.MAG_APER<c_d)]# & (df_i0.MAG_APER>c_c)];# print len(df_i2), len(df_isq2)

    icut_per=np.percentile(df_i2.MAG_APER,35) #35
    df_i3=df_i2[df_i2.MAG_APER>icut_per]
    df_isq3=df_isq2[df_i2.MAG_APER>icut_per]

    fit=np.polyfit(df_i3.MAG_APER, df_i3.MAG_APER-df_isq3.MAG_APER, 1)
    f=np.poly1d(fit)
    ax[0,0].plot(df_i2.MAG_APER,f(df_i2.MAG_APER),'--')

    res=(df_i3.MAG_APER-df_isq3.MAG_APER)-f(df_i3.MAG_APER)
    aa=np.abs(res)<1.5*np.std(res)
    # outl=np.abs(res)>=1.5*np.std(res)
    fit=np.polyfit(df_i3.MAG_APER[aa], df_i3.MAG_APER[aa]-df_isq3.MAG_APER[aa], 1)
    f=np.poly1d(fit)

    ax[0,0].axvline(icut_per,color='blue',label='35th quantile')
    ax[0,0].errorbar(df_i2.MAG_APER,df_i2.MAG_APER-df_isq2.MAG_APER,yerr=np.sqrt(df_i2.MAGERR_APER**2+df_isq2.MAGERR_APER**2),fmt='o')
    ax[0,0].set_title('only for star')
    ax[0,0].plot(df_i2.MAG_APER,f(df_i2.MAG_APER),'--',label='no outlier')
    ax[0,0].set_ylabel('MAG_APER-MAG_APER_sq')
    ax[0,0].set_xlabel('MAG APER i')

    #--->  #0.1 default, 0.2
    c_c=df_i2[f(df_i2.MAG_APER)-(df_i2.MAG_APER-df_isq2.MAG_APER)<0.1]['MAG_APER'].values\
    [np.argmin(df_i2[f(df_i2.MAG_APER)-(df_i2.MAG_APER-df_isq2.MAG_APER)<0.1]['MAG_APER'].values)] #edit10/30 (previous 0.1)
    #--->

    ax[0,0].axvline(c_c,color='red',label='new upper cut')
    ax[0,0].legend(loc='best')

    # color_axis='CLASS_STAR'
    color_axis='SPREAD_MODEL_i'


    ax[0,1].scatter(df_i0.MAG_APER,df_i0.MAG_APER-df_isq0.MAG_APER,marker='.',c=df_i0[color_axis],vmin=0., vmax=0.005)
    ax[0,1].plot(df_i3.MAG_APER,df_i3.MAG_APER-df_isq3.MAG_APER,'x')
    ax[0,1].set_title('for all objects')
    ax[0,1].set_ylabel('MAG_APER-MAG_APER_sq')
    ax[0,1].set_xlabel('MAG APER i')
    ax[0,1].axvline(c_b,ls='--')
    ax[0,1].axvline(c_c,ls='--')

    delta=(df_i0.MAG_APER-df_isq0.MAG_APER) - f(df_i0.MAG_APER)

    ax[0,2].scatter(df_i0.MAG_APER,delta,marker='.',c=df_i0[color_axis],vmin=0., vmax=0.005)
    ax[0,2].axhline(0,ls='--')
    ax[0,2].axvline(c_c,ls='--')
    ax[0,2].axvline(c_b,ls='--')
    ax[0,2].set_ylabel('Delta')
    ax[0,2].set_xlabel('MAG APER i')
    ax[0,2].set_ylim(0.5,-1.2)

    df_i1=df_i0[(df_i0.MAG_APER>c_c)&(df_i0.MAG_APER<c_b)].copy()
    df_isq1=df_isq0[(df_i0.MAG_APER>c_c)&(df_i0.MAG_APER<c_b)].copy()
    delta1=(df_i1.MAG_APER-df_isq1.MAG_APER) - f(df_i1.MAG_APER)
    ax[1,0].scatter(df_i1.MAG_APER, delta1, marker='o', c=df_i1[color_axis],vmin=0., vmax=0.005)
    ax[1,0].axhline(0,ls='--')
    ax[1,0].axhline(c_delta, ls='--')
    ax[1,0].set_ylabel('Delta')
    ax[1,0].set_xlabel('MAG APER i')
    ax[1,0].set_ylim(0.5,-2)

    # deltag=delta1[delta1<c_delta] #galaxy  0.1, 0.2 (0.005), 0.5 ()
    deltas=delta1[(delta1>=c_delta)&(delta1<3.)] #star

    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    p0 = [1., 0., 0.1]
    # def gauss(x, *p):
    #     A, sigma = p
    #     return A*np.exp(-(x-0)**2/(2.*sigma**2))
    # p0 = [1., 0.1]

    #galaxy
    # hist, bin_edges = np.histogram(deltag,bins=np.arange(-1.2,0.5,0.02))
    hist, bin_edges = np.histogram(delta1,bins=np.arange(-1.2,0.5,0.02))
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    ax[1,1].plot(bin_centres, hist, label='galaxies',linestyle='steps')

    #stars
    hist, bin_edges = np.histogram(deltas,bins=np.arange(-1,0.5,0.02)) #(0 vs -1,0.5,0.02)
    # hist, bin_edges = np.histogram(delta1, bins=np.arange(c_delta, 0.5, 0.02))
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    coeff2, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    ax[1,1].plot(bin_centres, hist, label='stars',linestyle='steps')

    # hist, bin_edges = np.histogram(delta1,bins=np.arange(-1.2,0.5,0.02)) #added for right gaussian fitting
    # bin_centres = (bin_edges[:-1] + bin_edges[1:])/2  # added for right gaussian fitting
    x=np.arange(-1.25,0.5,0.02)
    # hist_fit2 = gauss(x, *coeff2)
    hist_fit2 = gauss(x, *coeff2)
    hist_fit3 = gauss(x, *coeff2)/np.max(gauss(x, *coeff2)) #added for right gaussian fitting
    ax[1,1].plot(x, hist_fit2, label='stars_fit')
    ax[1,1].plot(x, hist_fit3, label='stars_fit_norm') #added for right gaussian fitting
    ax[1,1].axvline(x[hist_fit3>star_cut][0],c='tab:pink',label='cut:%.3f'%x[hist_fit3>star_cut][0])  #added for right gaussian fitting
    ax[1,1].legend(loc='best')
    ax[1,1].set_xlabel('Delta')
    ax[1,1].set_ylabel('Histogram')

    ax[0,2].axhline(x[hist_fit3>star_cut][0],c='tab:pink') #added for right gaussian fitting
    ax[1,0].axhline(x[hist_fit3>star_cut][0],c='tab:pink') #added for right gaussian fitting
    ax[1,2].axhline(star_cut, c='tab:red')  # added for right gaussian fitting

    maxi=np.max(gauss(delta,*coeff2))
    def prob_SG(delta,maxi,*coeff2):
        if delta>0.:
            return 0.
        elif delta<=0.:
            return 1. - (gauss(delta, *coeff2) / maxi)

    vprob_SG= np.vectorize(prob_SG)
    SG=1.-vprob_SG(delta1,maxi,*coeff2)
    df_i1.loc[:,'SG']=SG
    param_izp=read_param_izp('psf')
    mag0=param_izp['i_zp_day9']#%dir_dict[find_fits_dir(field)[-9:]]]
    axi = ax[1, 2].scatter(df_i1.MAG_APER + mag0, SG,
                           marker='.', c=df_i1[color_axis], vmin=0., vmax=0.005)
    ax[1,2].axvline(aper_cut, ls='--', c='tab:blue')
    ax[1,2].axhline(SG_upper, ls='--', c='tab:blue')
    ax[1,2].set_ylim(-0.02,1.02)
    ax[1,2].set_xlabel('MAG APER i')
    ax[1,2].set_ylabel('SG (probability to be a star)')
    plt.suptitle(field+' seeing vs true_seeing: '+str(seeing)+','+str(true_seeing))
    fig.colorbar(axi)
    plt.tight_layout(rect=[0, 0., 1, 0.98])
    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/star_galaxy_sep_12_all%s.png' % field, dpi=120)
    plt.close(fig)
    return df_i_cut, df_i1

def pisco_cut_frame(field):
    # df_i=Table(fits.open('/Users/taweewat/Documents/pisco_code/star_galaxy/'+
    #                      '%s_catalog.fits'%field)[1].data).to_pandas()
    """
    c_a: CLASS_STAR lower limit for stars used for the linear fit
    c_b, c_c: upper and lower limit for all objects selection
    c_c can be moved with the for loop to include more objects until the confusion limit
    c_d: Faintest magnitude for stars used for the linear fit
    c_delta: lower limit for Delta to consider stars before fitting the gaussian and find SG (Star/Galaxy) factor 
    """

    global star_cut
    global aper_cut
    global SG_upper
    star_cut=0.95
    aper_cut=21.5
    # aper_cut=22.0
    SG_upper=0.02
    seeing=find_seeing_fits(field)
    true_seeing=find_seeing(field,'i')
    if field=='Field179':
        true_seeing=1.12
 
    if field=='CHIPS1011-0505':
        # c_a,c_b,c_c,c_d,c_delta=[0.95,-9.,-12.5,-9.2,-0.25]
        c_a,c_b,c_c,c_d,c_delta=[0.95,-12.,-12.5,-13,-0.1]
    elif field=='CHIPS2317-1443':
        c_a,c_b,c_c,c_d,c_delta=[0.95,-11.,-12.5,-14,-0.05]
    elif (field == 'Field137') or (field == 'Field071') or (field == 'Field109'):
        # c_a,c_b,c_c,c_d,c_delta=[0.95,-8.,-11.5,-8.5,-0.2]
        c_a,c_b,c_c,c_d,c_delta=[0.95,-12.,-12.5,-13,-0.15]
    elif (field[0:3]=='PKS') or (field[0:4]=='SDSS'):
        # c_a,c_b,c_c,c_d,c_delta=[0.95,-7,-12.2,-8.5,-0.1]
        c_a, c_b, c_c, c_d, c_delta = [0.95, -11, -13, -13, -0.1]
    elif field[0:5]=='CHIPS':
        # c_a,c_b,c_c,c_d,c_delta=[0.95,-9.,-12.5,-11.8,-0.1] #-0.15 (default), -10
        c_a,c_b,c_c,c_d,c_delta=[0.95,-11.,-12.2,-12.2,-0.1] #-0.15 (default), -10
    elif field[0:5]=='Field':
        # c_a,c_b,c_c,c_d,c_delta=[0.95,-8.,-11.5,-8.5,-0.1] #-8.5, -0.1 (Default)
        c_a,c_b,c_c,c_d,c_delta=[0.95,-11,-12.5,-12.5,-0.1]

    df_i_cut, df_i1=pisco_cut_star(field,c_a,c_b,c_d,c_delta)
    while len(df_i1[(df_i1['MAG_APER']<c_b)&(df_i1['MAG_APER']>c_b-0.5)\
                    &(df_i1['SG']>0.2)&(df_i1['SG']<0.8)])<25: #8, 10 (default)
        len_df_i1=len(df_i1)
        c_b=c_b+0.5
        df_i_cut, df_i1 = pisco_cut_star(field,c_a,c_b,c_d,c_delta)
        if len_df_i1==len(df_i1):
            break

    def SG_cut(SG, aper, aper0):
        if aper<aper0:
            if SG <= SG_upper:
                return True
            else:
                return False
        else:
            return True
    param_izp=read_param_izp('psf')
    mag0=param_izp['i_zp_day9']#%dir_dict[find_fits_dir(field)[-9:]]]
    vSG_cut = np.vectorize(SG_cut)
    vSGcut = vSG_cut(df_i1['SG'].values, df_i1['MAG_APER'].values+mag0, aper_cut)
    dff=df_i1[vSGcut]
    # dff=df_i1[df_i1['CLASS_STAR']<0.5]
    # dff=df_i1[df_i1['SG']<star_cut] #0.8

    vSGcut_star = vSG_cut(0.01, df_i1['MAG_APER'].values+mag0, aper_cut)
    # dff_star = df_i1[(df_i1['SG'] > 0.9)]  # no_2mass cut

    ##After running pisco_photometry_v4.fits, but before running
    ##pisco_photometry_psf_v4.fits, and then 19_pisco_tilt_resequence.ipynb
    # fname = "/Users/taweewat/Documents/pisco_code/slr_output/total_psf_%s.csv"%field
    # df0 = pd.read_csv(fname)
    # df0['NUMBER']=np.arange(0,len(df0),1).tolist()  #add NUMBER parameter to match between cut catalog and the total catalog
    # df0.rename(columns={'ALPHA_J2000': 'ALPHA_J2000_i'}, inplace=True)
    # df0.rename(columns={'DELTA_J2000': 'DELTA_J2000_i'}, inplace=True)
    # print len(df0), len(df_i_cut), len(dff), len(dff_star), '=', seeing, 'vs [NEW]', true_seeing
    # if len(df0)!=len(df_i_cut):
    #     raise ValueError('the two tables do not have the same length')
    # dff0=pd.merge(dff,df0,on='NUMBER')

    ##Using SPREAD_MODEL to seperate star/galaxies
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/total_psf_%s.csv"%field
    df0 = pd.read_csv(fname)
    df0['NUMBER'] = np.arange(0, len(df0), 1).tolist()
    df0.rename(columns={'ALPHA_J2000': 'ALPHA_J2000_i'}, inplace=True)
    df0.rename(columns={'DELTA_J2000': 'DELTA_J2000_i'}, inplace=True)
    #EXTENDED_COADD: 0 star, 1 likely star, 2 mostly galaxies, 3 galaxies
    # df0['EXTENDED_COADD']=np.array(((df0['SPREAD_MODEL_i']+ 3*df0['SPREADERR_MODEL_i'])>0.005).values, dtype=int)+\
    # np.array(((df0['SPREAD_MODEL_i']+df0['SPREADERR_MODEL_i'])>0.003).values, dtype=int)+\
    # np.array(((df0['SPREAD_MODEL_i']-df0['SPREADERR_MODEL_i'])>0.003).values, dtype=int)
    # dff=df0[df0['EXTENDED_COADD']>1]
    # dff_star=df0[df0['EXTENDED_COADD']<2]

    dff=df0[(df0['SPREAD_MODEL_i'])>0.005]
    dff_star=df0[np.abs(df0['SPREAD_MODEL_i'])<0.004]
    # dff_star=df0[np.abs(df0['SPREAD_MODEL_i'])<0.004] #+5/3.*df0['SPREADERR_MODEL_i'] <0.002

    dff0=dff
    dff0.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/"+\
                "galaxy_psf_total_%s.csv"%field)
    # dff_star0=pd.merge(dff_star, df0, on='NUMBER') # for non-SPREAD_MODEL
    dff_star0=dff_star #for SPREAD_MODEL
    dff_star0.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/"+\
                     "star_psf_total_%s.csv"%field)

def panstar_cut_star(field):
    ##Using SPREAD_MODEL to seperate star/galaxies
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/total_psf_%s.csv"%field
    df0 = pd.read_csv(fname)
    df0['NUMBER'] = np.arange(0, len(df0), 1).tolist()
    df0.rename(columns={'ALPHA_J2000': 'ALPHA_J2000_i'}, inplace=True)
    df0.rename(columns={'DELTA_J2000': 'DELTA_J2000_i'}, inplace=True)

    dfi=df0[df0['MAG_AUTO_i']<-16]
    x=dfi['MAG_AUTO_i']
    y=dfi['SPREAD_MODEL_i']
    p_spread=np.poly1d(np.polyfit(x,y,1))
    xs=np.arange(np.min(df0['MAG_AUTO_i']),np.max(df0['MAG_AUTO_i']),0.01)

    fig=plt.figure(figsize=(8,4))
    plt.subplot(1,2,1)
    plt.plot(df0['MAG_AUTO_i'],df0['SPREAD_MODEL_i'],'.',alpha=0.5)
    plt.plot(x,y,'.',alpha=0.5)
    plt.plot(xs,p_spread(xs))
    plt.axhline(0.005,color='tab:orange')
    plt.ylim(-0.1, 0.1)

    plt.subplot(1,2,2)
    plt.plot(df0['MAG_AUTO_i'],df0['SPREAD_MODEL_i']-p_spread(df0['MAG_AUTO_i']),'.',alpha=0.5)
    plt.axhline(0.005,color='tab:orange')
    plt.ylim(-0.1, 0.1)
    plt.tight_layout()
    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/spread_model_i_fit_%s_%s.png' %
            (mode, field), dpi=120)
    plt.close(fig)

    df0['SPREAD_MODEL_i2']=df0['SPREAD_MODEL_i']-p_spread(df0['MAG_AUTO_i'])

    #EXTENDED_COADD: 0 star, 1 likely star, 2 mostly galaxies, 3 galaxies
    # df0['EXTENDED_COADD']=np.array(((df0['SPREAD_MODEL_i']+ 3*df0['SPREADERR_MODEL_i'])>0.005).values, dtype=int)+\
    # np.array(((df0['SPREAD_MODEL_i']+df0['SPREADERR_MODEL_i'])>0.003).values, dtype=int)+\
    # np.array(((df0['SPREAD_MODEL_i']-df0['SPREADERR_MODEL_i'])>0.003).values, dtype=int)
    # dff=df0[df0['EXTENDED_COADD']>1]
    # dff_star=df0[df0['EXTENDED_COADD']<2]
    df1=df0[df0['FLAGS_i']<4].copy()
    dff=df1[(df1['SPREAD_MODEL_i2'])>0.0035]

    # dff_star=df0[(df0['MAG_AUTO_i']<-8) & (df0['SPREAD_MODEL_i']<0.10)] #+5/3.*df0['SPREADERR_MODEL_i'] <0.002
    dff_star=df0[(df0['SPREAD_MODEL_i2']<0.004)&(df0['MAG_AUTO_i']<-16)&(df0['MAG_AUTO_i']>-18.5)]
    # dff_star=df0[np.abs(df0['SPREAD_MODEL_i2'])<0.003]
    # dff_star=df0[df0['CLASS_STAR_i']>0.9]

    dff.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/galaxy_psf_total_%s.csv"%field)
    dff_star.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/star_psf_total_%s.csv"%field)


def pisco_photometry_psf_v4(field, mode='psf', mode2mass='', slr=True):  #mode2mass: '' vs '_no2mass'
    def slr_running_psf(field, infile="None", mode="psf", mode2mass='', bigmacs="pisco_pipeline/big-macs-calibrate-master"):
        """
        slr_running: running SLR script from github.com/patkel/big-macs-calibrate to get a calibrated magnitude
        INPUT:
        - field: object of interset e.g., 'Field026'
        - bigmacs: the location for "big-macs-calibrate" directoty
        OUTPUT:
        - a new table with added columns with name MAG_g,...,MAGERR_g,...
        """
        slrdir = 'slr_output'
        pyfile = os.path.join(bigmacs, 'fit_locus.py')
        # cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 15 -l -r ALPHA_J2000_i -d DELTA_J2000_i -j --plot=PLOTS_%s_%s" \
        #     % (pyfile, infile, os.path.join(bigmacs, "coadd_mag_sex_%s%s.columns"%(mode,'')), mode, field) 
        if mode2mass=='':
            cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 15 -l -r ALPHA_J2000_i -d DELTA_J2000_i -j --plot=PLOTS_%s_%s" \
                % (pyfile, infile, os.path.join(bigmacs, "coadd_mag_sex_%s%s.columns"%(mode,mode2mass)), mode, field)  #'' vs '_no2mass'
        elif mode2mass=='_no2mass':
            cmd = "python %s --file %s --columns %s --extension 1 --bootstrap 15 -l -r ALPHA_J2000_i -d DELTA_J2000_i --plot=PLOTS_%s_%s" \
                    % (pyfile, infile, os.path.join(bigmacs, "coadd_mag_sex_%s%s.columns"%(mode,mode2mass)), mode, field)  #'' vs '_no2mass'
        print cmd
        subprocess.check_call(shlex.split(cmd))

    def update_color(fname, table, mode='psf'):
        """
        update_color: using the output from SLR, update to the correct magnitude
        INPUT:
        - fname: input file from SLR output (...offsets.list)
        - table: the table that we want to update the value (from column magg,etc to MAG_g,etc)
        OUTPUT:
        - a new table with added columns with name MAG_g,...,MAGERR_g,...
        """
        print fname
        with open(fname) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        # print content
        # if len(content)==8:
        #     red_content=content[4:]
        # elif len(content)==10:    
        #     red_content=content[5:-1]
        if len(content)==7:
            red_content=content[4:]
        elif len(content)==9:    
            red_content=content[5:-1]
        band = [x.split(' ')[0][-1] for x in red_content]
        corr = [float(x.split(' ')[1]) for x in red_content]
        ecorr = [float(x.split(' ')[3]) for x in red_content]    
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

        table['MAG_' + band[0]] = table['MAG_%s_'%MODE1 + band[0]] + corr[0]
        table['MAG_' + band[1]] = table['MAG_%s_'%MODE1 + band[1]] + corr[1]
        table['MAG_' + band[2]] = table['MAG_%s_'%MODE1 + band[2]] + corr[2]
        # table['MAG_' + band[3]] = table['MAG_%s_'%MODE1 + band[3]] + corr[3]
        table['MAGERR_' + band[0]] = (table['MAGERR_%s_'%MODE1 + band[0]]**2)**0.5# + ecorr[0]**2)**0.5
        table['MAGERR_' + band[1]] = (table['MAGERR_%s_'%MODE1 + band[1]]**2)**0.5# + ecorr[1]**2)**0.5
        table['MAGERR_' + band[2]] = (table['MAGERR_%s_'%MODE1 + band[2]]**2)**0.5# + ecorr[2]**2)**0.5
        # table['MAGERR_' + band[3]] = (table['MAGERR_%s_'%MODE1 + band[3]]**2)# + ecorr[3]**2)**0.5
        # table['MAGERR_' + band[0]] = (table['MAGERR_%s_'%MODE1 + band[0]]**2)**0.5
        # table['MAGERR_' + band[1]] = (table['MAGERR_%s_'%MODE1 + band[1]]**2)**0.5
        # table['MAGERR_' + band[2]] = (table['MAGERR_%s_'%MODE1 + band[2]]**2)**0.5
        # table['MAGERR_' + band[3]] = (table['MAGERR_%s_'%MODE1 + band[3]]**2)**0.5
        return table

    slrdir = 'slr_output'
    total3 = Table.from_pandas(pd.read_csv(
        "/Users/taweewat/Documents/pisco_code/slr_output/star_psf_total_%s.csv" % field))  # star_psf_total_gaia
    total3=total3[['NUMBER','ALPHA_J2000_i','DELTA_J2000_i','XWIN_IMAGE_i','YWIN_IMAGE_i',\
            'MAG_APER_i','MAGERR_APER_i','MAG_APER_g','MAGERR_APER_g','MAG_APER_r',\
            'MAGERR_APER_r','MAG_AUTO_i','MAGERR_AUTO_i',\
            'MAG_AUTO_g','MAGERR_AUTO_g','MAG_AUTO_r','MAGERR_AUTO_r','MAG_SPHEROID_i','MAGERR_SPHEROID_i','MAG_SPHEROID_g',\
            'MAGERR_SPHEROID_g','MAG_SPHEROID_r','MAGERR_SPHEROID_r','CLASS_STAR_i','CLASS_STAR_g','CLASS_STAR_r',\
            'FLAGS_g','FLAGS_r','FLAGS_i','MAG_PSF_g',\
            'MAG_PSF_r','MAG_PSF_i','MAGERR_PSF_g','MAGERR_PSF_r',\
            'MAGERR_PSF_i','MAG_MODEL_g','MAG_MODEL_r',\
            'MAG_MODEL_i','MAGERR_MODEL_g','MAGERR_MODEL_r',\
            'MAGERR_MODEL_i','SPREAD_MODEL_g','SPREAD_MODEL_r',\
            'SPREAD_MODEL_i','SPREADERR_MODEL_g','SPREADERR_MODEL_r',\
            'SPREADERR_MODEL_i']]

    print 'number of stars =', len(total3)

    if (mode2mass==''):# and (mode=='psf'):
        starpsfmode = '_psf'
    # elif (mode2mass=='') and (mode=='model'):
    #     starpsfmode = '_psf' #'_model'
    elif (mode2mass=='_no2mass'):# and (mode=='model'):
        starpsfmode ='_no2mass' #'_model_no2mass'
    # elif (mode2mass == '_no2mass') and (mode=='psf'):
    #     starpsfmode = '_no2mass'

    # total3.write(slrdir+'/star_psf%s_%s_%i.fits' % ('_psf',field,0), overwrite=True) #with 2MASS stars: star_psf_psf_%s_%i.fits
    total3.write(slrdir + '/star_psf%s_%s_%i.fits' % (starpsfmode, field, 0),
                 overwrite=True)  # no 2MASS star mode vs , '_psf' vs '_no2mass'

    if slr:
        slr_running_psf(field, infile=slrdir + '/star_psf%s_%s_%i.fits' %
                        (starpsfmode, field, 0), mode='psf', mode2mass=mode2mass)  # '_psf' vs '_no2mass'

    total_gal=Table.from_pandas(pd.read_csv("/Users/taweewat/Documents/pisco_code/slr_output/galaxy_psf_total_%s.csv"%(field)))
    print 'mode=', mode, '/star_psf%s_%s_%i.fits.offsets.list' % (starpsfmode, field, 0)
    ntotal_gal = update_color(slrdir + '/star_psf%s_%s_%i.fits.offsets.list' %
                              (starpsfmode, field, 0), total_gal, mode=mode)  # '' vs '_no2mass
    ntotal_gal.write(os.path.join(
        slrdir, 'galaxy_%s%s_ntotal_%s.csv' % (mode, mode2mass, field)), overwrite=True)  # '' vs '_no2mass'

def make_images(field,ax=None):
    dir='/Users/taweewat/Documents/pisco_code/Chips_images/'
    try:
        ax.imshow(image.imread(dir+"aplpy_panstar_%s_img4.jpeg"%field))
    except:
        ax.imshow(image.imread(dir+"aplpy_panstar_%s_img.jpeg"%field))
    # ax.imshow(image.imread(dir+"aplpy4_%s_img4.jpeg"%field))
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axis('off')
    return None

# def sur_pro(r): #Mpc
#     def fn(x):
#         if x>=1:
#             return 1.-(2/np.sqrt(x**2-1)*np.arctan(np.sqrt((x-1.)/(x+1.))))
#         elif x<1:
#             return 1.-(2/np.sqrt(1-x**2)*np.arctanh(np.sqrt((1.-x)/(x+1.))))
#     rs=0.15/0.71 #Mpc
#     if r>=(0.1/0.71):
#         return 1/((r/rs)**2-1)*fn(r/rs)
#     elif r<(0.1/0.71):
#         return 1./(((0.1/0.71)/rs)**2-1)*fn((0.1/0.71)/rs)
# def k_NFW():
#     def integrated(y):
#         return 1./integrate.quad(lambda r: 2*np.pi*r*sur_pro(r),0,y)[0]
#     xy=np.logspace(-3,3,num=30)
#     X = np.log(xy)
#     Y = np.log([integrated(np.e**(y)) for y in X])
#     Z=np.polyfit(X,Y,6)
#     k_NFW = np.poly1d(Z)
#     return k_NFW
# def sur_pro_prob(r,rc,k_NFW): #(Mpc,Mpc) # Weighted based on the distance from the center (Rykoff+12)
#     return np.e**(k_NFW(np.log(rc)))*sur_pro(r)


def sur_pro(r,rc): #(arcmin)
    def fn(x):
        if x>=1:
            return 1.-(2/np.sqrt(x**2-1)*np.arctan(np.sqrt((x-1.)/(x+1.))))
        elif x<1:
            return 1.-(2/np.sqrt(1-x**2)*np.arctanh(np.sqrt((1.-x)/(x+1.))))
    rs=0.15/0.71 #Mpc
    if r>=rc:
        return 1/((r/rs)**2-1)*fn(r/rs)
    elif r<rc:
        return 1./((rc/rs)**2-1)*fn(rc/rs)
def sur_pro_prob_ang(r,rc):
    return sur_pro(r,rc)/sur_pro(0.2,rc)

# def sur_pro_prob_ang(r,r_c):
#     if r < r_c:
#         return 1.
#     else: 
#         return 0. 
# 'ezmodel2_bc03_zf3.0_chab_0.02_exp_0.1.txt', 'ezmodel2_c09_zf3.0_chab_0.02_exp_0.1.txt'
name=['z','dist','age','mass','Abs_g','App_g','kcorr_g','Abs_r',\
      'App_r','kcorr_r','Abs_i','App_i','kcorr_i','Abs_z','App_z','kcorr_z']
df=pd.read_csv('/Users/taweewat/Documents/red_sequence/rsz/model/'+\
# 'ezmodel2_bc03_zf2.5_chab_0.016_exp_0.1.txt',
'ezmodel2_bc03_zf2.5_chab_0.02_exp_0.1.txt',
            #    'ezmodel2_c09_zf3.0_chab_0.02_exp_0.1.txt',
               skiprows=27,delim_whitespace=True,names=name)
df=df[(df.z>=0.1) & (df.z<1.)]
z_new=np.arange(0.1, 0.95, 0.0025)
Appi_new = interpolate.splev(z_new, interpolate.splrep(df.z, df.App_i, s=0), der=0)
Appi_f = interpolate.interp1d(df.z, df.App_i, kind='cubic')

#all extra options
extra_name= 'gremove_silk_zf3_c09_noebv_model_complete_gaia' #'gremove_lum_silk_zf2.5_c09_11', 'gremove_silk_zf3_c09_noebv_model_complete_no2mass'
gremove = False # remove non-detect g objects from the list
duplicate = False # remove duplicate redshift (uncertain)
colorerr = True  # add redshift with color_error taken into account
transparent = True # make transparent plot for flip book
img_filp = False  # make image flip from transparent
img_redshift = True # make image with redshift for each object

def linear_rmi(x0,redshift):
    x=df.z[:-11] #-12
    y=(df.App_r-df.App_i)[:-11] #-12
    yhat = np.polyfit(x, y, 5) #5 vs 9
    f_rmi = np.poly1d(yhat)
    slope=-0.0222174237562*1.007
    # Appi0=Appi_new[np.where(abs(z_new-redshift)<=1e-9)[0][0]]
    Appi0=Appi_f(redshift)
    return slope*(x0-Appi0)+f_rmi(redshift)

def linear_gmr(x0,redshift):
    x=df.z[:-24] #-25
    y=(df.App_g-df.App_r)[:-24] #-25
    yhat = np.polyfit(x, y, 5)
    f_gmr = np.poly1d(yhat)
    slope=-0.0133824600874*1.646
    # Appi0=Appi_new[np.where(abs(z_new-redshift)<=1e-9)[0][0]]
    Appi0=Appi_f(redshift)
    return slope*(x0-Appi0)+f_gmr(redshift)

def linear_gmi(x0,redshift):
    x=df.z[:-9] 
    y=(df.App_g-df.App_i)[:-9] 
    yhat = np.polyfit(x, y, 5)
    f_gmi = np.poly1d(yhat) 
    Appi0=Appi_f(redshift)
    slope = -0.04589707934164738 * 1.481
    return slope*(x0-Appi0)+f_gmi(redshift)

def find_fits_dir(field):
    home = '/Users/taweewat/Documents/pisco_code/'
    dirs = ['ut170103/', 'ut170104/', 'ut170619/', 'ut170621/',\
            'ut170624/', 'ut171208/', 'ut171209/', 'ut171212/']
    myReg = re.compile(r'(%s_A).*' % field)
    for di in dirs:
        diri = home + di
        for text in os.listdir(diri):
            if myReg.search(text) != None:
                # filename = myReg.search(text).group()
                allfilename = diri
    allfilename='ut170103/'
    return allfilename

dir_dict = dict(zip(['ut170103/','ut170104/','ut170619/',\
'ut170621/','ut170624/','ut171208/','ut171209/','ut171212/'], np.arange(1, 9)))

def find_ra_dec(field):
    if field == 'PKS1353':
        RA = 209.0225
        DEC = -34.3530556
        redshift = 0.223
    elif field == 'CHIPS2249-2808':
        RA = 336.99975202151825
        DEC = -43.57623068466675
        redshift = -1
    elif field == 'CHIPS2246-2854':
        RA = 335.7855174238757
        DEC = -34.934569299688185
        redshift = -1
    elif field[0:5] == 'Field':
        base = pd.read_csv(
            '/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs.csv')
        RA = base[base.name == field].ra.values[0]
        DEC = base[base.name == field].dec.values[0]
        redshift = base[base.name == field].redshift.values[0]
    elif field[0:5] == 'CHIPS':
        base = pd.read_csv(
            '/Users/taweewat/Documents/red_sequence/chips_all_obj.csv', index_col=0)
        RA = base[base.chips == field].ra.values[0]
        DEC = base[base.chips == field].dec.values[0]
        redshift = base[base.chips == field].redshift.values[0]
    elif field[0:4] == 'SDSS':
        base = pd.read_csv(
            '/Users/taweewat/Documents/xray_project/ned-result/final_sdss_cut5.csv', index_col=0)
        RA = base[base.name == field].RA.values[0]
        DEC = base[base.name == field].DEC.values[0]
        redshift = base[base.name == field].redshift.values[0]
    return RA, DEC, redshift

def pisco_tilt_resequence(field, mode='psf', mode2mass=''):
    RA, DEC, redshift = find_ra_dec(field)

    if redshift!=-1:
        qso_redshift=redshift
    else:
        qso_redshift=0.2

    print 'RA', RA
    print 'DEC', DEC

    ebv = ebvpy.calc_ebv(ra=[RA],dec=[DEC]); print 'ebv:', ebv[0]
    # ebv_g=ebvpy.calc_color_correction('g', ebv)[0]
    # ebv_r=ebvpy.calc_color_correction('r', ebv)[0]
    # ebv_i=ebvpy.calc_color_correction('i', ebv)[0]
    # ebv_z=0.0
    ebv_g,ebv_r,ebv_i,ebv_z=0.0,0.0,0.0,0.0  #no longer use reddening correction because it is already included in SLR
    print 'ebv_g:', ebv_g, 'ebv_r:', ebv_r, 'ebv_i:', ebv_i

    param_izp=read_param_izp(mode) #i zero point

    # fname = "/Users/taweewat/Documents/pisco_code/slr_output/galaxy_ntotal_%s.csv"%field
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/galaxy_%s%s_ntotal_%s.csv" % (
        mode, mode2mass, field)  # '' vs '_no2mass'
    df0 = pd.read_csv(fname,index_col=0)
    # gremove=True
    if gremove:
        nog=len(df0[df0['MAG_PSF_g'] >= 50.]); print "no g detected:", nog
        df0 = df0[df0['MAG_PSF_g'] < 50.].copy()  # cut out not detected objects in g band
    else: 
        nog=0

    c5 = SkyCoord(ra=df0['ALPHA_J2000_i'].values*u.degree, dec=df0['DELTA_J2000_i'].values*u.degree)
    c0 = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    sep = c5.separation(c0)
    df0['sep(deg)']=sep
    df0['sep(Mpc)']=sep*60.*cosmo.kpc_proper_per_arcmin(qso_redshift).value/1e3
    cut=df0

    dfi = cut#.drop_duplicates(subset=['XWIN_WORLD', 'YWIN_WORLD'], keep='first').copy()
    print 'duplicates:', len(df0), len(dfi)
    
    # Added Galactic Reddening (6/16/18)
    if mode2mass == '':
        dfi['MAG_i']=dfi['MAG_i']-ebv_i
        dfi['MAG_g']=dfi['MAG_g']-ebv_g
        dfi['MAG_r']=dfi['MAG_r']-ebv_r
    # Use i Zero Point from each day and g,r zero point fron the color (6/22/18)
    elif mode2mass == '_no2mass':
        dfi['MAG_i']=dfi['MAG_i']-ebv_i+param_izp['i_zp_day9']#%dir_dict[find_fits_dir(field)[-9:]]]
        dfi['MAG_g']=dfi['MAG_g']-ebv_g+param_izp['i_zp_day9']#%dir_dict[find_fits_dir(field)[-9:]]]
        dfi['MAG_r']=dfi['MAG_r']-ebv_r+param_izp['i_zp_day9']#%dir_dict[find_fits_dir(field)[-9:]]]
        # dfi['MAG_z']=dfi['MAG_z']-ebv_z+param_izp['i_zp_day%i'%dir_dict[find_fits_dir(field)[-9:]]]
        # dfi['MAGERR_i']=np.sqrt(dfi['MAGERR_i']**2-(99**2))


    # dfi=dfi[dfi['MAG_i']<21.5].copy()

    # dfi=dfi[dfi.MAGERR_g<0.5]
    # dfi=dfi[(dfi.MAG_g<100)&(dfi.MAG_i<100)&(dfi.MAG_r<100)]
    # dfi=dfi[(dfi.FLAGS_g<5)&(dfi.FLAGS_r<5)&(dfi.FLAGS_i<5)&(dfi.FLAGS_z<5)]
    print field, qso_redshift, df0.shape, cut.shape, dfi.shape, dfi['sep(deg)'].max(), dfi['sep(Mpc)'].max()

    norm = matplotlib.colors.Normalize(vmin=0.15,vmax=0.675)
    c_m = matplotlib.cm.cool
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    I=np.arange(16,24,0.01)
    dfi.loc[:,"z_gmr"] = np.nan
    dfi.loc[:,"z_rmi"] = np.nan
    dfi.loc[:,"w_gmr"] = np.nan
    dfi.loc[:,"w_rmi"] = np.nan
    dfi.loc[:,"w_col_gmr"] = np.nan
    dfi.loc[:,"w_col_rmi"] = np.nan
    # dfi.loc[:,"z_gmi"] = np.nan
    # dfi.loc[:,"w_gmi"] = np.nan
    # dfi.loc[:,"w_col_gmi"] = np.nan

    # k_NFW0=k_NFW()

    bin_width=0.035 #0.025

    bins_gmr_cen = np.arange(0.12315, 0.33315+0.01, bin_width)
    bins_gmr_edge = np.arange(0.10565, 0.35065 + 0.01, bin_width)
    # bins_gmr_cen = np.arange(0.15815, 0.33315+0.01, bin_width) # bins_gmr_cen = np.arange(0.15, 0.325+0.01, bin_width)
    # bins_gmr_edge = np.arange(0.14065, 0.35065+0.01, bin_width) # bins_gmr_edge = np.arange(0.1325, 0.3425+0.01, bin_width)
    
    bins_rmi_cen = np.arange(0.36815, 0.64815+0.01, bin_width) # bins_rmi_cen = np.arange(0.36, 0.675+0.01, bin_width)
    bins_rmi_edge = np.arange(0.35065, 0.66565+0.01, bin_width) # bins_rmi_edge = np.arange(0.3425, 0.6925+0.01, bin_width)

    z_rmi,w_rmi,w_col_rmi=[],[],[]
    for i, row in dfi.iterrows():
        for z in bins_rmi_cen:
            # if row['MAG_i'] < -18+5.*np.log10(ex.d_L(z)*1e6)-5.:
             # if row['MAG_i'] < magi_cut_rmi:
            # if np.sqrt(row['MAGERR_r']**2+row['MAGERR_i']**2)<0.134: #np.mean(f_rmi(x+0.07)-f_rmi(x))
            if np.sqrt(row['MAGERR_r']**2+row['MAGERR_i']**2)<0.067*4.5:#1.5: #0.067*1.5
            # if np.sqrt(row['MAGERR_r']**2+row['MAGERR_i']**2)<1:
                rmi=row['MAG_r']-row['MAG_i']
                # rmierr=np.sqrt(row['MAGERR_r']**2+row['MAGERR_i']**2)
                low_edge=linear_rmi(row['MAG_i'],round(z-0.0175,4)) #0.0125
                high_edge=linear_rmi(row['MAG_i'],round(z+0.0175,4)) #0.0125
                if (rmi > low_edge) & (rmi <= high_edge):
                    # if (np.sqrt(row['MAGERR_r']**2+row['MAGERR_i']**2) < 3.5*(high_edge-low_edge)):
                    z_rmi.append(round(z,3))
                    # wrmi0=sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0)
                    wrmi0=sur_pro_prob_ang(row['sep(deg)']*60, 1); w_rmi.append(wrmi0) #arcmin
                    # w_col_rmi0=scipy.stats.norm(rmi,rmierr).cdf(high_edge)-scipy.stats.norm(rmi,rmierr).cdf(low_edge); w_col_rmi.append(w_col_rmi0)
                    w_col_rmi0=1.; w_col_rmi.append(w_col_rmi0)
                    dfi.loc[i,"z_rmi"]=z
                    dfi.loc[i,"w_rmi"]=wrmi0
                    dfi.loc[i,"w_col_rmi"]=w_col_rmi0

    z_gmr,w_gmr,w_col_gmr=[],[],[]
    for i, row in dfi.iterrows():
        for z in bins_gmr_cen:
            # if row['MAG_i'] < -18+5.*np.log10(ex.d_L(z)*1e6)-5.:
            # if row['MAG_i'] < magi_cut_gmr:
            # if np.sqrt(row['MAGERR_g']**2+row['MAGERR_r']**2)<0.165: #np.mean(f_gmr(x+0.07)-f_gmr(x))
            if np.sqrt(row['MAGERR_g']**2+row['MAGERR_r']**2)<0.0825*4.5:#1.5: #0.0825*1.5
            # if np.sqrt(row['MAGERR_g']**2+row['MAGERR_r']**2)<1:
                gmr=row['MAG_g']-row['MAG_r']
                # gmrerr=np.sqrt((row['MAGERR_g'])**2+row['MAGERR_r']**2) #add factor 2.2 to reduce the g error to be similar to other bands
                low_edge=linear_gmr(row['MAG_i'],round(z-0.0175,4)) #0.0125
                high_edge=linear_gmr(row['MAG_i'],round(z+0.0175,4)) #0.0125
                if (gmr > low_edge) & (gmr <= high_edge):
                    # if (np.sqrt(row['MAGERR_g']**2+row['MAGERR_r']**2) < 3.5*(high_edge-low_edge)):
                    z_gmr.append(round(z,3))
                    # w_col_gmr0=scipy.stats.norm(gmr,gmrerr).cdf(high_edge)-scipy.stats.norm(gmr,gmrerr).cdf(low_edge); w_col_gmr.append(w_col_gmr0)
                    w_col_gmr0=1.; w_col_gmr.append(w_col_gmr0)
                    # wgmr0=sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0); w_gmr.append(wgmr0)
                    wgmr0 = sur_pro_prob_ang(row['sep(deg)'] * 60, 1); w_gmr.append(wgmr0)  # arcmin
                    dfi.loc[i,"z_gmr"]=z
                    dfi.loc[i,"w_gmr"]=wgmr0
                    dfi.loc[i,"w_col_gmr"]=w_col_gmr0

    # z_gmi,w_gmi,w_col_gmi=[],[],[]
    # for i, row in dfi.iterrows():
    #     # for z in np.arange(0.15,0.35,0.025):
    #     for z in np.arange(0.15,0.7,0.035):
    #         if row['MAG_i'] < -18+5.*np.log10(ex.d_L(z)*1e6)-5.:
    #             gmi=row['MAG_g']-row['MAG_i']
    #             gmierr=np.sqrt((row['MAGERR_g']/2.2)**2+row['MAGERR_i']**2) #add factor 2.2 to reduce the g error to be similar to other bands
    #             low_edge=linear_gmi(row['MAG_i'],round(z-0.0175,4)) #0.0125
    #             high_edge=linear_gmi(row['MAG_i'],round(z+0.0175,4)) #0.0125
    #             if (gmi > low_edge) & (gmi <= high_edge):
    #                 # if (np.sqrt(row['MAGERR_g']**2+row['MAGERR_r']**2) < 3.5*(high_edge-low_edge)):
    #                 z_gmi.append(round(z,3))
    #                 # w_col_gmi0=scipy.stats.norm(gmi,gmierr).cdf(high_edge)-scipy.stats.norm(gmi,gmierr).cdf(low_edge); w_col_gmi.append(w_col_gmi0)
    #                 w_col_gmi0=1.; w_col_gmi.append(w_col_gmi0)
    #                 # wgmi0=sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0); w_gmi.append(wgmi0)
    #                 wgmi0 = sur_pro_prob_ang(row['sep(deg)'] * 60, 1.); w_gmi.append(wgmi0)  # arcmin
    #                 dfi.loc[i,"z_gmi"]=z
    #                 dfi.loc[i,"w_gmi"]=wgmi0
    #                 dfi.loc[i,"w_col_gmi"]=w_col_gmi0

    # ns1,xs1=np.histogram(z_gmr,bins=np.arange(0.125,0.35,0.025),weights=w_gmr)
    ns1,xs1=np.histogram(z_gmr,bins=bins_gmr_edge,weights=np.array(w_gmr)*np.array(w_col_gmr)) #0.15-0.325
    bin_cen1 = (xs1[:-1] + xs1[1:])/2
    # ns2,xs2=np.histogram(z_rmi,bins=np.arange(0.325,0.7,0.025),weights=w_rmi)
    ns2,xs2=np.histogram(z_rmi,bins=bins_rmi_edge,weights=np.array(w_rmi)*np.array(w_col_rmi)) #0.36-0.675
    bin_cen2 = (xs2[:-1] + xs2[1:])/2
    # z_total=np.append(xs1[:-1],xs2[:-1])
    z_total=np.append(bin_cen1, bin_cen2)
    n_total=np.append(ns1,ns2)
    z_max=z_total[np.where(n_total==np.max(n_total))[0][0]]
    n_median = np.median(n_total[n_total != 0])
    n_mean = np.mean(n_total)

    n_bkg = np.mean(sorted(n_total)[2:-2]);
    z_total_added = np.insert(
        np.append(z_total, z_total[-1] + bin_width), 0, z_total[0] - bin_width)
    n_total_added = np.insert(np.append(n_total, 0), 0, 0) - n_bkg
    # print 'n_total_added', n_total_added
    
    lumfn=pd.read_csv('/Users/taweewat/Documents/red_sequence/coma_cluster_luminosity_function/schecter_fn.csv',\
    names=['M_r','theta(M)Mpc^-3']) 
    h=0.7
    x=lumfn['M_r']+5*np.log10(h);
    y=lumfn['theta(M)Mpc^-3']*(h**3)
    f1d=interp1d(x, y,kind='cubic')
    def lum_function(M):
        alpha = -1.20
        Nb = np.log(10) / 2.5 * 0.002 * (70 / 50.)**3
        Mb_s = -21. + 5 * np.log10(70 / 50.)
        return Nb * (10.**(0.4 * (alpha + 1) * (Mb_s - M))) * np.exp(-10.**(0.4 * (Mb_s - M)))
    def distance(z):
        return cosmo.luminosity_distance(z).value
    def abs_mag(m, z):
        return m - 5 * np.log10(distance(z) * 1e6) + 5
    def NFW_profile(r):
        rs = 1.  # Mpc
        rho0 = 500.
        return rho0 / (r / rs * (1 + r / rs)**2)
    lum_fn = lambda z: integrate.quad( f1d, -23.455, abs_mag(22.25, z))[0]
    lum_vfn = np.vectorize(lum_fn)
    dense_fn = lambda z: integrate.quad(NFW_profile,0.001,cosmo.kpc_proper_per_arcmin(z).value/1e3)[0]
    dense_vfn = np.vectorize(dense_fn)
    n_total_adj=n_total_added #/(lum_vfn(z_total_added)*dense_vfn(z_total_added)) (adjusted the peak before picking it)
    print 'n_total_added:', n_total_added
    print 'n_total_adj:', n_total_adj

    indi = np.where(n_total_adj == np.max(n_total_adj))[0][0]
    # indi = np.where(n_total_added == np.max(n_total_added))[0][0]
    z_fit = z_total_added[[indi - 1, indi, indi + 1]]; print 'z_fit', z_fit
    n_fit = n_total_added[[indi - 1, indi, indi + 1]]; print 'n_fit', n_fit
    def gaussian_func(x, a, mu):
        sigma=0.025
        return a * np.exp(-(x-mu)**2/(2*(sigma**2)))
    
    if (n_fit[0]<0.) and (n_fit[2]<0.):
        popt, pcov = curve_fit(gaussian_func, z_fit, [0,n_fit[1],0], p0=[n_fit[1],z_fit[1]])
    else:
        popt, pcov = curve_fit(gaussian_func, z_fit,
                               n_fit, p0=[n_fit[1], z_fit[1]])
    # signal=tuple(popt)[0]

    # def v_func(z):
    #     return (z**2+2*z)/(z**2+2*z+2)
    # signal=((np.max(n_total)-np.mean(n_total))*(v_func(z_max)*(4000))**2)/5.3e6 #normalization for r~1 at z~0.15
    # signal = (
    #     (tuple(popt)[0]) * (cosmo.luminosity_distance(tuple(popt)[1]).value)**1.5) / 5.3e5  # normalization for r~1 at z~0.15

    def lum_function(M):
        alpha = -1.20
        Nb = np.log(10) / 2.5 * 0.002 * (70 / 50.)**3
        Mb_s = -21. + 5 * np.log10(70 / 50.)
        return Nb * (10.**(0.4 * (alpha + 1) * (Mb_s - M))) * np.exp(-10.**(0.4 * (Mb_s - M)))
    def distance(z):
        return cosmo.luminosity_distance(z).value
    def abs_mag(m, z):
        return m - 5 * np.log10(distance(z) * 1e6) + 5
    def NFW_profile(r):
        rs = 1.  # Mpc
        rho0 = 500.
        return rho0 / (r / rs * (1 + r / rs)**2)

    lumfn=pd.read_csv('/Users/taweewat/Documents/red_sequence/coma_cluster_luminosity_function/schecter_fn.csv',\
    names=['M_r','theta(M)Mpc^-3']) 
    h=0.7
    x=lumfn['M_r']+5*np.log10(h);
    y=lumfn['theta(M)Mpc^-3']*(h**3)
    f1d=interp1d(x, y,kind='cubic')

    # lum_factor = integrate.quad(lum_function, -24, abs_mag(21.60, tuple(popt)[1]))[0]
    # lum_factor = cosmo.luminosity_distance(tuple(popt)[1]).value**-1.5*100
    lum_factor = integrate.quad( f1d, -23.455, abs_mag(22.25, tuple(popt)[1]))[0]  
    #-23.455: min abs Mag from schecter_fn.csv, 22.25: median of Mag r
    density_factor=integrate.quad(NFW_profile, 0.001, cosmo.kpc_proper_per_arcmin(tuple(popt)[1]).value/1e3)[0]
    signal = tuple(popt)[0] / (lum_factor * density_factor)
    z_max_fit = tuple(popt)[1]

    print 'z_max_fit', z_max_fit
    print 'lum_factor:', lum_factor
    print 'density_factor', density_factor

    # duplicate=False ## set duplication
    if duplicate:
        dff = dfi.copy()
        fig=plt.figure(figsize=(10,3.5))
        plt.subplot(1,3,1)
        bin_width=0.035
        ns1, xs1 = np.histogram(dff['z_gmr'].dropna(), bins=bins_gmr_edge, weights=np.array(dff['w_gmr'].dropna()))
        bin_cen1 = (xs1[:-1] + xs1[1:])/2
        ns2,xs2=np.histogram(dff['z_rmi'].dropna(),bins=bins_rmi_edge,weights=np.array(dff['w_rmi'].dropna()))
        bin_cen2 = (xs2[:-1] + xs2[1:])/2
        plt.bar(bin_cen1, ns1, width=0.035, color='#ff7f0e')
        plt.bar(bin_cen2, ns2, width=0.035, color='#1f77b4')
        if np.max(np.append(ns1,ns2))<30:
            plt.ylim(0,30)
        plt.xlabel('all objects')
        plt.xlim(0.1, 0.7)
        plt.subplot(1,3,2)
        dff_dup=dff.dropna()
        ns_dup1,xs_dup1=np.histogram(dff_dup['z_gmr'],bins=bins_gmr_edge,weights=np.array(dff_dup['w_gmr']))
        bin_cen1 = (xs1[:-1] + xs1[1:])/2
        ns_dup2,xs_dup2=np.histogram(dff_dup['z_rmi'],bins=bins_rmi_edge,weights=np.array(dff_dup['w_rmi']))
        bin_cen2 = (xs2[:-1] + xs2[1:])/2
        plt.bar(bin_cen1, ns_dup1, width=0.035, color='#ff7f0e')
        plt.bar(bin_cen2, ns_dup2, width=0.035, color='#1f77b4')
        if np.max(np.append(ns_dup1,ns_dup2))<30:
            plt.ylim(0,30)
        plt.xlabel('duplicate')
        plt.xlim(0.1,0.7)
        
        n_total_dup=np.append(ns1-ns_dup1,ns2-ns_dup2)
        n_bkg_dup = np.mean(sorted(n_total_dup)[2:-2])
        z_total_added = np.insert(
            np.append(z_total, z_total[-1] + bin_width), 0, z_total[0] - bin_width)
        n_total_added = np.insert(np.append(n_total_dup, 0), 0, 0) - n_bkg_dup
        indi = np.where(n_total_added == np.max(n_total_added))[0][0]
        z_fit_dup = z_total_added[[indi - 1, indi, indi + 1]]
        n_fit_dup = n_total_added[[indi - 1, indi, indi + 1]]
        popt_dup, pcov_dup = curve_fit(gaussian_func, z_fit_dup, n_fit_dup, p0=[n_fit_dup[1],z_fit_dup[1]])
        signal_dup= tuple(popt_dup)[0] / (lum_factor * density_factor)
        z_max_fit=tuple(popt_dup)[1]

        plt.subplot(1, 3, 3)
        plt.bar(bin_cen1, ns1 - ns_dup1, width=0.035, color='#ff7f0e')
        plt.bar(bin_cen2, ns2 - ns_dup2, width=0.035, color='#1f77b4')
        if np.max(n_total_dup)<30:
            plt.ylim(0,30)
        plt.xlabel('all objects-duplicate')
        plt.xlim(0.1, 0.7)
        plt.axvline(z_max_fit,ls='--',color='purple',label='z_max:%.2f'%z_max_fit)
        plt.axvline(redshift,color='red',label='z:%.2f'%redshift)
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_dup_%s_all_%.3f_%s_tilted.png' %
                    (mode, signal_dup, field), dpi=120)
        plt.close(fig)
    else:
        n_total_dup=0

    ## Plot the figure
    cmap=matplotlib.cm.RdYlGn
    if duplicate or colorerr:
        fig,ax=plt.subplots(1,5,figsize=(25,5))
    else: 
        fig,ax=plt.subplots(1,4,figsize=(20,5))
    make_images(field,ax[0])

    norm = matplotlib.colors.Normalize(vmin=0.01,vmax=2)
    dfi_ri=dfi.loc[dfi['z_rmi'].dropna().index]
    ax[1].scatter(dfi['MAG_i'],dfi['MAG_r']-dfi['MAG_i'],c='black',alpha=0.1)#dfi['w_rmi'],cmap=cmap)
    ax[1].scatter(dfi_ri['MAG_i'],dfi_ri['MAG_r']-dfi_ri['MAG_i'],c=dfi_ri['w_rmi'],cmap=cmap)#,norm=norm)
    ax[1].errorbar(dfi_ri['MAG_i'],dfi_ri['MAG_r']-dfi_ri['MAG_i'],xerr=dfi_ri['MAGERR_i'],yerr=np.sqrt(dfi_ri['MAGERR_r']**2+dfi_ri['MAGERR_i']**2),fmt='none',c='k',alpha=0.05)
    # plt.plot(df.App_i,df.App_r-df.App_i,'.')
    # ax[1].axhline(xs[:-1][(xs[:-1]<1.33) & (xs[:-1]>0.6)][0],lw=0.7,color='green')
    for z in bins_rmi_cen:
        ax[1].plot(I,linear_rmi(I,round(z,4)),color=s_m.to_rgba(z))
    ax[1].set_ylim(0.25,1.5)
    ax[1].set_xlim(16,24)
    # cbar=plt.colorbar(s_m)
    ax[1].set_xlabel('I')
    ax[1].set_ylabel('R-I')
    ax[1].set_title('z=0.35-0.675')#, icut:'+str(magi_cut_rmi))

    # plt.plot([corr_f(z) for z in df.z.values[5:-12]],df.App_r[5:-12]-df.App_i[5:-12],'-')
    dfi_gr=dfi.loc[dfi['z_gmr'].dropna().index]
    ax[2].scatter(dfi['MAG_i'],dfi['MAG_g']-dfi['MAG_r'],c='black',alpha=0.1)#,c=dfi['w_gmr'],cmap=cmap)
    ax[2].scatter(dfi_gr['MAG_i'],dfi_gr['MAG_g']-dfi_gr['MAG_r'],c=dfi_gr['w_gmr'],cmap=cmap)#,norm=norm)
    ax[2].errorbar(dfi_gr['MAG_i'],dfi_gr['MAG_g']-dfi_gr['MAG_r'],xerr=dfi_gr['MAGERR_i'],yerr=np.sqrt(dfi_gr['MAGERR_g']**2+dfi_gr['MAGERR_r']**2),fmt='none',c='k',alpha=0.05)
    # plt.plot(df.App_i,df.App_g-df.App_r,'.')
    # ax[2].axhline(xs[:-1][(xs[:-1]<1.65) & (xs[:-1]>np.min(x2))][0],lw=0.7,color='green')
    for z in bins_gmr_cen:
        ax[2].plot(I,linear_gmr(I,round(z,4)),color=s_m.to_rgba(z))
    ax[2].set_ylim(0.75,2)
    ax[2].set_xlim(16,24)
    # cbar=plt.colorbar(s_m)
    ax[2].set_xlabel('I')
    ax[2].set_ylabel('G-R')
    ax[2].set_title('z=0.15-0.325')
    # plt.plot([corr_f(z) for z in df.z.values[:-25]],df.App_g[:-25]-df.App_r[:-25],'-')

    xs=np.arange(np.min(z_fit)-0.1,np.max(z_fit)+0.1,0.001)
    ax[3].bar(bin_cen2, ns2, width=bin_width, color='#1f77b4') #widht = 0.025
    ax[3].bar(bin_cen1, ns1, width=bin_width, color='#ff7f0e') #width = 0.025
    ax[3].axvline(z_max,ls='--',color='purple',label='z_max:%.2f'%z_max)
    ax[3].axvline(redshift,color='red',label='z:%.2f'%redshift)
    ax[3].plot(z_fit,n_fit+n_bkg,'o',c='tab:purple')
    ax[3].plot(xs, gaussian_func(xs, *popt)+n_bkg, c='tab:green', ls='--', label='fit: a=%.2f, mu=%.4f'% tuple(popt))
    ax[3].axhline(n_median,color='tab:green',label='median:%.2f'%n_median)
    ax[3].axhline(n_mean,color='tab:red',label='mean:%.2f'%n_mean)
    ax[3].legend(loc='best')
    ax[3].set_xlabel('z')
    ax[3].set_xlim(0.1,0.7)
    ax[3].set_title('ebv:%.3f,ebv_g-r:-%.3f,ebv_r-i:-%.3f'%(ebv[0],ebv_g-ebv_r,ebv_r-ebv_i))
    if np.max(n_total)<30:
        ax[3].set_ylim(0,30)

    if duplicate:
        xs = np.arange(np.min(z_fit_dup) - 0.1, np.max(z_fit_dup) + 0.1, 0.001)
        ax[4].bar(bin_cen2, ns2-ns_dup2, width=bin_width, color='#1f77b4') #widht = 0.025
        ax[4].bar(bin_cen1, ns1-ns_dup1, width=bin_width, color='#ff7f0e') #width = 0.025
        ax[4].axvline(z_max,ls='--',color='purple',label='z_max:%.2f'%z_max)
        ax[4].axvline(redshift,color='red',label='z:%.2f'%redshift)
        ax[4].plot(z_fit_dup,n_fit_dup+n_bkg_dup,'o',c='tab:purple')
        ax[4].plot(xs, gaussian_func(xs, *popt_dup)+n_bkg_dup, c='tab:green', ls='--', label='fit: a=%.2f, mu=%.4f'% tuple(popt))
        ax[4].legend(loc='best')
        ax[4].set_xlabel('z')
        ax[4].set_xlim(0.1,0.7)
        if np.max(n_total)<30:
            ax[4].set_ylim(0,30)
        
    if colorerr:
        dfi_rmi = dfi[~np.isnan(dfi['z_rmi'])]
        dfi_gmr = dfi[~np.isnan(dfi['z_gmr'])]
        zs_gmr = np.arange(0.11, 0.3425, 0.002)
        zs_rmi = np.arange(0.3425, 0.65, 0.002)

        ntot_rmi = np.repeat(0, len(zs_rmi))
        ntot_gmr = np.repeat(0, len(zs_gmr))
        for i, row in dfi_rmi.iterrows():
            # for i, row in dfi.iterrows():
            i0 = row['MAG_i']
            rmi = row['MAG_r'] - row['MAG_i']
            rmierr = np.sqrt((row['MAGERR_r'])**2 + row['MAGERR_i']**2)
            ntot_rmi0 = scipy.stats.norm(rmi, rmierr).pdf(
                linear_rmi(i0, zs_rmi))
            ntot_rmi = ntot_rmi + ntot_rmi0 * row['w_rmi']
            ax[4].plot(zs_rmi,ntot_rmi0*row['w_rmi'],'-',color='tab:red',alpha=0.2)

        for i, row in dfi_gmr.iterrows():
            # for i, row in dfi.iterrows():
            i0 = row['MAG_i']
            gmr = row['MAG_g'] - row['MAG_r']
            gmrerr = np.sqrt((row['MAGERR_g'])**2 + row['MAGERR_r']**2)
            ntot_gmr0 = scipy.stats.norm(gmr, gmrerr).pdf(
                linear_gmr(i0, zs_gmr))
            ntot_gmr = ntot_gmr + ntot_gmr0 * row['w_gmr']
            ax[4].plot(zs_gmr,ntot_gmr0*row['w_gmr'],'-',color='tab:cyan',alpha=0.2)
        
        ax[4].plot(zs_gmr, ntot_gmr, '.')
        ax[4].plot(zs_rmi, ntot_rmi, '.')
        ax[4].axvline(z_max,ls='--',color='purple',label='z_max:%.2f'%z_max)
        ax[4].axvline(redshift,color='red',label='z:%.2f'%redshift)
        ax[4].legend(loc='best')
        ax[4].set_xlabel('z')
        ax[4].set_xlim(0.1, 0.7)
        if np.max(np.append(ntot_gmr,ntot_rmi)) < 200:
            ax[4].set_ylim(0, 200)
        n_total_cerr=np.append(ntot_gmr,ntot_rmi)
    else:
        n_total_cerr=0

    signal_final = signal_dup if duplicate else signal

    plt.tight_layout(rect=[0, 0., 1, 0.98])
    purge('/Users/taweewat/Documents/red_sequence/pisco_color_plots/',
          'redsq_richg%s_%s_all_.*_%s_tilted.png' % ('', mode, field))
    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_richg%s_%s_all_%.3f_%s_tilted.png' % ('',mode,signal_final,field), dpi=120)
    plt.close(fig)

    # fig,ax=plt.subplots(1,4,figsize=(20,5))
    # make_images(field,ax[0])
    # dfi_gmi=dfi[~np.isnan(dfi['z_gmi'])]
    # zs_gmi=np.arange(0.115,0.69,0.002)
    # ntot_gmi=np.repeat(0,len(zs_gmi))
    # for i, row in dfi_gmi.iterrows():
    #     i0 = row['MAG_i']
    #     gmi = row['MAG_g'] - row['MAG_i']
    #     gmierr = np.sqrt((row['MAGERR_g'])**2 + row['MAGERR_i']**2)
    #     ntot_gmi0 = scipy.stats.norm(gmi, gmierr).pdf(
    #         linear_gmi(i0, zs_gmi))
    #     ntot_gmi = ntot_gmi + ntot_gmi0 * row['w_gmi']
    #     ax[3].plot(zs_gmi,ntot_gmi0*row['w_gmi'],'-',color='tab:cyan',alpha=0.2)
    # ax[1].scatter(dfi['MAG_i'],dfi['MAG_g']-dfi['MAG_i'],c='black',alpha=0.1)#dfi['w_rmi'],cmap=cmap)
    # ax[1].scatter(dfi_gmi['MAG_i'],dfi_gmi['MAG_g']-dfi_gmi['MAG_i'],c=dfi_gmi['w_gmi'],cmap=cmap)
    # ax[1].errorbar(dfi_gmi['MAG_i'], dfi_gmi['MAG_g'] - dfi_gmi['MAG_i'], xerr=dfi_gmi['MAGERR_i'],
    #                yerr=np.sqrt(dfi_gmi['MAGERR_g']**2 + dfi_gmi['MAGERR_i']**2), fmt='none', c='k', alpha=0.05)
    # for z in np.arange(0.15, 0.71, bin_width):
    #     ax[1].plot(I,linear_gmi(I,z),color=s_m.to_rgba(z))
    # ax[1].set_ylim(1.0,3.5)   
    # ax[1].set_xlim(16,24)   
    # ax[1].set_xlabel('I')
    # ax[1].set_ylabel('G-I')
    # ax[1].set_title('z=0.15-0.675')
    # ns3,xs3=np.histogram(z_gmi,bins=np.arange(0.1325,0.7,0.035),weights=np.array(w_gmi)*np.array(w_col_gmi))
    # bin_cen3 = (xs3[:-1] + xs3[1:])/2
    # z_max_gmi = bin_cen3[np.where(ns3 == np.max(ns3))[0][0]]
    # n_bkg = np.mean(sorted(ns3)[2:-2]);
    # z_total_added = np.insert(
    #     np.append(bin_cen3, bin_cen3[-1] + bin_width), 0, bin_cen3[0] - bin_width)
    # n_total_added = np.insert(np.append(ns3, 0), 0, 0) - n_bkg
    # indi = np.where(n_total_added == np.max(n_total_added))[0][0]
    # z_fit = z_total_added[[indi - 1, indi, indi + 1]]; print 'z_fit', z_fit
    # n_fit = n_total_added[[indi - 1, indi, indi + 1]]; print 'n_fit', n_fit    
    # if (n_fit[0]<0.) and (n_fit[2]<0.):
    #     popt_gmi, pcov_gmi = curve_fit(gaussian_func, z_fit, [0,n_fit[1],0], p0=[n_fit[1],z_fit[1]])
    # else:
    #     popt_gmi, pcov_gmi = curve_fit(gaussian_func, z_fit,
    #                            n_fit, p0=[n_fit[1], z_fit[1]])
    # lum_factor2 = integrate.quad( f1d, -23.455, abs_mag(22.25, tuple(popt_gmi)[1]))[0]
    # density_factor2=integrate.quad(NFW_profile,0.001,cosmo.kpc_proper_per_arcmin(tuple(popt_gmi)[1]).value/1e3)[0]
    # signal_gmi = tuple(popt_gmi)[0] / (lum_factor2 * density_factor2)
    # z_max_fit_gmi = tuple(popt_gmi)[1]    
    # ax[2].bar(bin_cen3, ns3, width = 0.035, color='#1f77b4')#, alpha=0.5)
    # ax[2].axvline(z_max_gmi, ls='--', color='purple',
    #               label='z_max=%.3f'%z_max_gmi)
    # ax[2].axvline(z_max_fit_gmi, ls='--', color='tab:green',
    #               label='z_max_fit=%.3f'%z_max_fit_gmi)
    # ax[2].axvline(redshift,color='red',label='z:%.3f'%redshift)    
    # ax[2].plot(z_fit,n_fit+n_bkg,'o',c='tab:purple')
    # xs=np.arange(np.min(z_fit)-0.1,np.max(z_fit)+0.1,0.001)
    # ax[2].plot(xs, gaussian_func(xs, *popt_gmi) + n_bkg, c='tab:green',
    #            ls='--', label='fit: a=%.2f, mu=%.4f' % tuple(popt_gmi))
    # ax[2].legend(loc='best')
    # ax[2].set_xlabel('z')
    # ax[2].set_xlim(0.1,0.7)
    # if np.max(n_total)<30:
    #     ax[2].set_ylim(0,30)
    # ax[3].plot(zs_gmi,ntot_gmi,'.')
    # ax[3].set_xlabel('z')
    # ax[3].set_xlim(0.1,0.7)
    # ax[3].axvline(z_max_fit_gmi,ls='--',color='purple',label='z_max_fit:%.2f'%z_max_fit_gmi)
    # ax[3].axvline(redshift,color='red',label='z:%.2f'%redshift)
    # if np.max(ntot_gmi)<70:
    #     ax[3].set_ylim(0,70)
    # ntot_gmi_max=np.max(ntot_gmi)
    # zs_gmi_max=zs_gmi[np.argmax(ntot_gmi)]
    # ax[3].axvline(zs_gmi_max,ls='--',color='pink',label='zs_gmi_max:%.2f'%zs_gmi_max)
    # plt.tight_layout(rect=[0, 0., 1, 0.98])
    # plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_gmi_%s_all_%.3f_%s_tilted.png' %
    #             (mode, signal_gmi, field), dpi=120)
    # plt.close(fig)

    # transparent=False
    if transparent:
        fig,ax=plt.subplots(figsize=(7,4))
        ax.bar(bin_cen2, ns2, width=0.035, color='#1f77b4') #widht = 0.025
        ax.bar(bin_cen1, ns1, width = 0.035, color='#ff7f0e') #width = 0.025
        ax.axvline(z_max,ls='--',color='purple',label='z_max:%.2f'%z_max)
        ax.set_xlabel('z')
        ax.set_xlim(0.1,0.7)
        if np.max(n_total)<30:
            ax.set_ylim(0,30)
        for axp in ax.spines:
            ax.spines[axp].set_color('white') 
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        purge('/Users/taweewat/Documents/red_sequence/pisco_color_plots/',
            'redsq_transparent_%.3f_%s_tilted.png' % (signal_final,field)) 
        plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_transparent_%.3f_%s_tilted.png' % (signal_final,field), dpi=120, transparent=True)
        plt.close(fig)

    red_dir='/Users/taweewat/Documents/red_sequence/'
    rich_filename = 'all_richness_%s.csv'%extra_name
    if not os.path.isfile(red_dir + rich_filename):
        os.system("cp %s %s"%(red_dir+'all_richness_gremove_lum_silk_zf2.5.csv',red_dir+rich_filename))

    df_richness=pd.read_csv(red_dir+rich_filename)
    df_richness=df_richness.copy()
    df_richness.loc[df_richness['name'] == field, 'Nmax'] = np.max(n_total)
    df_richness.loc[df_richness['name'] == field, 'Nbkg_mean'] = np.mean(n_total)
    df_richness.loc[df_richness['name'] == field, 'Nbkg_median'] = np.median(n_total)
    df_richness.loc[df_richness['name'] == field, 'zmax'] = z_max
    df_richness.loc[df_richness['name'] == field, 'amp'] = signal_final
    df_richness.loc[df_richness['name'] == field, 'zmax_fit'] = z_max_fit
    df_richness.loc[df_richness['name'] == field, 'gremove'] = nog
    df_richness.loc[df_richness['name'] == field, 'lum_factor'] = lum_factor
    df_richness.loc[df_richness['name'] == field, 'density_factor'] = density_factor
    # df_richness.loc[df_richness['name'] == field, 'amp_gmi'] = signal_gmi
    # df_richness.loc[df_richness['name'] == field, 'z_max_fit_gmi'] = z_max_fit_gmi
    # df_richness.loc[df_richness['name']==field,'distance[Mpc]']=z_max*(4000)
    # df_richness.loc[df_richness['name']==field,'R']=(np.max(n_total)-np.mean(n_total))*(z_max*4000)**2
    df_richness.to_csv(red_dir+rich_filename,index=0)

    dfi.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/galaxy_%s_final_%s.csv"%(mode,field))

    # get member redshfit in the figure
    if img_redshift:
        image_redshift(field,signal,tuple(popt)[1],mode)
    # get total images with red-sequence
    if img_filp:
        image_flip(field,signal,tuple(popt)[1],mode)

    if colorerr:
        return z_total, n_total, n_total_cerr
    else:
        return z_total, n_total, n_total_dup


def pisco_combine_imgs(fields, mode='psf', mode2mass=''):
    dir1='/Users/taweewat/Documents/red_sequence/pisco_color_plots/psf_est/'
    dir2='/Users/taweewat/Documents/red_sequence/pisco_color_plots/'
    dir3='/Users/taweewat/Documents/red_sequence/pisco_color_plots/'
    dirout='/Users/taweewat/Documents/red_sequence/pisco_all/'

    myReg = re.compile(r'(redsq_richg_%s_all_.*%s.*png)' % (mode, field))
    myReg2=re.compile(r'(\d{1,3}\.\d{1,3})')
    names=[]
    for text in os.listdir(dir3):
        if myReg.search(text) != None:
            names.append(myReg.search(text).group())

    if names==[]:
        print 'no files', field

    img1=dir1+'psf_est3_'+field+'_i.png'
    img2=dir2+'star_galaxy_sep_12_all'+field+'.png'
    img3=dir3+names[0]

    images_list=[img1, img2, img3]
    imgs=[]
    try:
        imgs = [ Image_PIL.open(i) for i in images_list ]
    except:
        print 'no image file', field
    mw = imgs[2].width/2
    h = imgs[0].height+imgs[1].height/2+imgs[2].height/2
    result = Image_PIL.new("RGBA", (mw, h))
    y,index=0,0
    for i in imgs:
        if (index==2) or (index==1):
            i=i.resize((i.width/2,i.height/2))
        result.paste(i, (0, y))
        y += i.size[1]
        index+=1
    # result.save(dirout + 'all_combine%s_%s_%s_%s.png' %
    #             (mode2mass, field, mode, myReg2.search(names[0]).group())) 
    # result.save(dirout + 'all_combine_%s_%s_%s_%s_%s.png' %
    #             (extra_name, mode2mass, myReg2.search(names[0]).group(), mode, field)) 

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))

def image_redshift(field,signal,redshift,mode):
    df_total=pd.read_csv('/Users/taweewat/Documents/pisco_code/slr_output/galaxy_%s_final_%s.csv'%(mode,field),index_col=0)
    df_star=pd.read_csv('/Users/taweewat/Documents/pisco_code/slr_output/star_psf_total_%s.csv'%field,index_col=0)
    # df_star=df_star[df_star['SG']>0.95]
    hdu=fits.open('/Users/taweewat/Documents/red_sequence/panstar/coadd_panstar_%s_i.fits'%field)
    img=hdu[0].data.astype(float)
    img-=np.median(img.ravel()[~np.isnan(img.ravel())])    
    def redshift_f(row):
        if not np.isnan(row['z_gmr']):
            redshift=row['z_gmr']
        if not np.isnan(row['z_rmi']):
            redshift=row['z_rmi']
        if np.isnan(row['z_rmi']) and np.isnan(row['z_gmr']):
            redshift=0
        return redshift
        
    df_total['redshift_m']=df_total.apply(lambda row: redshift_f(row), axis=1)
    def size_f(row):
        if not np.isnan(row['w_gmr']):
            size=row['w_gmr']
        if not np.isnan(row['w_rmi']):
            size=row['w_rmi']
        if np.isnan(row['w_rmi']) and np.isnan(row['w_gmr']):
            size=0
        return size
    df_total['size_m']=df_total.apply(lambda row: size_f(row), axis=1)

    df_total0=df_total.copy()
    df_total=df_total[df_total['redshift_m'] > 0].copy()

    norm = matplotlib.colors.Normalize(vmin=-2,vmax=4)
    c_m = matplotlib.cm.Greys_r
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    normalize = matplotlib.colors.Normalize(vmin=0.1, vmax=0.7)


    fig, (a0, a1) = plt.subplots(1,2, figsize=(30,18), gridspec_kw = {'width_ratios':[0.8, 1]})
    # a0.imshow(img, cmap=c_m, norm=norm, origin='lower')
    # a0.scatter(df_star['XWIN_IMAGE_i'].values,df_star['YWIN_IMAGE_i'].values,s=100, marker='*', facecolors='none', edgecolors='yellow', label='star')
    # df1i=df_total[df_total['w_rmi']>0.1]
    # df2i=df_total[df_total['w_rmi']<=0.1]
    # # a0.scatter(df1i['XWIN_IMAGE_i'].values,df1i['YWIN_IMAGE_i'].values,s=100, facecolors='none', edgecolors='blue')
    # a0.scatter(df1i['XWIN_IMAGE_i'].values, df1i['YWIN_IMAGE_i'].values, s=100, c=df1i['size_m'].values, cmap='RdYlGn')
    # a0.scatter(df2i['XWIN_IMAGE_i'].values,df2i['YWIN_IMAGE_i'].values, s=100, facecolors='none', edgecolors='white')
    # a0.set_xlim(0,1600)
    # a0.set_ylim(0, 2250)
    try:
        img2 = mpimg.imread('/Users/taweewat/Documents/pisco_code/Chips_images/aplpy_panstar_%s_img.jpeg' % field)
    except:
        img2 = mpimg.imread('/Users/taweewat/Documents/pisco_code/Chips_images/aplpy_panstar_%s_img4.jpeg' % field)
    imgplot = a0.imshow(img2)
    a0.axis('off')
    a0.annotate('Redshift: %.3f\nRichness: %.2f' %
                (redshift, signal), xy=(150, 100), color='white')

    a1.imshow(img, cmap=c_m, norm=norm, origin='lower')

    a1.scatter(df_star['XWIN_IMAGE_i'].values,df_star['YWIN_IMAGE_i'].values, s=300,edgecolor='orange', facecolor='none',lw=3)
    a1.scatter(df_total0['XWIN_IMAGE_i'].values,df_total0['YWIN_IMAGE_i'].values, s=150,edgecolor='tab:blue', facecolor='none',lw=1,alpha=0.5)
    #,s=100, marker='*', facecolors='none', edgecolors='yellow', label='star')
    axi = a1.scatter(df_total['XWIN_IMAGE_i'].values, df_total['YWIN_IMAGE_i'].values,
                     s=(df_total['size_m'].values * 200)+30, c=df_total['redshift_m'].values, cmap='tab20b', norm=normalize)
    plt.colorbar(axi)  # df_total['size_m'].values*300
    # a1.set_xlim(0, 1600)
    # a1.set_ylim(0, 2250)
    plt.tight_layout()

    left, bottom, width, height = [0.05, 0.24, 0.3, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.imshow(mpimg.imread(
        '/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_transparent_%.3f_%s_tilted.png' % (signal, field)))
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.axis('off')

    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_image_redshift/img_redshift_%s_%.3f_%s.png' %
                (mode, signal,field), dpi=50)
    plt.close(fig)


def image_flip(field, signal, redshift, mode):
    img = mpimg.imread(
        '/Users/taweewat/Documents/pisco_code/Chips_images/aplpy4_%s_img.jpeg' % field)
    fig, ax = plt.subplots(figsize=(7, 7))
    imgplot = ax.imshow(img)
    ax.axis('off')
    ax.annotate('Redshift: %.3f\nRichness: %.2f' %
                (redshift, signal), xy=(150, 100), color='white')

    left, bottom, width, height = [0.2, 0.18, 0.3, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.imshow(mpimg.imread(
        '/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_transparent_%.3f_%s_tilted.png' % (signal, field)))
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.axis('off')

    # plt.tight_layout()
    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_image_redshift/image_flip_%s_%.3f_%s.png' %
                (mode, signal, field), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    """
    execute:
    python pisco_pipeline/pisco_photometry_all.py CHIPS111 psf slr
    #updated version with no2mass option for no more comparison with known 2mass stars
    python pisco_pipeline/pisco_photometry_all.py CHIPS111 psf allslr no2mass
    """
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    field = str(sys.argv[1])
    mode = str(sys.argv[2])   #aper, psf, auto, hybrid
    all_argv=sys.argv[3:]  #allslr, slr, noslr
    if (all_argv[0]=='allslr') | (all_argv[0]=='slr'):
        slr=str(all_argv[0])
        slr_param=True
    elif all_argv[0]=='noslr':
        slr='no_slr'
        slr_param=False
    if all_argv[1]=='2mass':
        mode2mass=''
    elif all_argv[1]=='no2mass':
        mode2mass='_no2mass'

    home='/Users/taweewat/Documents/pisco_code/' #09, 171208
    # dirs=['ut170103/','ut170104/','ut170619/','ut170621/','ut170624/','ut171208/','ut171209/','ut171212/']
    # 'ut171208/', 'ut171209/','ut171212/', 'ut170621/', 'ut170624/'
    dirs = ['ut170619/']
    # dirs = ['ut171209/']
    names=[]
    myReg=re.compile(r'(CHIPS\d{4}[+-]\d{4})|(Field\d{3})')
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                names.append(myReg.search(text).group())
    all_fields=list(set(names))

    exception = ['CHIPS0525-6938', 'Field234']
    z_total_all,n_total_all,n_total_dup_all=[],[],[]

    # all_fields=['CHIPS1911+4455']
    all_fields=['CHIPS1011-0505']
    all_fields_cut = all_fields[:]
    notgoflag=True
    for index, field in enumerate(all_fields_cut):
        print field, '%i/%i' % (index, len(all_fields_cut))

        # if field == 'CHIPS0122-2646':
        #     notgoflag = False; continue
        # if notgoflag:
        #     continue 
        if field in exception:
            continue

        if slr=='allslr':
            # star_galaxy_bleem(field)
            pisco_photometry_v4(field)
            # pisco_cut_frame(field)
        elif slr=='slr':
            # star_galaxy_bleem(field)
            # pisco_photometry_v4(field)
            panstar_cut_star(field)
            # pisco_cut_frame(field)
            pisco_photometry_psf_v4(field, mode=mode, mode2mass=mode2mass, slr=slr_param)

            purge('/Users/taweewat/Documents/red_sequence/pisco_color_plots/'\
                ,r'(redsq_%s_all_.*%s.*png)'%(mode,field))
            z_total,n_total,n_total_dup=pisco_tilt_resequence(field, mode=mode, mode2mass=mode2mass)
            z_total_all.append(z_total)
            n_total_all.append(n_total)
            n_total_dup_all.append(n_total_dup)
            # pisco_combine_imgs(field, mode=mode, mode2mass=mode2mass)
            pickle.dump( [z_total_all,n_total_all,n_total_dup_all], open( "pickle_all_richness_%s.pickle"%extra_name, "wb" ) )
            print 'save pickle fie at', "pickle_all_richness_%s.pickle" % extra_name
        elif slr == 'no_slr':
            # pisco_cut_frame(field)
            panstar_cut_star(field)
            pisco_photometry_psf_v4(field, mode=mode, mode2mass=mode2mass, slr=slr_param)
            purge('/Users/taweewat/Documents/red_sequence/pisco_color_plots/'\
                ,r'(redsq_%s_all_.*%s.*png)'%(mode,field))
            z_total,n_total,n_total_dup=pisco_tilt_resequence(field, mode=mode, mode2mass=mode2mass)
            z_total_all.append(z_total)
            n_total_all.append(n_total)
            n_total_dup_all.append(n_total_dup)
            # pisco_combine_imgs(field, mode=mode, mode2mass=mode2mass)
            pickle.dump( [z_total_all,n_total_all,n_total_dup_all], open( "pickle_all_richness_%s.pickle"%extra_name, "wb" ) )
            print 'save pickle fie at', "pickle_all_richness_%s.pickle" % extra_name
    
        purge('final', "proj_coadd_c%s_.*\.fits" % field)
        purge('.', "proto_psf_%s_.*\.fits" % field)
        purge('.', "samp_psf_%s_.*\.fits" % field)
        purge('.', "resi_psf_%s_.*\.fits" % field)
        purge('.', "snap_psf_%s_.*\.fits" % field)
        purge('.', "chi_psf_%s_.*\.fits" % field)
        # purge('psfex_output', "psf_%s_.*\.fits" % field)
        # purge('slr_output', "a_psf_%s_.*\.fits" % field)
        purge('final', "coadd_c%s_sq_.*\.fits" % field)

        
