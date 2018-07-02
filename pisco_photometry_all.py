import sys, os, re, yaml, subprocess, shlex, FITS_tools
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import image
import matplotlib.cm as cm

from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy import interpolate

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
Example: python pisco_pipeline/pisco_photometry_all.py PKS1353 psf allslr 2mass
field: name of the fields
mode: psf, auto, aper, hybrid
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
    return seeing
def read_param():
    with open("pisco_pipeline/params.yaml", 'r') as stream:
        try:
            param=yaml.load(stream)
            return param
        except yaml.YAMLError as exc:
            print(exc)
def read_param_izp():
    with open("/Users/taweewat/Documents/pisco_code/pisco_pipeline/params_izeropoint.yaml", 'r') as stream:
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
    fits.writeto('final/coadd_c%s_sq_i.fits'%field, data2, header=header, overwrite=True)
    cmd='sex final/coadd_c%s_i.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -CATALOG_TYPE FITS_1.0 -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE 0.22 -DETECT_MINAREA %s -CHECKIMAGE_NAME checki.fits,segmenti.fits'%\
    (field,'sex.param',sg_dir+'/%s_catalog.fits'%(field),str(seeing),str(param['satur_level_i']),str(1.1/minarea*np.pi*(seeing/0.22)**2)); print cmd
    subprocess.check_call(shlex.split(cmd))
    cmd='sex final/coadd_c%s_i.fits,final/coadd_c%s_sq_i.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -CATALOG_TYPE FITS_1.0 -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE 0.22 -DETECT_MINAREA %s'%\
    (field,field,'sex.param',sg_dir+'/%s_sq_catalog.fits'%(field),str(seeing),str(param['satur_level_i']),str(1.1/minarea*np.pi*(seeing/0.22)**2)); print cmd
    subprocess.check_call(shlex.split(cmd))

def pisco_photometry_v4(field):
    def aperature_proj(field,band):
        param=read_param()
        seeing=find_seeing(field,band)
        # seeing=1.5

        slrdir = 'slr_output'
        to_be_projected = 'final/coadd_c%s_%s.fits'%(field,band)
        reference_fits  = 'final/coadd_c%s_i.fits'%field
        im1,im2, header = FITS_tools.match_fits(to_be_projected,reference_fits,return_header=True)
        outname = 'final/proj_coadd_c%s_%s.fits'%(field,band)
        print 'projecting from %s band to i band the fits file '%band + outname
        fits.writeto(outname, im1, header, overwrite=True)

        minarea=1.7
        cmd='sex final/coadd_c%s_%s.fits -c pisco_pipeline/config.sex -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PHOT_APERTURES 15 -PIXEL_SCALE 0.22 -DETECT_MINAREA %s'%\
        (field,band,'sex_psf.param','psfex_output/psf_%s_%s.fits'%(field,band),str(seeing),str(param['satur_level_%s'%band]),str(1.1/minarea*np.pi*(seeing/0.22)**2))
        print cmd
        subprocess.check_call(shlex.split(cmd))

        cmd='psfex %s -c pisco_pipeline/pisco.psfex' % ('psfex_output/psf_%s_%s.fits'%(field,band))
        print cmd
        subprocess.check_call(shlex.split(cmd))

        cmd='sex final/coadd_c%s_i.fits,final/proj_coadd_c%s_%s.fits -c pisco_pipeline/config.sex -PSF_NAME %s -PARAMETERS_NAME pisco_pipeline/%s -CATALOG_NAME %s -SEEING_FWHM %s -SATUR_LEVEL %s -PIXEL_SCALE 0.22 -CATALOG_TYPE FITS_1.0 -PHOT_APERTURES 15 -DETECT_MINAREA %s -CHECKIMAGE_NAME check%s.fits,segment%s.fits'%\
        (field,field,band,'psfex_output/psf_%s_%s.psf'%(field,band),'sex_after_psf.param','%s/a_psf_%s_%s.fits'%(slrdir,field,band),str(seeing),str(param['satur_level_%s'%band]),str(1.1/1.7*np.pi*(seeing/0.22)**2),band,band)
        print cmd
        subprocess.check_call(shlex.split(cmd))

        table=Table.read('%s/a_psf_%s_%s.fits'%(slrdir,field,band))
        for name in table.colnames[1:]:
            table.rename_column(name, name + '_%s' % band)
        return table

    slrdir = 'slr_output'
    if not os.path.exists(slrdir):
        os.makedirs(slrdir)

    tableg=aperature_proj(field,'g')
    tablei=aperature_proj(field,'i')
    tabler=aperature_proj(field,'r')
    tablez=aperature_proj(field,'z')

    print 'len of all table', len(tableg), len(tablei), len(tabler), len(tablez)

    ci=SkyCoord(ra=np.array(tablei['ALPHA_J2000_i'])*u.degree, dec=np.array(tablei['DELTA_J2000_i'])*u.degree)# print len(ci)
    cg=SkyCoord(ra=np.array(tableg['ALPHA_J2000_g'])*u.degree, dec=np.array(tableg['DELTA_J2000_g'])*u.degree)# print len(cg)
    cr=SkyCoord(ra=np.array(tabler['ALPHA_J2000_r'])*u.degree, dec=np.array(tabler['DELTA_J2000_r'])*u.degree)# print len(cr)
    cz=SkyCoord(ra=np.array(tablez['ALPHA_J2000_z'])*u.degree, dec=np.array(tablez['DELTA_J2000_z'])*u.degree)# print len(cz)

    idxn, d2dn, d3dn=cg.match_to_catalog_sky(ci)
    Table_I=tablei[idxn][['ALPHA_J2000_i','DELTA_J2000_i','MAG_APER_i','MAGERR_APER_i','MAG_AUTO_i','MAGERR_AUTO_i','MAG_HYBRID_i','MAGERR_HYBRID_i',\
                  'CLASS_STAR_i','FLAGS_i','MAG_PSF_i','MAGERR_PSF_i','MAG_MODEL_i','MAGERR_MODEL_i','SPREAD_MODEL_i']]
    Table_I.rename_column('ALPHA_J2000_i','ALPHA_J2000')
    Table_I.rename_column('DELTA_J2000_i','DELTA_J2000')

    idxn, d2dn, d3dn=cg.match_to_catalog_sky(cr)
    Table_R=tabler[idxn][['ALPHA_J2000_r','DELTA_J2000_r','MAG_APER_r','MAGERR_APER_r','MAG_AUTO_r','MAGERR_AUTO_r','MAG_HYBRID_r','MAGERR_HYBRID_r',\
                  'CLASS_STAR_r','FLAGS_r','MAG_PSF_r','MAGERR_PSF_r','MAG_MODEL_r','MAGERR_MODEL_r','SPREAD_MODEL_r']]
    Table_R.rename_column('ALPHA_J2000_r','ALPHA_J2000')
    Table_R.rename_column('DELTA_J2000_r','DELTA_J2000')

    idxn, d2dn, d3dn=cg.match_to_catalog_sky(cz)
    Table_Z=tablez[idxn][['ALPHA_J2000_z','DELTA_J2000_z','MAG_APER_z','MAGERR_APER_z','MAG_AUTO_z','MAGERR_AUTO_z','MAG_HYBRID_z','MAGERR_HYBRID_z',\
                  'CLASS_STAR_z','FLAGS_z','MAG_PSF_z','MAGERR_PSF_z','MAG_MODEL_z','MAGERR_MODEL_z','SPREAD_MODEL_z']]
    Table_Z.rename_column('ALPHA_J2000_z','ALPHA_J2000')
    Table_Z.rename_column('DELTA_J2000_z','DELTA_J2000')

    Table_G=tableg[['ALPHA_J2000_g','DELTA_J2000_g','MAG_APER_g','MAGERR_APER_g','MAG_AUTO_g','MAGERR_AUTO_g','MAG_HYBRID_g','MAGERR_HYBRID_g',\
                  'CLASS_STAR_g','FLAGS_g','MAG_PSF_g','MAGERR_PSF_g','MAG_MODEL_g','MAGERR_MODEL_g','SPREAD_MODEL_g']]
    Table_G.rename_column('ALPHA_J2000_g','ALPHA_J2000')
    Table_G.rename_column('DELTA_J2000_g','DELTA_J2000')

    print 'len of all new table', len(Table_G), len(Table_I), len(Table_R), len(Table_Z)

    total=join(join(join(Table_I,Table_G,keys=['ALPHA_J2000','DELTA_J2000']),Table_R,keys=['ALPHA_J2000','DELTA_J2000']),\
         Table_Z,keys=['ALPHA_J2000','DELTA_J2000'])

    # total=join(join(join(mag_ii,mag_ig,keys='NUMBER'), mag_ir,keys='NUMBER'),\
    #            mag_iz,keys='NUMBER')
    #'NUMBER'
    total2=total[['ALPHA_J2000','DELTA_J2000',\
                  'MAG_APER_i','MAGERR_APER_i','MAG_APER_g','MAGERR_APER_g','MAG_APER_r',\
                  'MAGERR_APER_r','MAG_APER_z','MAGERR_APER_z','MAG_AUTO_i','MAGERR_AUTO_i',\
                  'MAG_AUTO_g','MAGERR_AUTO_g','MAG_AUTO_r','MAGERR_AUTO_r','MAG_AUTO_z',\
                  'MAGERR_AUTO_z','MAG_HYBRID_i','MAGERR_HYBRID_i','MAG_HYBRID_g',\
                  'MAGERR_HYBRID_g','MAG_HYBRID_r','MAGERR_HYBRID_r','MAG_HYBRID_z',\
                  'MAGERR_HYBRID_z','CLASS_STAR_i','CLASS_STAR_g','CLASS_STAR_r',\
                  'CLASS_STAR_z','FLAGS_g','FLAGS_r','FLAGS_i','FLAGS_z','MAG_PSF_g',\
                  'MAG_PSF_r','MAG_PSF_i','MAG_PSF_z','MAGERR_PSF_g','MAGERR_PSF_r',\
                  'MAGERR_PSF_i','MAGERR_PSF_z','MAG_MODEL_g','MAG_MODEL_r',\
                  'MAG_MODEL_i','MAG_MODEL_z','MAGERR_MODEL_g','MAGERR_MODEL_r',\
                  'MAGERR_MODEL_i','MAGERR_MODEL_z','SPREAD_MODEL_g','SPREAD_MODEL_r',\
                  'SPREAD_MODEL_i','SPREAD_MODEL_z',]]

    total2.write(os.path.join(slrdir, 'total_psf_%s.csv' % field), overwrite=True)
    total2.write(slrdir+'/all_psf_%s.fits' % field, overwrite=True)

def pisco_cut_star(field,c_a,c_b,c_d,c_delta):
    seeing=find_seeing_fits(field)
    true_seeing=find_seeing(field,'i')

    df_i=Table(fits.open('/Users/taweewat/Documents/pisco_code/star_galaxy/%s_catalog.fits'%field)[1].data).to_pandas()
    df_isq=Table(fits.open('/Users/taweewat/Documents/pisco_code/star_galaxy/%s_sq_catalog.fits'%field)[1].data).to_pandas()

    #cut the object out so that it has the same number of object between the sq catalog list and the psf mag list.
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/total_psf_%s.csv"%field
    df0 = pd.read_csv(fname)
    cf_i=SkyCoord(ra=np.array(df_i['ALPHA_J2000'])*u.degree, dec=np.array(df_i['DELTA_J2000'])*u.degree)
    cf_isq=SkyCoord(ra=np.array(df_isq['ALPHA_J2000'])*u.degree, dec=np.array(df_isq['DELTA_J2000'])*u.degree)
    cf0=SkyCoord(ra=np.array(df0['ALPHA_J2000'])*u.degree, dec=np.array(df0['DELTA_J2000'])*u.degree)
    idxn, d2dn, d3dn=cf0.match_to_catalog_sky(cf_i)
    df_i_cut=df_i.loc[idxn].copy()
    df_i_cut['NUMBER']=np.arange(0,len(df0),1).tolist()
    idxn, d2dn, d3dn=cf0.match_to_catalog_sky(cf_isq)
    df_isq_cut=df_isq.loc[idxn].copy()
    df_isq_cut['NUMBER']=np.arange(0,len(df0),1).tolist()

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
    [np.argmin(df_i2[f(df_i2.MAG_APER)-(df_i2.MAG_APER-df_isq2.MAG_APER)<0.1]['MAG_APER'].values)]
    #--->

    ax[0,0].axvline(c_c,color='red',label='new upper cut')
    ax[0,0].legend(loc='best')

    ax[0,1].scatter(df_i0.MAG_APER,df_i0.MAG_APER-df_isq0.MAG_APER,marker='.',c=df_i0.CLASS_STAR)
    ax[0,1].plot(df_i3.MAG_APER,df_i3.MAG_APER-df_isq3.MAG_APER,'x')
    ax[0,1].set_title('for all objects')
    ax[0,1].set_ylabel('MAG_APER-MAG_APER_sq')
    ax[0,1].set_xlabel('MAG APER i')
    ax[0,1].axvline(c_b,ls='--')
    ax[0,1].axvline(c_c,ls='--')

    delta=(df_i0.MAG_APER-df_isq0.MAG_APER) - f(df_i0.MAG_APER)

    ax[0,2].scatter(df_i0.MAG_APER,delta,marker='.',c=df_i0.CLASS_STAR)
    ax[0,2].axhline(0,ls='--')
    ax[0,2].axvline(c_c,ls='--')
    ax[0,2].axvline(c_b,ls='--')
    ax[0,2].set_ylabel('Delta')
    ax[0,2].set_xlabel('MAG APER i')
    ax[0,2].set_ylim(0.5,-1.2)

    df_i1=df_i0[(df_i0.MAG_APER>c_c)&(df_i0.MAG_APER<c_b)].copy()
    df_isq1=df_isq0[(df_i0.MAG_APER>c_c)&(df_i0.MAG_APER<c_b)].copy()
    delta1=(df_i1.MAG_APER-df_isq1.MAG_APER) - f(df_i1.MAG_APER)
    ax[1,0].scatter(df_i1.MAG_APER,delta1,marker='o',c=df_i1.CLASS_STAR)
    ax[1,0].axhline(0,ls='--')
    ax[1,0].axhline(-0.1,ls='--')
    ax[1,0].set_ylabel('Delta')
    ax[1,0].set_xlabel('MAG APER i')
    ax[1,0].set_ylim(0.5,-2)

    deltag=delta1[delta1<c_delta] #galaxy  0.1, 0.2 (0.005), 0.5 ()
    deltas=delta1[(delta1>=c_delta)&(delta1<3.)] #star

    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))

    #galaxy
    hist, bin_edges = np.histogram(deltag,bins=np.arange(-1.2,0.5,0.02))#, density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    ax[1,1].plot(bin_centres, hist, label='galaxies',linestyle='steps')

    #stars
    hist, bin_edges = np.histogram(deltas,bins=np.arange(-1,0.5,0.02))
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    p0 = [1., 0., 0.1]
    coeff2, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    x=np.arange(min(bin_edges),max(bin_edges),0.02)
    hist_fit2 = gauss(x, *coeff2)
    ax[1,1].plot(bin_centres, hist, label='stars',linestyle='steps')
    ax[1,1].plot(x, hist_fit2, label='stars_fit')
    ax[1,1].legend(loc='best')
    ax[1,1].set_xlabel('Delta')
    ax[1,1].set_ylabel('Histogram')

    maxi=np.max(gauss(delta,*coeff2))
    def prob_SG(delta,maxi,*coeff2):
        if delta>0.:
            return 0.
        elif delta<=0.:
            return 1.-(gauss(delta,*coeff2)/maxi)

    vprob_SG= np.vectorize(prob_SG)
    SG=1.-vprob_SG(delta1,maxi,*coeff2)
    df_i1.loc[:,'SG']=SG
    ax[1,2].scatter(df_i1.MAG_APER,SG,marker='.',c=df_i1.CLASS_STAR)
    ax[1,2].set_ylim(-0.02,1.02)
    ax[1,2].set_xlabel('MAG APER i')
    ax[1,2].set_ylabel('SG (probability to be a star)')
    plt.suptitle(field+' seeing vs true_seeing: '+str(seeing)+','+str(true_seeing))
    plt.tight_layout(rect=[0, 0., 1, 0.98])
    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/star_galaxy_sep_12_all%s.png' % field, dpi=72)
    plt.close(fig)
    return df_i_cut, df_i1

def pisco_cut_frame(field):
    # df_i=Table(fits.open('/Users/taweewat/Documents/pisco_code/star_galaxy/'+
    #                      '%s_catalog.fits'%field)[1].data).to_pandas()
    seeing=find_seeing_fits(field)
    true_seeing=find_seeing(field,'i')

    if field=='Field179':
        true_seeing=1.12

    if field=='CHIPS1011-0505':
        c_a,c_b,c_c,c_d,c_delta=[0.95,-9.,-14,-10,-0.3]
    elif (field[0:3]=='PKS') or (field[0:4]=='SDSS'):
        c_a,c_b,c_c,c_d,c_delta=[0.95,-7,-12.2,-8.5,-0.1]
    elif field[0:5]=='CHIPS':
        c_a,c_b,c_c,c_d,c_delta=[0.95,-9.,-12.5,-10,-0.15] #-0.15 (default)
    elif field[0:5]=='Field':
        c_a,c_b,c_c,c_d,c_delta=[0.95,-8.,-11.5,-8.5,-0.1] #-8.5, -0.1 (Default)

    df_i_cut, df_i1=pisco_cut_star(field,c_a,c_b,c_d,c_delta)
    # print type(df_i1)
    while len(df_i1[(df_i1['MAG_APER']<c_b)&(df_i1['MAG_APER']>c_b-0.5)\
                    &(df_i1['SG']>0.2)&(df_i1['SG']<0.8)])<15: #8, 10 (default)
        len_df_i1=len(df_i1)
        c_b=c_b+0.5
        df_i_cut, df_i1 = pisco_cut_star(field,c_a,c_b,c_d,c_delta)
        if len_df_i1==len(df_i1):
            break
    dff=df_i1[df_i1['SG']<0.1] #0.8
    dff_star=df_i1[df_i1['SG']>0.1] #0.9 vs 0.2

    #After running pisco_photometry_v4.fits, but before running
    #pisco_photometry_psf_v4.fits, and then 19_pisco_tilt_resequence.ipynb
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/total_psf_%s.csv"%field
    df0 = pd.read_csv(fname)
    df0['NUMBER']=np.arange(0,len(df0),1).tolist()  #add NUMBER parameter to match between cut catalog and the total catalog
    df0.rename(columns={'ALPHA_J2000': 'ALPHA_J2000_i'}, inplace=True)
    df0.rename(columns={'DELTA_J2000': 'DELTA_J2000_i'}, inplace=True)
    print len(df0), len(df_i_cut), len(dff), len(dff_star), '=', seeing, 'vs [NEW]', true_seeing
    if len(df0)!=len(df_i_cut):
        raise ValueError('the two tables do not have the same length')
    dff0=pd.merge(dff,df0,on='NUMBER')
    dff0.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/"+\
                "galaxy_psf_total_%s.csv"%field)
    dff_star0=pd.merge(dff_star, df0, on='NUMBER')
    dff_star0.to_csv("/Users/taweewat/Documents/pisco_code/slr_output/"+\
                     "star_psf_total_%s.csv"%field)

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
        with open(fname) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        if len(content)==8:
            red_content=content[4:]
        elif len(content)==10:    
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
        table['MAG_' + band[3]] = table['MAG_%s_'%MODE1 + band[3]] + corr[3]
        table['MAGERR_' + band[0]] = (table['MAGERR_%s_'%MODE1 + band[0]]**2 + ecorr[0]**2)**0.5
        table['MAGERR_' + band[1]] = (table['MAGERR_%s_'%MODE1 + band[1]]**2 + ecorr[1]**2)**0.5
        table['MAGERR_' + band[2]] = (table['MAGERR_%s_'%MODE1 + band[2]]**2 + ecorr[2]**2)**0.5
        table['MAGERR_' + band[3]] = (table['MAGERR_%s_'%MODE1 + band[3]]**2 + ecorr[3]**2)**0.5
        return table

    slrdir = 'slr_output'
    total3=Table.from_pandas(pd.read_csv("/Users/taweewat/Documents/pisco_code/slr_output/star_psf_total_%s.csv"%field))
    total3=total3[['NUMBER','ALPHA_J2000_i','DELTA_J2000_i',\
                  'MAG_APER_i','MAGERR_APER_i','MAG_APER_g','MAGERR_APER_g','MAG_APER_r',\
                  'MAGERR_APER_r','MAG_APER_z','MAGERR_APER_z','MAG_AUTO_i','MAGERR_AUTO_i',\
                  'MAG_AUTO_g','MAGERR_AUTO_g','MAG_AUTO_r','MAGERR_AUTO_r','MAG_AUTO_z',\
                  'MAGERR_AUTO_z','MAG_HYBRID_i','MAGERR_HYBRID_i','MAG_HYBRID_g',\
                  'MAGERR_HYBRID_g','MAG_HYBRID_r','MAGERR_HYBRID_r','MAG_HYBRID_z',\
                  'MAGERR_HYBRID_z','CLASS_STAR_i','CLASS_STAR_g','CLASS_STAR_r',\
                  'CLASS_STAR_z','FLAGS_g','FLAGS_r','FLAGS_i','FLAGS_z','MAG_PSF_g',\
                  'MAG_PSF_r','MAG_PSF_i','MAG_PSF_z','MAGERR_PSF_g','MAGERR_PSF_r',\
                  'MAGERR_PSF_i','MAGERR_PSF_z','MAG_MODEL_g','MAG_MODEL_r',\
                  'MAG_MODEL_i','MAG_MODEL_z','MAGERR_MODEL_g','MAGERR_MODEL_r',\
                  'MAGERR_MODEL_i','MAGERR_MODEL_z','SPREAD_MODEL_g','SPREAD_MODEL_r',\
                  'SPREAD_MODEL_i','SPREAD_MODEL_z',]]

    print 'number of stars =', len(total3)

    if mode2mass=='':
        starpsfmode = '_psf'
    elif mode2mass=='_no2mass':
        starpsfmode = '_no2mass'

    # total3.write(slrdir+'/star_psf%s_%s_%i.fits' % ('_psf',field,0), overwrite=True) #with 2MASS stars: star_psf_psf_%s_%i.fits
    total3.write(slrdir + '/star_psf%s_%s_%i.fits' % (starpsfmode, field, 0),
                 overwrite=True)  # no 2MASS star mode vs , '_psf' vs '_no2mass'


    if slr:
        slr_running_psf(field, infile=slrdir + '/star_psf%s_%s_%i.fits' %
                        (starpsfmode, field, 0), mode='psf', mode2mass=mode2mass)  # '_psf' vs '_no2mass'

    total_gal=Table.from_pandas(pd.read_csv("/Users/taweewat/Documents/pisco_code/slr_output/galaxy_psf_total_%s.csv"%(field)))
    print 'mode=', mode
    ntotal_gal = update_color(slrdir + '/star_psf%s_%s_%i.fits.offsets.list' %
                              (starpsfmode, field, 0), total_gal, mode=mode)  # '' vs '_no2mass
    ntotal_gal.write(os.path.join(
        slrdir, 'galaxy_psf%s_ntotal_%s.csv' % (mode2mass, field)), overwrite=True)  # '' vs '_no2mass'

def make_images(field,ax=None):
    dir='/Users/taweewat/Documents/pisco_code/Chips_images/'
    ax.imshow(image.imread(dir+"aplpy4_%s_img.jpeg"%field))
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axis('off')
    return None

def sur_pro(r): #Mpc
    def fn(x):
        if x>=1:
            return 1.-(2/np.sqrt(x**2-1)*np.arctan(np.sqrt((x-1.)/(x+1.))))
        elif x<1:
            return 1.-(2/np.sqrt(1-x**2)*np.arctanh(np.sqrt((1.-x)/(x+1.))))
    rs=0.15/0.71 #Mpc
    if r>=(0.1/0.71):
        return 1/((r/rs)**2-1)*fn(r/rs)
    elif r<(0.1/0.71):
        return 1./(((0.1/0.71)/rs)**2-1)*fn((0.1/0.71)/rs)
def k_NFW():
    def integrated(y):
        return 1./integrate.quad(lambda r: 2*np.pi*r*sur_pro(r),0,y)[0]
    xy=np.logspace(-3,3,num=30)
    X = np.log(xy)
    Y = np.log([integrated(np.e**(y)) for y in X])
    Z=np.polyfit(X,Y,6)
    k_NFW = np.poly1d(Z)
    return k_NFW
def sur_pro_prob(r,rc,k_NFW): #(Mpc,Mpc)
    # Weighted based on the distance from the center (Rykoff+12)
    # rc: the cutoff radius (truncated at the cluster radius)
    return np.e**(k_NFW(np.log(rc)))*sur_pro(r)

name=['z','dist','age','mass','Abs_g','App_g','kcorr_g','Abs_r',\
      'App_r','kcorr_r','Abs_i','App_i','kcorr_i','Abs_z','App_z','kcorr_z']
df=pd.read_csv('/Users/taweewat/Documents/red_sequence/rsz/model/'+\
               'ezmodel2_bc03_zf1.5_chab_0.016_exp_0.1.txt',\
               skiprows=27,delim_whitespace=True,names=name)
df=df[(df.z>0.1) & (df.z<1.)]
z_new=np.arange(0.1, 0.95, 0.0025)
Appi_new = interpolate.splev(z_new, interpolate.splrep(df.z, df.App_i, s=0), der=0)

def linear_rmi(x0,redshift):
    x=df.z[:-12]
    y=(df.App_r-df.App_i)[:-12]
    yhat = np.polyfit(x, y, 5)
    f_rmi = np.poly1d(yhat)
    slope=-0.0222174237562*1.007
    Appi0=Appi_new[np.where(abs(z_new-redshift)<=1e-9)[0][0]]
    return slope*(x0-Appi0)+f_rmi(redshift)

def linear_gmr(x0,redshift):
    x=df.z[:-25]
    y=(df.App_g-df.App_r)[:-25]
    yhat = np.polyfit(x, y, 5)
    f_gmr = np.poly1d(yhat)
    slope=-0.0133824600874*1.646
    Appi0=Appi_new[np.where(abs(z_new-redshift)<=1e-9)[0][0]]
    return slope*(x0-Appi0)+f_gmr(redshift)

def find_fits_dir(field):
    home = '/Users/taweewat/Documents/pisco_code/'
    dirs = ['ut170103/', 'ut170104/', 'ut170619/', 'ut170621/',\
            'ut170624/', 'ut171208/', 'ut171209/', 'ut171212/']
    myReg = re.compile(r'(%s_A).*' % field)
    for di in dirs:
        diri = home + di
        for text in os.listdir(diri):
            if myReg.search(text) != None:
                filename = myReg.search(text).group()
                allfilename = diri
    return allfilename

def pisco_tilt_resequence(field, mode='psf', mode2mass=''):
    if field[0:5]=='Field':
        base = pd.read_csv('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs.csv')
        RA = base[base.name == field].ra.values[0]
        DEC = base[base.name == field].dec.values[0]
        redshift = base[base.name == field].redshift.values[0]
    elif field[0:5]=='CHIPS':
        base = pd.read_csv('/Users/taweewat/Documents/red_sequence/chips_all_obj.csv', index_col=0)
        RA = base[base.chips == field].ra.values[0]
        DEC = base[base.chips == field].dec.values[0]
        redshift = base[base.chips == field].redshift.values[0]
    elif field[0:4]=='SDSS':
        base = pd.read_csv('/Users/taweewat/Documents/xray_project/ned-result/final_sdss_cut5.csv', index_col=0)
        RA = base[base.name == field].RA.values[0]
        DEC = base[base.name == field].DEC.values[0]
        redshift = base[base.name == field].redshift.values[0]
    elif field=='PKS1353':
        RA = 209.0225
        DEC = -34.3530556
        redshift = 0.223
    elif field=='CHIPS2249-2808':
        RA = 336.99975202151825
        DEC = -43.57623068466675
        redshift = -1
    elif field=='CHIPS2246-2854':
        RA = 335.7855174238757
        DEC = -34.934569299688185
        redshift = -1

    if redshift!=-1:
        qso_redshift=redshift
    else:
        qso_redshift=0.2
    
    ebv = ebvpy.calc_ebv(ra=[RA],dec=[DEC]); print 'ebv:', ebv[0]
    ebv_g=ebvpy.calc_color_correction('g', ebv)[0]; print 'ebv_g:', ebv_g
    ebv_r=ebvpy.calc_color_correction('r', ebv)[0]; print 'ebv_r:', ebv_r
    ebv_i=ebvpy.calc_color_correction('i', ebv)[0]; print 'ebv_i:', ebv_i

    param_izp=read_param_izp() #i zero point
    dir_dict = dict(zip(['ut170103/','ut170104/','ut170619/',\
    'ut170621/','ut170624/','ut171208/','ut171209/','ut171212/'], np.arange(1, 9)))

    # fname = "/Users/taweewat/Documents/pisco_code/slr_output/galaxy_ntotal_%s.csv"%field
    fname = "/Users/taweewat/Documents/pisco_code/slr_output/galaxy_psf%s_ntotal_%s.csv" % (
        mode2mass, field)  # '' vs '_no2mass'
    df0 = pd.read_csv(fname,index_col=0)

    c5 = SkyCoord(ra=df0['ALPHA_J2000_i'].values*u.degree, dec=df0['DELTA_J2000_i'].values*u.degree)
    c0 = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    sep = c5.separation(c0)
    df0['sep(deg)']=sep
    df0['sep(Mpc)']=sep*60.*cosmo.kpc_proper_per_arcmin(qso_redshift).value/1e3
    cut=df0

    dfi=cut
    # Added Galactic Reddening (6/16/18)
    if mode2mass == '':
        dfi['MAG_i']=dfi['MAG_i']-ebv_i
        dfi['MAG_g']=dfi['MAG_g']-ebv_g
        dfi['MAG_r']=dfi['MAG_r']-ebv_r
    # Use i Zero Point from each day and g,r zero point fron the color (6/22/18)
    elif mode2mass == '_no2mass':
        dfi['MAG_i']=dfi['MAG_i']-ebv_i+param_izp['i_zp_day%i'%dir_dict[find_fits_dir(field)[-9:]]]
        dfi['MAG_g']=dfi['MAG_g']-ebv_g+param_izp['i_zp_day%i'%dir_dict[find_fits_dir(field)[-9:]]]
        dfi['MAG_r'] = dfi['MAG_r'] - ebv_r + param_izp['i_zp_day%i' %
                                                        dir_dict[find_fits_dir(field)[-9:]]]
        dfi['MAGERR_i']=0

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
    dfi.loc[:,"w_gmr"] = np.nan# [sur_pro_prob(r,1.) for r in dfi['sep(Mpc)']]
    dfi.loc[:,"w_rmi"] = np.nan# [sur_pro_prob(r,1.) for r in dfi['sep(Mpc)']]

    k_NFW0=k_NFW()

    z_rmi,w_rmi=[],[]
    for i, row in dfi.iterrows():
        for z in np.arange(0.35,0.7,0.025):
            if row['MAG_i'] < -18+5.*np.log10(ex.d_L(z)*1e6)-5.:
             # if row['MAG_i'] < magi_cut_rmi:
                rmi=row['MAG_r']-row['MAG_i']
                low_edge=linear_rmi(row['MAG_i'],round(z-0.0125,4))
                high_edge=linear_rmi(row['MAG_i'],round(z+0.0125,4))
                if (rmi > low_edge) & (rmi < high_edge):
                    # if (np.sqrt(row['MAGERR_r']**2+row['MAGERR_i']**2) < 3.5*(high_edge-low_edge)):
                    z_rmi.append(round(z,3))
                    w_rmi.append(sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0))
                    dfi.loc[i,"z_rmi"]=z
                    dfi.loc[i,"w_rmi"]=sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0)

    z_gmr,w_gmr=[],[]
    for i, row in dfi.iterrows():
        for z in np.arange(0.15,0.35,0.025):
            if row['MAG_i'] < -18+5.*np.log10(ex.d_L(z)*1e6)-5.:
            # if row['MAG_i'] < magi_cut_gmr:
                gmr=row['MAG_g']-row['MAG_r']
                low_edge=linear_gmr(row['MAG_i'],round(z-0.0125,4))
                high_edge=linear_gmr(row['MAG_i'],round(z+0.0125,4))
                if (gmr > low_edge) & (gmr < high_edge):
                    # if (np.sqrt(row['MAGERR_g']**2+row['MAGERR_r']**2) < 3.5*(high_edge-low_edge)):
                    z_gmr.append(round(z,3))
                    w_gmr.append(sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0))
                    dfi.loc[i,"z_gmr"]=z
                    dfi.loc[i,"w_gmr"]=sur_pro_prob(row['sep(Mpc)'],1.,k_NFW0)

    ns1,xs1=np.histogram(z_gmr,bins=np.arange(0.125,0.35,0.025),weights=w_gmr)
    ns2,xs2=np.histogram(z_rmi,bins=np.arange(0.325,0.7,0.025),weights=w_rmi)
    z_total=np.append(xs1[:-1],xs2[:-1])
    n_total=np.append(ns1,ns2)
    z_max=z_total[np.where(n_total==np.max(n_total))[0][0]]
    # n_max=np.max(n_total)

    hist2, bin_edges2 = np.histogram(z_gmr,bins=np.arange(0.125,0.35,0.025),weights=w_gmr)
    hist, bin_edges = np.histogram(z_rmi,bins=np.arange(0.325,0.7,0.025),weights=w_rmi)

    cmap=matplotlib.cm.RdYlGn

    fig,ax=plt.subplots(1,4,figsize=(20,5))
    make_images(field,ax[0])

    dfi_ri=dfi.loc[dfi['z_rmi'].dropna().index]
    ax[1].scatter(dfi['MAG_i'],dfi['MAG_r']-dfi['MAG_i'],c='black',alpha=0.1)#dfi['w_rmi'],cmap=cmap)
    ax[1].scatter(dfi_ri['MAG_i'],dfi_ri['MAG_r']-dfi_ri['MAG_i'],c=dfi_ri['w_rmi'],cmap=cmap)
    ax[1].errorbar(dfi_ri['MAG_i'],dfi_ri['MAG_r']-dfi_ri['MAG_i'],xerr=dfi_ri['MAGERR_i'],yerr=np.sqrt(dfi_ri['MAGERR_r']**2+dfi_ri['MAGERR_i']**2),fmt='none',c='k',alpha=0.05)
    # plt.plot(df.App_i,df.App_r-df.App_i,'.')
    # ax[1].axhline(xs[:-1][(xs[:-1]<1.33) & (xs[:-1]>0.6)][0],lw=0.7,color='green')
    for z in np.arange(0.35,0.7,0.025):
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
    ax[2].scatter(dfi_gr['MAG_i'],dfi_gr['MAG_g']-dfi_gr['MAG_r'],c=dfi_gr['w_gmr'],cmap=cmap)
    ax[2].errorbar(dfi_gr['MAG_i'],dfi_gr['MAG_g']-dfi_gr['MAG_r'],xerr=dfi_gr['MAGERR_i'],yerr=np.sqrt(dfi_gr['MAGERR_g']**2+dfi_gr['MAGERR_r']**2),fmt='none',c='k',alpha=0.05)
    # plt.plot(df.App_i,df.App_g-df.App_r,'.')
    # ax[2].axhline(xs[:-1][(xs[:-1]<1.65) & (xs[:-1]>np.min(x2))][0],lw=0.7,color='green')
    for z in np.arange(0.15,0.35,0.025):
        ax[2].plot(I,linear_gmr(I,round(z,4)),color=s_m.to_rgba(z))
    ax[2].set_ylim(0.75,2)
    ax[2].set_xlim(16,24)
    # cbar=plt.colorbar(s_m)
    ax[2].set_xlabel('I')
    ax[2].set_ylabel('G-R')
    ax[2].set_title('z=0.15-0.325')#, icut:'+str(magi_cut_gmr))
    # plt.plot([corr_f(z) for z in df.z.values[:-25]],df.App_g[:-25]-df.App_r[:-25],'-')

    ax[3].bar(bin_edges[:-1], hist, width = 0.025, color='#1f77b4')#, alpha=0.5)
    ax[3].bar(bin_edges2[:-1], hist2, width = 0.025, color='#ff7f0e')#, alpha=0.5)
    ax[3].axvline(z_max,ls='--',color='purple',label='z_max='+str(z_max))
    ax[3].axvline(redshift,color='red',label='z: '+str(round(redshift,3)))
    ax[3].legend(loc='best')
    ax[3].set_xlabel('z')
    ax[3].set_xlim(0.1,0.7)
    ax[3].set_title('ebv:%.3f,ebv_g-r:-%.3f,ebv_r-i:-%.3f'%(ebv[0],ebv_g-ebv_r,ebv_r-ebv_i))

    if np.max(n_total)<30:
        ax[3].set_ylim(0,30)

    plt.tight_layout(rect=[0, 0., 1, 0.98])
    plt.savefig('/Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq%s_%s_all_%.2f_%s_tilted.png' % ('',mode,np.max(n_total),field), dpi=120)
    plt.close(fig)

def pisco_combine_imgs(fields, mode='psf', mode2mass=''):
    dir1='/Users/taweewat/Documents/red_sequence/pisco_color_plots/psf_est/'
    dir2='/Users/taweewat/Documents/red_sequence/pisco_color_plots/'
    dir3='/Users/taweewat/Documents/red_sequence/pisco_color_plots/'
    dirout='/Users/taweewat/Documents/red_sequence/pisco_all/'

    myReg = re.compile(r'(redsq_%s_all_.*%s.*png)' % (mode, field))
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
    h = imgs[0].height+imgs[1].height+imgs[2].height/2
    result = Image_PIL.new("RGBA", (mw, h))
    y,index=0,0
    for i in imgs:
        if index==2:
            i=i.resize((i.width/2,i.height/2))
        result.paste(i, (0, y))
        y += i.size[1]
        index+=1
    # result.save(dirout + 'all_combine%s_%s_%s_%s.png' %
    #             (mode2mass, field, mode, myReg2.search(names[0]).group())) 
    result.save(dirout + '%s_all_combine%s_%s_%s.png' %
                (myReg2.search(names[0]).group(), mode2mass, field, mode)) 

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            print 'remove', f
            os.remove(os.path.join(dir, f))

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

    # all_fields=['Field042'] 'Field047','CHIPS0150-3338'
    # all_fields=['CHIPS1034-2837']

    # all_fields=['Field084','Field085','Field086','Field087',\
    #  'Field082','Field083','CHIPS0146-3648','CHIPS0821-0815','CHIPS0127-4016','CHIPS0106-2358','CHIPS0137-1248']
    # all_fields=['Field143','CHIPS2245-4931','CHIPS0107-1310','CHIPS0300-3413','CHIPS2223-3455','Field046','CHIPS2251-3827']
    # all_fields=['Field228','CHIPS0310-3723','CHIPS0022-5140']

    # all_fields=['Field226','Field227','Field022','Field229'] #'Field029' (bad color)
    # all_fields=['CHIPS0022-5140']

    home='/Users/taweewat/Documents/pisco_code/' #09, 171208
    dirs=['ut170103/','ut170104/','ut170619/','ut170621/','ut170624/','ut171208/','ut171209/','ut171212/']
    names=[]
    myReg=re.compile(r'(CHIPS\d{4}[+-]\d{4})|(Field\d{3})')
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                names.append(myReg.search(text).group())
    all_fields=list(set(names))
    # print all_fields

    exception=[]
    #'CHIPS0411-5149','CHIPS1011-0505','CHIPS0222-4159','CHIPS0040-2902','CHIPS0025-5427',\
    #     'CHIPS2243-3034','CHIPS0012-1628','CHIPS0302-2758','CHIPS2251-3210','CHIPS0005-2758','CHIPS0132-1608',\
    #     'CHIPS2340-2302','CHIPS0304-3556','CHIPS0018-1840','CHIPS1141-1407','CHIPS1142-1422','CHIPS2311-4718',\
    #     'CHIPS0325-4926','Field029','Field042','CHIPS0150-3338','CHIPS0335-4715','CHIPS0118-1430','CHIPS1034-2837',\
    #     'Field237','Field234','Field037','CHIPS0140-1533','CHIPS0409-2839','Field228','CHIPS0022-5140',\
    #     'CHIPS0512-3257','CHIPS2218-2947','CHIPS0536-3401','CHIPS0152-5028','Field060','CHIPS0355-6645','CHIPS0514-5046',\
    #     'CHIPS0050-5249','CHIPS0153-3143','CHIPS0157-1043','CHIPS2227-4333','CHIPS0316-2247','CHIPS0724-0715','Field116',\
    #     'CHIPS0253-5441','CHIPS2228-3220','CHIPS2306-3439','CHIPS2307-4236','CHIPS1205-2633','Field105',\
    #     'Field103','CHIPS0004-4736','CHIPS0024-6820','CHIPS2303-6807','CHIPS0342-3703','CHIPS0827-2026',\
    #     'CHIPS2348-3831','CHIPS0525-6938','CHIPS0847-0703','CHIPS0449-4350','CHIPS0449-3910','CHIPS2325-4800','CHIPS2254-3635','CHIPS0219-3626']
    #'Field166', 'Field054', 'CHIPS0229-5232','CHIPS0552-2336'

    # all_fields = ['Field092', 'Field036', 'Field056', 'Field055', 'Field204', 'Field205', 'Field201',\
    #               'Field198', 'Field219', 'Field218', 'Field071', 'Field073', 'Field075', 'Field182', 'Field151']
    
    # all_fields = ['Field137']
    # all_fields=['Field101','Field027','Field024','Field102']
    # ['Field166','Field137','Field048','Field039']  #['Field166']
    # all_fields = ['CHIPS0745-0714']
    # all_fields = ['Field187']
    # ['CHIPS0003-2521', 'CHIPS2349-2352', 'CHIPS2340-2302','CHIPS0106-1149']
    # all_fields = ['Field039','Field233','Field084','Field210'] 
    # ['CHIPS0824-3020', 'CHIPS0936-3342','CHIPS1009-3015', 'CHIPS1036-3513']

    #Day8 2mass:  u'CHIPS0106-2358' (1 2mass stars)
    # all_fields=[u'CHIPS2317-1443',u'CHIPS2333-2407',u'CHIPS2349-2352',u'CHIPS2357-1125',u'CHIPS0003-2521',u'CHIPS0015-1513', u'CHIPS0050-1412',\
    # u'CHIPS0106-1149',u'CHIPS0107-1310',u'CHIPS0116-1136',u'CHIPS0122-2646',u'CHIPS0137-1248']
    # Day7 2mass
    # all_fields=[u'CHIPS2101-6132',u'CHIPS2127-4725',u'CHIPS2133-2815',u'CHIPS2139-3911',u'CHIPS2141-3729',u'CHIPS2148-2735',u'CHIPS2148-3715',u'CHIPS2210-5508',\
    # u'CHIPS2211-3707',u'CHIPS2216-2803',u'CHIPS2217-3034',u'CHIPS2221-2804',u'CHIPS2246-2854',u'CHIPS2249-2808',u'CHIPS0004-2902',u'CHIPS0112-2919',\
    # u'CHIPS0115-3047',u'CHIPS0206-7148',u'CHIPS0209-6810',u'CHIPS0300-3413',u'CHIPS0310-3723',u'CHIPS0148-2238',u'CHIPS0522-1745',u'CHIPS0303-2407',\
    # u'CHIPS0609-0247',u'CHIPS0745-0714',u'CHIPS0821-0815',u'CHIPS0849-1721',u'CHIPS0920-2257',u'CHIPS0934-1721',u'CHIPS1102-3031',u'CHIPS1147-1252']

    # Day6: 2mass CHIPS2251-3827 (1 2mass star)
    # all_fields=[u'CHIPS2333-6144',\
    # u'CHIPS0049-4457',u'CHIPS0127-4016',u'CHIPS0133-4358',u'CHIPS0146-3648',u'CHIPS0146-3711',u'CHIPS0423-3953',\
    # u'CHIPS0449-2859',u'CHIPS0532-3917',u'CHIPS0535-6602',u'CHIPS0824-3020',u'CHIPS0936-3342',u'CHIPS0957-7554',\
    # u'CHIPS1009-3015',u'CHIPS1036-3513']

    # Day3: 
    # all_fields = [u'Field074',u'Field045',u'Field058',u'Field059',u'Field210',u'Field212',u'Field047',u'Field225',u'Field226',u'Field046',u'Field072',u'Field076',u'Field084',u'Field085',u'Field048',u'Field038',u'Field087',
    #               u'Field233',u'Field025',u'Field020',u'Field039',u'Field018',
    #               u'Field021',
    #               u'Field022',
    #               u'Field077',
    #               u'Field121',
    #               u'Field266',
    #               u'Field269',
    #               u'Field115',
    #               u'Field088',
    #               u'Field274',
    #               u'Field279',
    #               u'Field292',
    #               u'Field124']

    # all_fields = ['CHIPS0118-1430', 'CHIPS0003-2521']
    # all_fields = ['Field074', 'Field044', 'Field136']
    # all_fields = ['SDSS123','SDSS501','SDSS603']
    notgoflag=True
    for index, field in enumerate(all_fields[:]):
        print field, '%i/%i'%(index,len(all_fields))

        if field in exception:
            continue

        if field == 'CHIPS0525-6938':
            notgoflag=False
            continue
        if notgoflag:
            continue 

        if slr=='allslr':
            star_galaxy_bleem(field)
            pisco_photometry_v4(field)
        pisco_cut_frame(field)
        pisco_photometry_psf_v4(field, mode=mode, mode2mass=mode2mass, slr=slr_param)
        purge('/Users/taweewat/Documents/red_sequence/pisco_color_plots/'\
              ,r'(redsq_%s_all_.*%s.*png)'%(mode,field))
        pisco_tilt_resequence(field, mode=mode, mode2mass=mode2mass)
        pisco_combine_imgs(field, mode=mode, mode2mass=mode2mass)
    
        purge('final', "proj_coadd_c%s_.*\.fits" % field)
        purge('.', "proto_psf_%s_.*\.fits" % field)
        purge('.', "samp_psf_%s_.*\.fits" % field)
        purge('.', "resi_psf_%s_.*\.fits" % field)
        purge('.', "snap_psf_%s_.*\.fits" % field)
        purge('.', "chi_psf_%s_.*\.fits" % field)
        purge('psfex_output', "psf_%s_.*\.fits" % field)
        purge('slr_output', "a_psf_%s_.*\.fits" % field)
        purge('final', "coadd_c%s_sq_.*\.fits" % field)
        