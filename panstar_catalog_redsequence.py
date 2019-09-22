## python pisco_pipeline/panstar_catalog_redsequence.py

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image
import matplotlib
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import interp1d
import scipy.integrate as integrate

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=71, Om0=0.3, Tcmb0=2.725)

import extra_program as ex
##----------

name=['z','dist','age','mass','Abs_g','App_g','kcorr_g','Abs_r',\
      'App_r','kcorr_r','Abs_i','App_i','kcorr_i','Abs_z','App_z','kcorr_z']
df=pd.read_csv('/Users/taweewat/Documents/red_sequence/rsz/model/'+\
                'ezmodel2_bc03_zf2.5_chab_0.02_exp_0.1.txt',
            #    'ezmodel2_c09_zf3.0_chab_0.02_exp_0.1.txt',
               skiprows=27,delim_whitespace=True,names=name)
df=df[(df.z>=0.1) & (df.z<1.)]
z_new=np.arange(0.1, 0.95, 0.0025)
# Appi_new = interpolate.splev(z_new, interpolate.splrep(df.z, df.App_i, s=0), der=0)
Appi_f = interpolate.interp1d(df.z, df.App_i, kind='cubic')
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
##---------


norm = matplotlib.colors.Normalize(vmin=0.15,vmax=0.675)
c_m = matplotlib.cm.cool
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

bins_gmr_cen = np.arange(0.15815, 0.33315+0.01, 0.035)
bins_rmi_cen = np.arange(0.36815, 0.64815+0.01, 0.035) 
I=np.arange(16,24,0.01)


## PARAMETER ##
core_radius=0.25
extra_name='panstar_cat2_sur0.25'

# field='CHIPS1011-0505'

# total=pd.read_csv('/Users/taweewat/Documents/red_sequence/chips_all_obj.csv',index_col=0)
total=pd.read_csv('/Users/taweewat/Documents/xray_project/ned-result/total_cut_des_panstar.csv',index_col=0)
panstar=total[total['panstar']==True]
panstar_cut=panstar#[panstar.chips=='CHIPS1011-0505'].copy()
panstar_cut2=panstar_cut['chips'][0:1]

infile = open('/Users/taweewat/Documents/xray_project/code_github/allremove_chips.txt', 'r')
exception = [i.strip() for i in infile.readlines()]

for i, field in enumerate(panstar_cut2):

    if field in exception:
        continue
        
    print "%i/%i: %s"%(i, len(panstar_cut2), field)
    # df0=pd.read_csv('/Users/taweewat/Documents/red_sequence/panstarr_csv/pnew3_{}.csv'.format(field),index_col=0)
    df0=pd.read_csv('/Users/taweewat/Documents/red_sequence/panstarr_csv/pnew_dr2_{}.csv'.format(field),index_col=0)
    if df0.shape[0]==0:
        continue
    df0=df0.drop_duplicates(subset=['objID'], keep="first").copy()
    df0=df0[df0['qualityFlag']<100].copy()
    df0=df0[~((df0['iMeanKronMag']==-999) | (df0['rMeanKronMag']==-999))].copy() #had i and r band
    df0=df0[~((df0['gMeanKronMag']==-999) & (df0['rMeanKronMag']==-999) & (df0['iMeanKronMag']==-999))].copy() #had magnitude

    base = pd.read_csv('/Users/taweewat/Documents/red_sequence/chips_all_obj.csv', index_col=0)
    RA = base[base.chips == field].ra.values[0]
    DEC = base[base.chips == field].dec.values[0]
    redshift = base[base.chips == field].redshift.values[0]

    # dfi=df_all_gal_cut.copy()
    # dfi=df0.copy()
    df0['MAG_g']=df0['gMeanKronMag']
    df0['MAG_r']=df0['rMeanKronMag']
    df0['MAG_i']=df0['iMeanKronMag']
    df0['MAGERR_g']=df0['gMeanKronMagErr']
    df0['MAGERR_r']=df0['rMeanKronMagErr']
    df0['MAGERR_i']=df0['iMeanKronMagErr']

    c5 = SkyCoord(ra=df0['raMean'].values*u.degree, dec=df0['decMean'].values*u.degree)
    c0 = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    sep = c5.separation(c0)
    df0['sep(deg)']=sep
    df0=df0[df0['sep(deg)']*60.<=10].copy()

    df0['STAR_CUT']=np.array(((df0['iMeanPSFMag']-df0['iMeanKronMag'])<0.05)&\
                            ((df0['iMeanPSFMag']-df0['iMeanKronMag'])<100.), dtype=int)\
        +np.array(((df0['rMeanPSFMag']-df0['rMeanKronMag'])<0.05)&\
                            ((df0['rMeanPSFMag']-df0['rMeanKronMag'])<100.), dtype=int)\
        +np.array(((df0['gMeanPSFMag']-df0['gMeanKronMag'])<0.05)&\
                            ((df0['gMeanPSFMag']-df0['gMeanKronMag'])<100.), dtype=int)

    df0['GALAXY_CUT']=np.array(((df0['iMeanPSFMag']-df0['iMeanKronMag'])>0.05)&\
                            ((df0['iMeanPSFMag']-df0['iMeanKronMag'])<100.), dtype=int)\
        +np.array(((df0['rMeanPSFMag']-df0['rMeanKronMag'])>0.05)&\
                            ((df0['rMeanPSFMag']-df0['rMeanKronMag'])<100.), dtype=int)\
        +np.array(((df0['gMeanPSFMag']-df0['gMeanKronMag'])>0.05)&\
                            ((df0['gMeanPSFMag']-df0['gMeanKronMag'])<100.), dtype=int)

    df_all_gal_cut=df0[df0['GALAXY_CUT']>0].copy()
    df_all_star=df0[df0['STAR_CUT']==3].copy()
    print df0.shape, df_all_gal_cut.shape, df_all_star.shape

    dfi=df_all_gal_cut.copy()

    dfi2=dfi[np.sqrt(dfi['MAGERR_r']**2+dfi['MAGERR_i']**2)<0.5].copy()
    x=dfi2[dfi2['MAG_i']>13]['MAG_i']
    y=np.sqrt(dfi2[dfi2['MAG_i']>13]['MAGERR_r']**2+dfi2[dfi2['MAG_i']>13]['MAGERR_i']**2)
    p=np.poly1d(np.polyfit(x,np.log(y+1e-9),1, w=np.sqrt(y)))
    Mag_cut=(p-np.log(0.067*3)).roots;
    print "Mag_cut: %.2f"%(Mag_cut)
    xs=np.arange(np.min(x),np.max(x),0.01)

    # fig,ax=plt.subplots(figsize=(5,5))
    # plt.plot(x,y,'.',label='r-i')
    # plt.plot(xs,np.exp(p(xs)),label='exp({:.1f}+{:.1f}x)'.format(p[0],p[1]))
    # plt.xlabel('Mag_i'); plt.ylabel('$\Delta r-i$ err')
    # plt.ylim(-0.05,0.35)
    # plt.axvline(Mag_cut,label='Mag_cut')
    # plt.legend(loc='best')
    # plt.close()

    I=np.arange(16,24,0.01)
    dfi.loc[:,"z_gmr"] = np.nan
    dfi.loc[:,"z_rmi"] = np.nan
    dfi.loc[:,"w_gmr"] = np.nan
    dfi.loc[:,"w_rmi"] = np.nan
    dfi.loc[:,"w_col_gmr"] = np.nan
    dfi.loc[:,"w_col_rmi"] = np.nan

    bin_width=0.035 #0.025

    # bins_gmr_cen = np.arange(0.15815, 0.33315+0.01, bin_width) # bins_gmr_cen = np.arange(0.15, 0.325+0.01, bin_width)
    bins_gmr_cen = np.arange(0.12315, 0.33315+0.01, bin_width)
    bins_rmi_cen = np.arange(0.36815, 0.64815+0.01, bin_width) # bins_rmi_cen = np.arange(0.36, 0.675+0.01, bin_width)
    # bins_gmr_edge = np.arange(0.14065, 0.35065+0.01, bin_width) # bins_gmr_edge = np.arange(0.1325, 0.3425+0.01, bin_width)
    bins_gmr_edge = np.arange(0.10565, 0.35065 + 0.01, bin_width)
    bins_rmi_edge = np.arange(0.35065, 0.66565+0.01, bin_width) # bins_rmi_edge = np.arange(0.3425, 0.6925+0.01, bin_width)

    z_rmi,w_rmi,w_col_rmi=[],[],[]
    for i, row in dfi.iterrows():
        for z in bins_rmi_cen:
            if row['MAG_i'] < Mag_cut:
                rmi=row['MAG_r']-row['MAG_i']
                low_edge=linear_rmi(row['MAG_i'],round(z-0.0175,4)) 
                high_edge=linear_rmi(row['MAG_i'],round(z+0.0175,4))
                if (rmi > low_edge) & (rmi <= high_edge):
    #                 if row['sep(deg)']*60.<4.: #less than 4 arcmin
                    z_rmi.append(round(z,3))
                    wrmi0=ex.sur_pro_prob_ang(row['sep(deg)']*60, core_radius); w_rmi.append(wrmi0) #arcmin
                    w_col_rmi0=1.; w_col_rmi.append(w_col_rmi0)
                    dfi.loc[i,"z_rmi"]=z
                    dfi.loc[i,"w_rmi"]=wrmi0
                    dfi.loc[i,"w_col_rmi"]=w_col_rmi0

    z_gmr,w_gmr,w_col_gmr=[],[],[]
    for i, row in dfi.iterrows():
        for z in bins_gmr_cen:
            if row['MAG_i'] < Mag_cut:
                gmr=row['MAG_g']-row['MAG_r']
                low_edge=linear_gmr(row['MAG_i'],round(z-0.0175,4))
                high_edge=linear_gmr(row['MAG_i'],round(z+0.0175,4))
                if (gmr > low_edge) & (gmr <= high_edge):
    #                 if row['sep(deg)']*60.<4.:#less than 4 arcmin
                    z_gmr.append(round(z,3))
                    w_col_gmr0=1.; w_col_gmr.append(w_col_gmr0)
                    wgmr0 = ex.sur_pro_prob_ang(row['sep(deg)'] * 60, core_radius); w_gmr.append(wgmr0)  # arcmin
                    dfi.loc[i,"z_gmr"]=z
                    dfi.loc[i,"w_gmr"]=wgmr0
                    dfi.loc[i,"w_col_gmr"]=w_col_gmr0

    ns1,xs1=np.histogram(z_gmr,bins=bins_gmr_edge,weights=np.array(w_gmr)*np.array(w_col_gmr)) #0.15-0.325
    bin_cen1 = (xs1[:-1] + xs1[1:])/2
    ns2,xs2=np.histogram(z_rmi,bins=bins_rmi_edge,weights=np.array(w_rmi)*np.array(w_col_rmi)) #0.36-0.675
    bin_cen2 = (xs2[:-1] + xs2[1:])/2
    z_total=np.append(bin_cen1, bin_cen2)
    n_total=np.append(ns1,ns2)
    z_max=z_total[np.where(n_total==np.max(n_total))[0][0]]
    n_max=np.max(n_total)
    n_median = np.median(n_total[n_total != 0])
    n_mean = np.mean(n_total)

    n_bkg = np.mean(sorted(n_total)[2:-2]);
    z_total_added = np.insert(
        np.append(z_total, z_total[-1] + bin_width), 0, z_total[0] - bin_width)
    n_total_added = np.insert(np.append(n_total, 0), 0, 0) - n_bkg
    indi = np.where(n_total_added == np.max(n_total_added))[0][0]
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

    lumfn=pd.read_csv('/Users/taweewat/Documents/red_sequence/coma_cluster_luminosity_function/schecter_fn.csv',\
    names=['M_r','theta(M)Mpc^-3']) 
    h=0.7; 
    x=lumfn['M_r']+5*np.log10(h);
    y=lumfn['theta(M)Mpc^-3']*(h**3);
    f1d=interp1d(x, y,kind='cubic')
    z_max_fit = tuple(popt)[1]
    lum_factor = integrate.quad( f1d, -23.455, ex.abs_mag(22.25, z_max_fit))[0]  #-23.455: min abs Mag from schecter_fn.csv, 22.25: median of Mag r
    density_factor=integrate.quad(ex.NFW_profile, 0.001, core_radius*cosmo.kpc_proper_per_arcmin(z_max_fit).value/1e3)[0]
    signal = tuple(popt)[0] / (lum_factor * density_factor)
    

    print 'z_max_fit', z_max_fit

    dfi_ri=dfi.loc[dfi['z_rmi'].dropna().index]
    cmap=matplotlib.cm.RdYlGn
    norm = matplotlib.colors.Normalize(vmin=0.15,vmax=0.675)
    c_m = matplotlib.cm.cool
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    plt.figure(figsize=(18,4.5))

    plt.subplot(1,4,1)
    dir='/Users/taweewat/Documents/pisco_code/Chips_images/'
    try:
        plt.imshow(image.imread(dir+"aplpy_panstar_%s_img4.jpeg"%field))
    except:
        plt.imshow(image.imread(dir+"aplpy_panstar_%s_img.jpeg"%field))
    # ax.axes.get_xaxis().set_visible(False)
    # ax.axes.get_yaxis().set_visible(False)
    # ax.axis('off')

    plt.subplot(1,4,2)
    for z in bins_rmi_cen:
        plt.plot(I,linear_rmi(I,round(z,4)),color=s_m.to_rgba(z))
    plt.scatter(dfi['MAG_i'],dfi['MAG_r']-dfi['MAG_i'],c='black',alpha=0.05)#dfi['w_rmi'],cmap=cmap)
    plt.scatter(dfi_ri['MAG_i'],dfi_ri['MAG_r']-dfi_ri['MAG_i'],c=dfi_ri['w_rmi'],cmap=cmap)#,norm=norm)
    plt.errorbar(dfi_ri['MAG_i'],dfi_ri['MAG_r']-dfi_ri['MAG_i'],xerr=dfi_ri['MAGERR_i'],yerr=np.sqrt(dfi_ri['MAGERR_r']**2+dfi_ri['MAGERR_i']**2),fmt='none',c='k',alpha=0.05)
    plt.ylim(0.25,1.5)
    plt.xlim(16,24)
    plt.ylabel('R-I')
    plt.xlabel('I')

    dfi_gr=dfi.loc[dfi['z_gmr'].dropna().index]
    plt.subplot(1,4,3)
    for z in bins_gmr_cen:
        plt.plot(I,linear_gmr(I,round(z,4)),color=s_m.to_rgba(z))
    plt.scatter(dfi['MAG_i'],dfi['MAG_g']-dfi['MAG_r'],c='black',alpha=0.05)#,c=dfi['w_gmr'],cmap=cmap)
    plt.scatter(dfi_gr['MAG_i'],dfi_gr['MAG_g']-dfi_gr['MAG_r'],c=dfi_gr['w_gmr'],cmap=cmap)#,norm=norm)
    plt.errorbar(dfi_gr['MAG_i'],dfi_gr['MAG_g']-dfi_gr['MAG_r'],xerr=dfi_gr['MAGERR_i'],yerr=np.sqrt(dfi_gr['MAGERR_g']**2+dfi_gr['MAGERR_r']**2),fmt='none',c='k',alpha=0.05)
    plt.ylim(0.75,2)
    plt.xlim(16,24)
    plt.ylabel('G-R')
    plt.xlabel('I')

    plt.subplot(1,4,4)
    xs=np.arange(np.min(z_fit)-0.1,np.max(z_fit)+0.1,0.001)
    plt.bar(bin_cen2, ns2, width=bin_width, color='#1f77b4')
    plt.bar(bin_cen1, ns1, width=bin_width, color='#ff7f0e')
    plt.plot(z_fit,n_fit+n_bkg,'o',c='tab:purple')
    plt.axvline(z_max,ls='--',color='purple',label='z_max:%.2f'%z_max)
    plt.axvline(redshift,color='red',label='z:%.2f'%(redshift))
    plt.axhline(n_median,color='tab:green',label='median:%.2f'%n_median)
    plt.axhline(n_mean,color='tab:red',label='mean:%.2f'%n_mean)
    plt.plot(xs, gaussian_func(xs, *popt)+n_bkg, c='tab:green', ls='--', label='fit: a=%.2f, mu=%.4f'% tuple(popt))
    plt.legend(loc='best')
    plt.xlabel('z')
    plt.xlim(0.1,0.7)
    plt.ylim(0,30)

    plt.suptitle(field)
    plt.tight_layout()

    plt.savefig('/Users/taweewat/Documents/red_sequence/panstarr_color_plots/%s_%.2f_%s.png' %(extra_name,signal,field), dpi=120)
    plt.close()
        
    red_dir='/Users/taweewat/Documents/red_sequence/'
    rich_filename = 'all_richness_%s.csv'%extra_name
    if not os.path.isfile(red_dir + rich_filename):
        df_richness=pd.DataFrame(columns=['name','zmax','nmax','signal','zmax_fit','lum_factor','density_factor'])
        df_richness.to_csv(red_dir+rich_filename)
    
    df_richness=pd.read_csv(red_dir+rich_filename,index_col=0)
    dic={'name':field, 'zmax':z_max,'nmax':n_max,'signal':signal,'zmax_fit':z_max_fit,'lum_factor':lum_factor,'density_factor':density_factor}
    if field in df_richness['name'].values: 
        df_richness=df_richness[df_richness['name']!=field]
    df_richness=df_richness.append(pd.Series(dic),ignore_index=True).copy()
    df_richness.to_csv(red_dir+rich_filename)