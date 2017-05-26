import sys
import os
import pandas as pd
import numpy as np
import subprocess
import shlex

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib import image
import matplotlib

import extra_program as ex

import ezgal
from rsz import RSModel

##----
def make_images(field,ax=None):
    dir='final/'
    ax.imshow(image.imread(dir+"img%s_2.eps" % field))
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    return None

def red_seq_color_plot(color,df,mags,ax=None):
    if ax is None:
        ax = plt.gca()
    #https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Jordi2006
    n=1000 #repeat number for sampling

    slope_fit,i_band0,color_err,rs_models,band_1,band_2=color_sloan(color, mags)

    ysample=df[df_color[band_1]]-df[df_color[band_2]]
    ysample_err=np.sqrt(df[df_colorerr[band_1]]**2+df[df_colorerr[band_2]]**2)
    total=[]
    for i in ysample.index:
        total.append(np.random.normal(loc=ysample[i],scale=ysample_err[i],size=n))
        #total.append(0.1)
    total=np.array(total)

    band_x='sloan_i'
    all_x=np.repeat(df[df_color[band_x]],n)
    total=np.reshape(total, len(all_x))
    bp=ax.errorbar(df[df_color[band_x]],ysample,yerr=ysample_err,fmt='.',alpha=0.5)
    #bp=ax.errorbar(df[df_color[band_x]],ysample,fmt='.',alpha=0.5)

    red_band=np.arange(16,25,0.01)  #just for the line plot in the 3rd plot
    redshift_range=np.arange(0.10,0.8,0.05) #for the actual data
    number=[]
    if color=='sloan_g-sloan_r':
        redshift_range=np.arange(0.10,0.36,0.05)
    elif color=='sloan_r-sloan_i':
        redshift_range=np.arange(0.10,0.71,0.05)


    for redshift in redshift_range:
        if color=='sloan_g-sloan_r':
#            i_band_cut=20.5
            i_band_cut=i_band0+5.*np.log10(ex.d_L(redshift)*1e6)-5.
        elif color=='sloan_r-sloan_i':
            i_band_cut=i_band0+5.*np.log10(ex.d_L(redshift)*1e6)-5.
        aa=red_band<i_band_cut

        loc=[(all_x<i_band_cut)&\
             (total < rs_models[color][round(redshift+0.025,2)].rs_color(all_x))&\
             (total > rs_models[color][round(redshift-0.025,2)].rs_color(all_x))][0]
        number.append(np.sum(loc))
        ax.plot(red_band[aa],rs_models[color][round(redshift,2)].rs_color(red_band[aa]),\
                 color=s_m.to_rgba(round(redshift,2)),ls='-')
        ax.plot(red_band[aa],rs_models[color][round(redshift+0.025,2)].rs_color(red_band[aa]),\
                 color=s_m.to_rgba(round(redshift,2)),ls=':')
        ax.plot(red_band[aa],rs_models[color][round(redshift-0.025,2)].rs_color(red_band[aa]),\
                 color=s_m.to_rgba(round(redshift,2)),ls=':')

    ax.set_xlim(16,25)
    if color == 'sloan_g-sloan_i':
        ax.set_ylim(0,4)
    elif color == 'sloan_g-sloan_r':
        ax.set_ylim(0.0,2.5)
    else:
        ax.set_ylim(-0.5,1.75)
    ax.set_xlabel(band_x)
    ax.set_ylabel(color)

    return np.array(redshift_range),np.array(number)

def color_sloan(color, mags):
    if color=='sloan_r-sloan_z':
        slope_r_m_i=-0.0192138872893
        slope_r_m_z=(1.584 * slope_r_m_i)
        slope_fit=[slope_r_m_z, 0]
        i_band0=-20.
    elif color=='sloan_g-sloan_i':
        slope_v_m_i=-0.029
        slope_g_m_i=(1.481 * slope_v_m_i)
        slope_fit=[slope_g_m_i, 0]
        i_band0=-20.
    elif color=='sloan_r-sloan_i':
        slope_rc_m_ic=-0.0192138872893
        slope_r_m_i=(1.007 * slope_rc_m_ic)
        slope_fit=[slope_r_m_i, 0]
        i_band0=-20.5
        color_err=0.18
    elif color=='sloan_g-sloan_r':
        slope_v_m_r=-0.0133824600874
        slope_g_m_r=(1.646 * slope_v_m_r)
        slope_fit=[slope_g_m_r, 0]
        i_band0=-20.5
        color_err=0.15

    band_1, band_2 = color.split("-")
    band_1_idx=filters.index(band_1)
    band_2_idx=filters.index(band_2)
    rs_models=dict()
    rs_models[color]=dict()
    for z, m in zip(zs,mags):
        #mag_1=m[band_1_idx]
        mag_2=m[band_2_idx]
        mag_1=blue_model(color,mags,z,mag_2)
        this_model=RSModel(z, mag_1, mag_2, slope_fit)
        rs_models[color][this_model.z]=this_model

    return slope_fit,i_band0,color_err,rs_models,band_1,band_2

# adding the slope for different color set that we are interested in (01_rsz_test,fit_gr_ri01.ipyn)
def blue_model(color,mags,redshift,red_mag):
    #g-r
    if color=='sloan_g-sloan_r':
        blue_mag=(0.787302458781+2.9352*redshift)+red_mag
    elif color=='sloan_r-sloan_i':
        if redshift <= 0.36:
            blue_mag=(0.348871987852+0.75340856*redshift)+red_mag
        else:
            blue_mag=(-0.210727367027+2.2836974*redshift)+red_mag
    else:
        print 'This color has not been implemented.'
    return blue_mag

def histogram_plot(xranf,numberf,df,ax=None,line=False,cbar=False):
    l2=6
    ax.set_xlim(0,0.8)
    ic2,ic3=0,0

    numbers=numberf[:6]
    numbers2=numberf[l2:]

    ax.bar(xranf[:6],numbers,width=0.05,color='red',alpha=0.5,align='center')
    ax.bar(xranf[l2:],numbers2,width=0.05,alpha=0.5,align='center')
    if cbar:
        cbar=fig.colorbar(s_m, ax=ax)
        cbar.set_label("redshift")

    if line:
        if dff_sdss.loc[ind].redshift!=-1:
            ax.axvline(dff_sdss.redshift[ind],ls='--',color='#66cc00',lw=2.,label='qso z=%.2f'%dff_sdss.redshift[ind])
        ax.axvline(xranf[:6][ic2],ls='--',color='black',lw=2.,label='red_seq g-r z=%.2f'%xranf[:6][ic2])
        ax.axvline(xranf[l2:][ic3],ls='--',color='purple',lw=2.,label='red_seq r-i z=%.2f'%xranf[l2:][ic3])

        ax.legend(loc='best',frameon=False)

    sigma,sigma2,sigma3=0.,0.,0.
    if line:
        return np.array([xranf[:6][ic2],sigma2,xranf[l2:][ic3],sigma3,dff_sdss.redshift[ind],sigma])
    else:
        return np.array([xranf[:6][ic2],sigma2,xranf[l2:][ic3],sigma3])

def save_rgb_image_extra(field, f026):
    cmd = "ds9 -zscale -crosshair %f %f wcs fk5 -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom out -saveimage final/img%s_2.eps -exit" % \
        (f026.RA0.values[0], f026.DEC0.values[0], field, field, field, field)
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))

    cmd = "ds9 -rgb -red final/coadd_c%s_i.fits -green final/coadd_c%s_r.fits -blue final/coadd_c%s_g.fits -zoom out -saveimage final/img%s_3.eps -exit" % \
        (field, field, field, field)
    print cmd
    sub = subprocess.check_call(shlex.split(cmd))
    print 'finished saving final/img%s.eps' % field

def find_offset(fname):
    with open(fname) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    band=[x.split(' ')[0][-1] for x in content[5:-1]]
    corr=[float(x.split(' ')[1]) for x in content[5:-1]]
    ecorr=[float(x.split(' ')[3]) for x in content[5:-1]]
    return zip(band,corr,ecorr), corr

def find_num(fname):
    with open(fname) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    num_2mass=content[0].split(' ')[3]
    num_star=content[3].split(' ')[1]
    chisq=content[2].split(' ')[1]
    return num_2mass,num_star,chisq

##--------
if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    filters=['sloan_r','sloan_i','sloan_z','sloan_g']
    zfs = np.arange(1.0, 6.001, 0.05)
    zf = 3.0  #formation redshift
    spacing=0.01  #spacing of redshift for resolution (0.01 is high_res, 0.05 low_res)
    zs = np.arange(0.05, 2.500001, spacing)

    new_model = ezgal.model("pisco_pipeline/pisco_exp_chab_evolved.model")
    new_model.set_normalization(filter='ks', mag=10.9, apparent=True, vega=True,z=0.023) ##normalize to Coma
    new_mags = new_model.get_apparent_mags(zf, filters=filters, zs=zs, ab=True)

    df_color=dict()
    df_color['sloan_g']='MAG_g'
    df_color['sloan_r']='MAG_r'
    df_color['sloan_i']='MAG_i'
    df_color['sloan_z']='MAG_z'

    df_colorerr=dict()
    df_colorerr['sloan_g']='MAGERR_g'
    df_colorerr['sloan_r']='MAGERR_r'
    df_colorerr['sloan_i']='MAGERR_i'
    df_colorerr['sloan_z']='MAGERR_z'

    zss=zs[0:80:5]
    norm = matplotlib.colors.Normalize(vmin=np.min(zss),vmax=np.max(zss))
    c_m = matplotlib.cm.RdYlBu
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    # Pipeline to run PISCO reduction data
    #dir = str(sys.argv[1])
    field = str(sys.argv[1])
    slrdir = 'slr_output'

#    field = 'Field054'
    df_all = pd.read_csv("/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/all_objs_list_new.csv")
    f026 = df_all[df_all["name"]==field]
    redshift=f026.redshift.values[0]
    priority=f026.priority.values[0]

    seeing=Table.read('/Users/taweewat/Dropbox/Documents/MIT/Observation/2017_1/PISCO_Jan17_seeing.csv')
    see=seeing[seeing['Field']==int(field[-3:])]['Seeing'][0]

    offset=find_offset('slr_output/star_%s.fits.offsets.list' % field)
    num_2mass,num_star,chisq=find_num('../pisco_code/slr_output/star_%s.fits.offsets.list' % field)

    #save_rgb_image_extra(field, f026)

    df = pd.read_csv(os.path.join(slrdir,'ntotal_%s.csv' % field),index_col=0)
    c5 = SkyCoord(ra=df['XWIN_WORLD'].values*u.degree, dec=df['YWIN_WORLD'].values*u.degree)
    c0 = SkyCoord(ra=f026.RA0*u.degree, dec=f026.DEC0*u.degree)
    sep = c5.separation(c0)
    cut=df[(sep.arcmin<ex.rad_A(redshift,dist=1.5)) & (df["CLASS_STAR"]<0.75)] #CLASS_STAR < 0.75
    #ncut=df[(sep.arcmin>2.5) & (df["CLASS_STAR"]<0.8)]
    print see
    print offset[1]
    fig,ax=plt.subplots(1,4,figsize=(20,5));
    fig.suptitle(field+', Redshift='+str(redshift)+', Priority='+priority+', Seeing='+str(see)+', Offset(r,i,g,z)='+str(offset[1])+', #2mass='+str(num_2mass)+', #stars='+str(num_star)+', chisq='+str(chisq))

    make_images(field,ax[0])
    xran,numbers_gr=red_seq_color_plot('sloan_g-sloan_r',cut,new_mags,ax[1])
    xran2,numbers_ri=red_seq_color_plot('sloan_r-sloan_i',cut,new_mags,ax[2])
    total_sigma=histogram_plot(np.append(xran,xran2),np.append(numbers_gr,numbers_ri),cut,ax[3])
    ax[3].axvline(redshift, color='green')


    fig.tight_layout()

    fig.savefig('plots/plot_%s.png' % (field), dpi=200)
