import os, re
import subprocess
import shlex
import sys

from astropy.io import fits


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


def find_fits_dir(field):
    home = '/Users/taweewat/Documents/pisco_code/'
    dirs = ['ut170103/', 'ut170104/', 'ut170619/', 'ut170621/',
            'ut170624/', 'ut171208/', 'ut171209/', 'ut171212/']
    myReg = re.compile(r'(%s_A).*' % field)
    for di in dirs:
        diri = home + di
        for text in os.listdir(diri):
            if myReg.search(text) != None:
                filename = myReg.search(text).group()
                allfilename = diri
    return allfilename

# dir='data2/'
# dir='ut171208/'  #09, 12
# dir='ut170104/' #03
# dir='ut170624/' #21, 24

if __name__ == "__main__":
    dirs = ['ut170103/', 'ut170104/', 'ut170619/', 'ut170621/',
            'ut170624/', 'ut171208/', 'ut171209/', 'ut171212/']

    # all_fields = list(set([i.split('_')[0].split('/')[-1]
    #                        for i in list_file_name(dir, 'CHIPS')]))

    names=[]
    # myReg=re.compile(r'(CHIPS\d{4}-\d{4})|(Field\d{3})')
    # # print myReg
    # for text in os.listdir(dir):
    #     if myReg.search(text) != None:
    #         names.append(myReg.search(text).group())
    # all_fields=list(set(names))

    home='/Users/taweewat/Documents/pisco_code/' #09, 171208
    # ,'ut170104/','ut170619/','ut170621/','ut170624/','ut171208/','ut171209/','ut171212/']
    dirs = ['ut171208/']
    names=[]
    myReg=re.compile(r'(CHIPS\d{4}[+-]\d{4})|(Field\d{3})')
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                names.append(myReg.search(text).group())
    all_fields=list(set(names))

    # all_fields = list(set([i.split('_')[0].split('/')[-1] for i in list_file_name(dir, 'Field')]))

    # all_fields = list(set([i.split('/')[-1].split('_')[1].split('.')[0] for i in list_file_name('slr_output/', 'ntotal_Field')]))
    # all_fields = ['Field234','Field237']#'Field037','Field042','Field045','Field084','Field058','Field059','Field060','Field074'] #'Field234','Field237'
    # all_fields = ['SDSS501','SDSS603','SDSS123','PKS1353','Field227','Field228','Field229','Field230']
    # all_fields = ['Field227','Field228','Field229','Field230']
# 
    # Didn't work yet
    # UT171208: 'CHIPS0525-6938'
    # UT171209: ,
    # UT171208
    # Already DONE
    # ['CHIPS0004-4736','CHIPS0005-2758','CHIPS0022-5140','CHIPS0024-6820','CHIPS0025-5427','CHIPS0040-2902',\
    #             'CHIPS0049-4457','CHIPS0050-5249','CHIPS0127-4016','CHIPS0133-4358','CHIPS0146-3648','CHIPS0146-3711',\
    #             'CHIPS0150-3338','CHIPS0335-4715','CHIPS0342-3703','CHIPS0355-6645','CHIPS0409-2839','CHIPS0411-5149',\
    #             'CHIPS0423-3953','CHIPS0449-2859','CHIPS0449-3910','CHIPS0449-4350','CHIPS0449-4350','CHIPS0512-3257',\
    #             'CHIPS0514-5046','CHIPS0532-3917','CHIPS0535-6602','CHIPS0536-3401','CHIPS0824-3020',\
    #             'CHIPS0936-3342','CHIPS0957-7554','CHIPS1009-3015','CHIPS1034-2837','CHIPS1036-3513','CHIPS2223-3455',\
    #             'CHIPS2227-4333','CHIPS2228-3220','CHIPS2240-5231','CHIPS2243-3034','CHIPS2245-4931','CHIPS2251-3210',\
    # 'CHIPS2251-3827','CHIPS2254-3635','CHIPS2303-6807','CHIPS2306-3439','CHIPS2307-4236','CHIPS2311-4718',\
    # 'CHIPS2325-4800','CHIPS2333-6144','CHIPS2348-3831']
    #
    # UT1712089
    # Already DONE
    # 'CHIPS0004-2902','CHIPS0112-2919','CHIPS0115-3047','CHIPS0148-2238','CHIPS0152-5028','CHIPS0153-3143',\
    #             'CHIPS0157-1043','CHIPS0206-7148','CHIPS0209-6810','CHIPS0219-3626','CHIPS0222-4159','CHIPS0229-5232',\
                # 'CHIPS0253-5441','CHIPS0300-3413','CHIPS0302-2758','CHIPS0303-2407','CHIPS0304-3556','CHIPS0310-3723'
                # 'CHIPS0316-2247','CHIPS0325-4926','CHIPS0522-1745','CHIPS0552-2336','CHIPS0609-0247'
                # 'CHIPS0821-0815','CHIPS0827-2026','CHIPS0847-0703','CHIPS0849-1721',\
                # 'CHIPS0920-2257','CHIPS0934-1721','CHIPS1011-0505','CHIPS1102-3031','CHIPS1141-1407','CHIPS1142-1422',\
                # 'CHIPS1147-1252','CHIPS1205-2633','CHIPS2101-6132','CHIPS2127-4725','CHIPS2133-2815','CHIPS2139-3911',\
                # 'CHIPS2141-3729','CHIPS2148-2735','CHIPS2148-3715','CHIPS2210-5508','CHIPS2211-3707','CHIPS2216-2803',\
                # 'CHIPS2217-3034','CHIPS2218-2947','CHIPS2221-2804', 'CHIPS2246-2854','CHIPS2249-2808','CHIPS0724-0715','CHIPS0745-0714'
    # Not Yet
    # UT171209
    # all_fields=['CHIPS2311-4718']
    # Already DONE
    # 'CHIPS0003-2521','CHIPS0012-1628','CHIPS0015-1513','CHIPS0018-1840','CHIPS0050-1412',\
                # 'CHIPS0106-1149','CHIPS0106-2358','CHIPS0107-1310','CHIPS2333-2407','CHIPS2340-2302','CHIPS2349-2352','CHIPS2357-1125','CHIPS0118-1430','CHIPS2317-1443'
    # UT171212
    # CHIPS1011-0505: still some problems
    #Still problem with CHIPS2357-1125

# 'Field018','Field020','Field021','Field022',
    # all_fields=['Field102']#,'Field025','Field026','Field027','Field028','Field029','Field030','Field033','Field034','Field036','Field037','Field038','Field039','Field040','Field042','Field044','Field045'] #
    # 32 Field087
    # Not done
    # 40 Field103

    # UT171208 13 [CHIPS0025-5427] (px=0.12), 22 [CHIPS0525-6938] (already: dense field), 49 [CHIPS0536-3401] (px=0.12)
    # UT171209 10 [CHIPS0302-2758] (px=0.12), 17 [CHIPS1011-0505] (already,x), 20 [CHIPS0153-3143] (px=0.12)
    # UT171212 0 [CHIPS0012-1628] (already), 4 [CHIPS0018-1840](px=0.12, bad corner), 7 [CHIPS0132-1608] (already: bright stars)
    # UT170104 40 Field103,

    ## PHOTOMETRY
    # no WCS solution: UT171208 22 CHIPS0525-6938
    # Too Bright: UT171209 17 CHIPS1011-0505, UT171212 7 CHIPS0132-1608
    # not enough star in 2MASS: UT171212 0 CHIPS0012-1628

    # Problem with PSF Photometry
    # UT171208 2 CHIPS2303-6807, 6 CHIPS0449-2859, 7 CHIPS0150-3338, 8 CHIPS0342-3703, 13 CHIPS0025-5427, 14 CHIPS2243-3034, 17 CHIPS0409-2839, 21 CHIPS1009-3015, 33 CHIPS0449-3910, 48 CHIPS0127-4016, 49 CHIPS0536-3401

    # all_fields=['Field026'] #['Field179']
    # all_fields=['PKS1353'
    # UT170103 = array(['Field143', 'Field091', 'Field135',
    #                     'Field136', 'Field132', 'Field123', 'Field222',
    #                     'Field206', 'Field202',
    #                     'Field068', 'Field066', 'Field063', 'Field190',
    #                     'Field197', 'Field194', 'Field211', 'Field213',
    #                     'Field217', 'Field216', 'Field117', 'Field173',
    #                     'Field174', 'Field187', 'Field109', 'Field154',
    #                     'Field089', 'Field082', 'Field083']
    
    # UT170103
    # not finished: 'Field143', 
    # not enough star: 'Field089'
    # bad calibration: 'Field082'
    # finished: Field099, Field052, Field053, Field091, Field136, Field132, 'Field123', 'Field222', 'Field135'
    #                     'Field206', 'Field202' Field068, Field066, Field063, Field190, Field197, Field211,
    # Field213, 'Field217', 'Field216', Field117, Field173, Field174, Field187, 'Field109', 'Field154', 'Field083'
    # working SLR:
    # all_fields = ['Field083']
    # 'Field044'
    #UT170104
    # all_fields = ['Field092', 'Field036', 'Field056', 'Field055', 'Field204', 'Field205', 'Field201',\
    #               'Field198', 'Field219', 'Field218', 'Field071', 'Field073', 'Field075', 'Field182', 'Field151']

    # UT170619
    all_fields = ['Field039','Field225','Field226','Field022','Field020','Field021','Field025','Field048','Field124',\
    'Field121','Field045','Field046','Field047','Field058','Field233','Field038','Field059','Field292',\
    'Field018','Field210','Field212','Field115','Field072','Field074','Field077','Field076','Field269',\
    'Field266','Field279','Field274','Field088','Field084','Field085','Field087']

    # ['Field039','Field233','Field084','Field210','Field074']

    # all_fields = ['Field137']
    # ,Field048'Field039'] #['Field201'] #['Field166']
    # CHIPS0423-3953o, CHIPS2325-4800o, CHIPS2251-3827o, CHIPS2307-4236o
    # CHIPS0335-4715o, CHIPS2348-3831o, CHIPS0535-6602o
    # Day8: CHIPS0018-1840, CHIPS2357-1125o, CHIPS2317-1443o, CHIPS0015-1513o, 
    # CHIPS0106-2358o, CHIPS0118-1430o, CHIPS00501412, CHIPS0122-2646,
    # 'CHIPS0003-2521o', 'CHIPS2349-2352o', 'CHIPS2340-2302o', 'CHIPS0106-1149o', 'CHIPS2333-2407o', 'CHIPS0132-1608o'
    # Day7: 'CHIPS0148-2238', 'CHIPS1011-0505', 'CHIPS0303-2407',
    # 'CHIPS0824-3020', 'CHIPS0936-3342', 'CHIPS1009-3015', 'CHIPS1036-3513','CHIPS0005-2758'
    # Day6: all_fields=[u'CHIPS2223-3455',u'CHIPS2240-5231',u'CHIPS2245-4931',u'CHIPS2251-3827',u'CHIPS2333-6144',\
    # u'CHIPS0049-4457',u'CHIPS0127-4016',u'CHIPS0133-4358',u'CHIPS0146-3648',u'CHIPS0146-3711',u'CHIPS0423-3953',\
    # u'CHIPS0449-2859',u'CHIPS0532-3917',u'CHIPS0535-6602',u'CHIPS0824-3020',u'CHIPS0936-3342',u'CHIPS0957-7554',\
    # u'CHIPS1009-3015',u'CHIPS1036-3513']
    # all_fields = ['CHIPS0118-1430', 'CHIPS0003-2521']
    # SDSS123, SDSS501, SDSS603  # Field237, Field234, Field103
    all_fields = ['SDSS603']
    print len(all_fields)
    for i, field in enumerate(all_fields[:]):#all_fields[:22]+all_fields[23:]):
        print i, field

        # if (field=='CHIPS0012-1628')|(field=='CHIPS0018-1840')|(field=='Field103')|(field=='CHIPS0525-6938')|(field=='Field089'):
        #     continue

        file_dir = find_fits_dir(field)[-9:]
        if (file_dir == 'ut170103/') | (file_dir == 'ut170104/'):
            cmd = "python pisco_pipeline/pisco_combine.py %s %s" % (
                file_dir, field)
        else:
            cmd = "python pisco_pipeline/pisco_combine.py %s %s 'twilight'" % (
                file_dir, field)
        print cmd
        sub = subprocess.check_call(shlex.split(cmd))

        # cmd = "python pisco_pipeline/pisco_combine_dec08.py %s %s 'twilight'" % (dir,field)

        # myReg=re.compile(r'%s_A_\d{1,4}\.fits'%field)
        # for root, dirs, files in os.walk('/Users/taweewat/Documents/pisco_code/'):
        #     for file in files:
        #         if myReg.search(file) != None:
        #             seeing=float(fits.open(root+'/'+myReg.search(file).group())[0].header['FWHM1'])
        # print seeing
        # raw_input("Press Enter to continue...")


        # for band in ['i','g','r','z']:
        #     cmd = "python pisco_pipeline/pisco_psf_extract.py %s %s" % (field,band)
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))


        # cmd = "python pisco_pipeline/pisco_photometry_v3.py %s" % field
        # print cmd
        # sub = subprocess.check_call(shlex.split(cmd))

        # cmd = "python pisco_pipeline/pisco_star_galaxy_bleem.py %s" % field
        # print cmd
        # sub = subprocess.check_call(shlex.split(cmd))

        # if not os.path.isfile(os.path.join('slr_output', 'all_psf_%s.fits'%field)):
        #     cmd = "python pisco_pipeline/pisco_photometry_v4.py %s" % field  #CHIPS0005-2758
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))

        # if not os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        # if os.path.exists(os.path.join('slr_output','all_psf_'+field+'.fits')):
        # if not os.path.exists(os.path.join('slr_output','ntotal_psf_'+field+'.csv')):
        #     cmd = "python pisco_pipeline/pisco_photometry_psf.py %s" % field
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))


        # if os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        #     cmd = "python pisco_pipeline/pisco_redsequence.py %s" % field
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))
