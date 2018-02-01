import os
import subprocess
import shlex
import sys


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



dir='data2/'
dir='ut171212/'
# dir='ut170624/'
# dir='ut171209/'

if __name__ == "__main__":
    all_fields = list(set([i.split('_')[0].split('/')[-1]
                           for i in list_file_name(dir, 'Field')]))
    # all_fields = list(set([i.split('/')[-1].split('_')[1].split('.')[0] for i in list_file_name('slr_output/', 'ntotal_Field')]))
    # all_fields = ['Field234','Field237']#'Field037','Field042','Field045','Field084','Field058','Field059','Field060','Field074'] #'Field234','Field237'
    # all_fields = ['SDSS501','SDSS603','SDSS123','PKS1353','Field227','Field228','Field229','Field230']
    # all_fields = ['Field227','Field228','Field229','Field230']

    # Didn't work yet
    # UT171208: 'CHIPS0525-6938'
    # UT171209: 'CHIPS0724-0715','CHIPS0745-0714'

    # all_fields=[]

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
                # 'CHIPS2217-3034','CHIPS2218-2947','CHIPS2221-2804', 'CHIPS2246-2854','CHIPS2249-2808'

    # Not Yet
    #
    # UT171209
    # all_fields=[,\
    #             ,,\
    #             ,\
    #             ]
    # DIDN't WORK: 'CHIPS0116-1136' (cCHIPS0116-1136_A_44_i get a wrong wcs solution),
    # 'CHIPS0118-1430' [A] one get a wrong wcs solution
    # 'CHIPS2317-1443' one [A] and one [B] get a wrong wcs solution
    # Already DONE
    # 'CHIPS0003-2521','CHIPS0012-1628','CHIPS0015-1513','CHIPS0018-1840','CHIPS0050-1412',\
                # 'CHIPS0106-1149','CHIPS0106-2358','CHIPS0107-1310','CHIPS2333-2407','CHIPS2340-2302','CHIPS2349-2352','CHIPS2357-1125'
    # UT171212
    all_fields=['CHIPS0118-1430','CHIPS2317-1443']

    # 32 Field087

    # Not done
    # 40 Field103

    print len(all_fields)
    for i, field in enumerate(all_fields[:]):
        print i, field
        # cmd = "python pisco_pipeline/pisco_combine.py %s %s" % (dir,field)
        cmd = "python pisco_pipeline/pisco_combine.py %s %s 'twilight'" % (dir,field)
        # cmd = "python pisco_pipeline/pisco_combine_dec08.py %s %s 'twilight'" % (dir,field)
        print cmd
        sub = subprocess.check_call(shlex.split(cmd))

        # cmd = "rm wcs/*new"
        # print cmd
        # sub = subprocess.check_call(shlex.split(cmd))
        #
        # cmd = "rm new_fits/*new.fits"
        # print cmd
        # sub = subprocess.check_call(shlex.split(cmd))


        # if not os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        # if os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        #     cmd = "python pisco_pipeline/pisco_photometry.py %s" % field
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))
        #
        # if os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        #     cmd = "python pisco_pipeline/pisco_redsequence.py %s" % field
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))
