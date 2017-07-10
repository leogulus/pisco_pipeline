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
dir='ut170619/'

if __name__ == "__main__":
    # all_fields = list(set([i.split('_')[0].split('/')[-1]
    #                        for i in list_file_name(dir, 'Field')]))
    # all_fields = list(set([i.split('/')[-1].split('_')[1].split('.')[0] for i in list_file_name('slr_output/', 'ntotal_Field')]))
    all_fields = ['Field234','Field237']#'Field037','Field042','Field045','Field084','Field058','Field059','Field060','Field074'] #'Field234','Field237'
    for i, field in enumerate(all_fields[:]):
        print i, field
        # cmd = "python pisco_pipeline/pisco_combine.py %s %s" % (dir,field)
        cmd = "python pisco_pipeline/pisco_combine.py %s %s 'twilight'" % (dir,field)
        print cmd
        sub = subprocess.check_call(shlex.split(cmd))

        # if not os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        # if os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        #     cmd = "python pisco_pipeline/pisco_photometry.py %s" % field
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))

        # if os.path.exists(os.path.join('slr_output','ntotal_'+field+'.csv')):
        #     cmd = "python pisco_pipeline/pisco_redsequence.py %s" % field
        #     print cmd
        #     sub = subprocess.check_call(shlex.split(cmd))
