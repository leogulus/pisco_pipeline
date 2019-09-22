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
            'ut170624/', 'ut171208/', 'ut171209/', 'ut171212/','ut190412/','ut190413/']
    myReg = re.compile(r'(%s_A).*' % field)
    for di in dirs:
        diri = home + di
        for text in os.listdir(diri):
            if myReg.search(text) != None:
                filename = myReg.search(text).group()
                allfilename = diri
    return allfilename

if __name__ == "__main__":
    dirs = ['ut170103/', 'ut170104/', 'ut170619/', 'ut170621/',
            'ut170624/', 'ut171208/', 'ut171209/', 'ut171212/']

    names=[]
    home='/Users/taweewat/Documents/pisco_code/'
    names=[]
    myReg=re.compile(r'(CHIPS\d{4}[+-]\d{4})|(Field\d{3})')

    dirs=['ut190413/']
    for di in dirs:
        dir=home+di
        for text in os.listdir(dir):
            if myReg.search(text) != None:
                names.append(myReg.search(text).group())
    all_fields=list(set(names))
    all_fields_cut = all_fields[:]
    # all_fields_cut = ['CHIPS1422-2728']

    for i, field in enumerate(all_fields_cut[:]): # all_fields[:22]+all_fields[23:]):
        print str(i) + '/' + str(len(all_fields_cut)), field

        file_dir = find_fits_dir(field)[-9:]
        cmd = "python pisco_pipeline/pisco_combine_2019.py %s %s 'twilight'" % (file_dir, field)    
        print cmd
        # sub = subprocess.check_call(shlex.split(cmd))

        for band in ['i','g','r','z']:
            cmd = "python pisco_pipeline/pisco_psf_extract.py %s %s" % (field,band)
            print cmd
            sub = subprocess.check_call(shlex.split(cmd))