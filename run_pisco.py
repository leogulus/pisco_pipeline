import os
import subprocess
import shlex
import sys

def list_file_name(dir,name,end=0):
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
    names=[]
    for file in os.listdir(dir):
        if file.startswith(name):
            if end==0:
                names.append(os.path.join(dir, file))
            else:
                if file.endswith(end):
                    names.append(os.path.join(dir, file))
    if len(names)==0:
        print 'Cannot find the files'
    return names

if __name__ == "__main__":
    all_fields=list(set([i.split('_')[0].split('/')[-1] for i in list_file_name('data','Field')]))
    for field in all_fields[0:2]:
        cmd="python pisco_pipeline/pisco_combine.py data/ %s" % field
        print cmd
        sub=subprocess.check_call(shlex.split(cmd))
