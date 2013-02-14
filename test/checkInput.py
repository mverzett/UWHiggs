#! /bin/env python

import os
import sys
import subprocess
import re
from glob import glob
from fnmatch import fnmatch
from hashlib import md5
import time
from progressbar import ETA, ProgressBar, FormatLabel, Bar

def f7(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set( x for x in seq if x in seen or seen_add(x) )
    # turn the set into a list (as requested)
    return list( seen_twice )


if len(sys.argv) < 1 or '-h' in sys.argv or '--help' in sys.argv:
    print 'Usage ./checkInput.py /scratch/dir/to/file.out'
    sys.exit(1)

scratch_path        = sys.argv[-1]

if not os.path.isdir(scratch_path):
    print "Error! %s: no such directory" % scrathch_path
    sys.exit(1)


JOBID   = filter(lambda x: x is not '', scratch_path.split('/'))[-1]
SAMPLES = [i.split('/')[-1] for i in glob(scratch_path+'/*')]

for sample in SAMPLES:
    RUN = [i.split('/')[-1] for i in glob(scratch_path+sample+'/submit/*/*.out')]
    samplesfiles=[]
    for run in RUN:
        print 'reading file '+scratch_path+sample+'/submit/'+run[:-4]+'/'+run
        myfile=open(scratch_path+sample+'/submit/'+run[:-4]+'/'+run)
        for line in myfile:
            if line.find('argv read by .cc')==0:
                INPUTFILES= line[line.find('/store'):].split(',')
                print 'input file list length: ' +str(len(INPUTFILES))
                if len(f7(INPUTFILES))!=0 :  print 'file(s) ' +str(f7(INPUTFILES)) +' used more than once'
                samplesfiles.extend(INPUTFILES)
    if len(f7(samplesfiles))!=0 : print  'file(s) '+str(f7(samplesfiles))+ ' used more than once'
