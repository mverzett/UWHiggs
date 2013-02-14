#! /bin/env python

import os
import sys
import subprocess
import re
from glob import glob
from fnmatch import fnmatch
from hashlib import md5
import time
import collections
from progressbar import ETA, ProgressBar, FormatLabel, Bar

if len(sys.argv) < 1 or '-h' in sys.argv or '--help' in sys.argv:
    print 'Usage ./saveGenNtuple hdfs_path'
    sys.exit(1)

#branches_to_keep = sys.argv[-2]
hdfs_path        = sys.argv[-1]
#if not os.path.isfile(branches_to_keep):
#    print "Error! %s: no such file" % branches_to_keep
#    sys.exit(1)

if not os.path.isdir(hdfs_path):
    print "Error! %s: no such directory" % hdfs_path
    sys.exit(1)

print os.getcwd()
print "looking in the files to find gen quantities..."

JOBID   = filter(lambda x: x is not '', hdfs_path.split('/'))[-1]
SAMPLES = [i.split('/')[-1] for i in glob(hdfs_path+'/*')]
sample_tfile = glob('/'.join([hdfs_path,SAMPLES[0],'','*.root']))[0]
print 'Using %s as sample file to find trees and matching branches...' % sample_tfile

import ROOT
def GetContent(dir):
    tempList = dir.GetListOfKeys()
    retList = []
    for it in range(0,tempList.GetSize()):
       retList.append(tempList.At(it).ReadObj())
    return retList

def FindTrees( directory, dirName, objectList ):
    dirContent = GetContent(directory)
    for entry in dirContent:
        if type(entry) is ROOT.TDirectory or type(entry) is ROOT.TDirectoryFile:
            subdirName = os.path.join(dirName,entry.GetName())
            FindTrees(entry, subdirName,objectList)
        elif entry.InheritsFrom('TTree'):
            pathname = os.path.join(dirName,entry.GetName())
            objectList.append(pathname)

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

#Open sample file
tfile = ROOT.TFile.Open(sample_tfile)

#find trees inside
forest = []
FindTrees( tfile, '', forest )
tot_branches  = []

for tree_name in forest:
    if 'Ntuple' in tree_name:
        tree = tfile.Get(tree_name)
        #Get All the branches
        matching_branches =[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('Gen') is not -1 ]
        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('ComesFromHiggs') is not -1 ]
        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('recoil') is not -1 ]
##        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('ToMETDPhi') is not -1 ]
        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('Phi') is not -1 and branch.GetName().find('Phi')+3 == len(branch.GetName()) ]
        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('tDecayMode') is not -1 ]
        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('eChargeIdMedium') is not -1 ]
        matching_branches +=[branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName().find('bCSVVeto') is not -1 ]
        tot_branches     += matching_branches
d = {}
for x in tot_branches:
        d[x] = 1
        kept_branches=list(d.keys()) 


#print kept_branches
with open('branchesToKeep.list', 'w') as f:
    f.write('\n'.join(kept_branches))
    
