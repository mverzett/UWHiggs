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

if len(sys.argv) < 3 or '-h' in sys.argv or '--help' in sys.argv:
    print 'Usage ./slimNTuples other_branch_to_keep.list /hdfs/dir/with/files/to/be/skimmed'
    sys.exit(1)

branches_to_keep = sys.argv[-2]
hdfs_path        = sys.argv[-1]
if not os.path.isfile(branches_to_keep):
    print "Error! %s: no such file" % branches_to_keep
    sys.exit(1)

if not os.path.isdir(hdfs_path):
    print "Error! %s: no such directory" % hdfs_path
    sys.exit(1)

print os.getcwd()
print "looking in the code to find used branches..."

#look for row.*
proc = subprocess.Popen("grep 'row\.' %s/src/UWHiggs/?h/*.py" % (os.environ['CMSSW_BASE']), shell=True, stdout=subprocess.PIPE)
stdout, stderr = proc.communicate()
exitc  = proc.wait()
regex = re.compile(r'row\.\w+')
usedBranches = list( set([i.split('.')[-1] for i in regex.findall(stdout)]))
#print usedBranches
#look for more complex statements (getVar and getattr)
proc = subprocess.Popen("grep getattr %s/src/UWHiggs/?h/*.py | grep row" % (os.environ['CMSSW_BASE']), shell=True, stdout=subprocess.PIPE)
stdout, stderr = proc.communicate()
exitc  = proc.wait()
regex = re.compile(r'getattr\( *row *,[^\'\"]*[\'\"](?P<branch>(?:\w|\%)+)')
otherObjs = [m.group('branch') for m in regex.finditer(stdout) if m]
regex = re.compile(r'getVar\( *\w+ *,[^\'\"]*[\'\"](?P<branch>(?:\w|\%)+)')
otherObjs.extend([m.group('branch') for m in regex.finditer(stdout) if m])
otherObjs = [ i.replace('%s','*') if '%s' in i else '*'+i for i in otherObjs ]

#free some memory
del stdout

#get rid of doubles
otherObjs = list( set(otherObjs) )
#check that our starred branches don't match any non starred, in case of match discard the non starred
usedBranches = filter(lambda x: not any( ( fnmatch(x,y) for y in otherObjs) ), usedBranches) #use generator to be faster
usedBranches.extend(otherObjs)
usedBranches.extend([line.strip() for line in open(branches_to_keep, 'r')])

JOBID   = filter(lambda x: x is not '', hdfs_path.split('/'))[-1]
SAMPLES = [i.split('/')[-1] for i in glob(hdfs_path+'/*')]
##SAMPLES = ['VH_120_HWW' , 'VH_130_HWW', 'VH_H2Tau_M-120', 'VH_H2Tau_M-130']
##print SAMPLES
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

timestamp = str(int(time.mktime(time.gmtime())))
#listsdir  = '/scratch/taroni/%s' % timestamp
listsdir = '/afs/hep.wisc.edu/cms/%s/%s' % (os.environ['USER'],timestamp)
make_dir(listsdir)

#define matching function, use generator to be faster (should stop at the first match)
match  = lambda x:  any( ( fnmatch(x,y) for y in usedBranches) )
with open(listsdir+'/trees_location.list','w') as f:
    f.write('\n'.join(forest))

tot_branches  = 0
kept_branches = 0
for tree_name in forest:
    if 'Ntuple' in tree_name:
        tree = tfile.Get(tree_name)
        #Get All the branches
        matching_branches = [branch.GetName() for branch in tree.GetListOfBranches()]
        tot_branches     += len(matching_branches)
        matching_branches = filter(match, matching_branches)
        kept_branches    += len(matching_branches)
        with open(listsdir+'/%s.list' % tree_name.replace('/','_'),'w') as f:
            f.write('\n'.join(matching_branches))
        

compression_ratio = float(kept_branches) / float(tot_branches)

##newhdfs_path = hdfs_path
##print newhdfs_path.replace("mverzett","taroni")

#os.path.getsize(
print 'Compression ratio: %s' % compression_ratio
print "Beginning multi threaded copy/complession/merge..."
approx_chunk_size = 70*10**6 #~50MB
procs   = {}
#threads = os.environ['megaworkers'] if 'megaworkers' in os.environ else 2
#with open("copyMergeSlim_"+timestamp+".status",'w') as status_file:
#    logdir = os.getcwd()+'/'+timestamp
#    make_dir(logdir)
#    make_dir('/'.join(['','scratch',os.environ['USER'],'data',JOBID]) )
files_dict   = dict( [(sample, glob('/'.join([hdfs_path,sample,'','*.root'])))  for sample in SAMPLES] )
#    tot_numfiles = reduce(lambda x, y: x+y, map(len,files_dict.itervalues()) )
#    pbar  = ProgressBar(
#        widgets = [
#            FormatLabel(
#                'Copied %(value)i/' + str(tot_numfiles) + ' files. '),
#                ETA(),
#                Bar('>')],
#        maxval = tot_numfiles ).start()
#    files_already_merged = 0
run = 'source /afs/hep.wisc.edu/cms/taroni/newHiggs/src/FinalStateAnalysis/environment.sh\n'
for sample in SAMPLES:
    print 'Merging %s...' % sample
    sample_dir = '/'.join(['','scratch',os.environ['USER'],'data',JOBID,sample])
    make_dir( sample_dir )
    chunk = []
    isize  = 0
    for fname in files_dict[sample]:
        isize = isize + os.path.getsize(fname)
        
    ifile = len(files_dict[sample])/(isize*compression_ratio/(approx_chunk_size))
        
    print 'number of file for farmoutAnalysisJob '+ str(int(ifile)) +''
#    job = """slimAndMergeNtuple """+str(listsdir)+""" test"""+str(sample)+""".root """+str('/'.join([hdfs_path,sample,'']))+"""*.root
#    """
#    job_file= open("slimNtuples"+str(sample)+"","w")
#    job_file.write(job)
#    job_file.close()
##devo fare un file che contenga i file di input per tutti i sample?

    output_dir = 'srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/'+'/'.join([os.environ['USER'],JOBID+'_light',sample])    
    submit_dir = '/'.join(['','scratch',os.environ['USER'],JOBID,sample])
##    run+="""mkdir -p """+submit_dir+"""/dags
##farmoutAnalysisJobs  --infer-cmssw-path \"--submit-dir="""+submit_dir+"""/submit\" \
##\"--output-dag-file="""+submit_dir+"""/dags/dag\" \
##\"--output-dir="""+output_dir+"""\" \
##--input-files-per-job="""+str(int(ifile))+""" --shared-fs \"--input-dir="""+'/'.join([hdfs_path,sample,''])+"""\" --fwklite  \
##"""+'-'.join([JOBID,sample])+""" run_slim_and_merge.sh """+str(listsdir)+""" '$outputFileName' '$inputFileNames' \n """
    input_path=hdfs_path[5:]
    run+="""mkdir -p """+submit_dir+"""/dags
farmoutAnalysisJobs  --infer-cmssw-path \"--submit-dir="""+submit_dir+"""/submit\" \
\"--output-dag-file="""+submit_dir+"""/dags/dag\" \
\"--output-dir="""+output_dir+"""\" \
--input-files-per-job="""+str(int(ifile))+""" --shared-fs \"--input-dir=root://cmsxrootd.hep.wisc.edu/"""+'/'.join([input_path,sample,''])+"""\" --fwklite  \
"""+'-'.join([JOBID,sample])+""" run_slim_and_merge.sh """+str(listsdir)+""" '$outputFileName' '$inputFileNames' \n """
####        farmoutAnalysisJobs --infer-cmssw-path \"--submit-dir="""+str('/'.join([listsdir,sample]))+"""submit\" \
##        \"--output-dag-file="""+str('/'.join([listsdir,sample]))+"""dags/dag\" \
##        \"--output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN="""+str('/'.join([hdfs_path,sample]))+"""\" \
##        --input-files-per-job="""+str(ifile)+""" --shared-fs \"--input-dir="""+str('/'.join([hdfs_path,sample]))+"""\" \
##        slimNTuples makeQuad=1 makeTNP=1 makeH2Tau=0 makeTrilepton=1 make4L=1 \
##        rerunFSA=1 noPhotons=1 'inputFiles=$inputFileNames' 'outputFile=$outputFileName'"""

run_file = open("slimAndMergeNtuples"+JOBID+".run","w")
run_file.write(run)
run_file.close()
os.popen("chmod a+x slimAndMergeNtuples"+JOBID+".run")

sh ="""#! /bin/bash
source $CMSSW_BASE/src/FinalStateAnalysis/environment.sh
echo $@
slimAndMergeNtuple $@  """
sh_file= open("run_slim_and_merge.sh","w")
sh_file.write(sh)
sh_file.close()
os.popen("chmod a+x run_slim_and_merge.sh")
            
                
