import os 
import glob
import string

# List of datasets to run over, could be either MC or data
# Enter dataset name exactly as listed in CMS DAS
datasets = []
datasets.append('/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM')
datasets.append('/SingleMuon/Run2015B-PromptReco-v1/MINIAOD')

# Directory to dump all condor-related logs, schell and cmd files
subdir = "sub/"
if not (os.path.exists(os.getcwd()+'/'+subdir)):
    os.mkdir(os.getcwd()+'/'+subdir)

for ds in datasets:
    # parse the dataset name and guess the path on hadoop to create the input file list
    isdata = False # in case we want to pass this as argument to the job later
    path = ''
    tags = string.split(ds,'/')
    dsname = tags[1]
    campaign = (string.split(tags[2],'-'))[0]
    reco = tags[2][len(campaign)+1:]
    filetype = tags[3]
    if 'PromptReco' in ds:
        path = '/'.join(['/mnt/hadoop/cms/store/data',campaign,dsname,filetype,reco,'*/*/*/*/*root'])
        isdata = True
    else:
        path = '/'.join(['/mnt/hadoop/cms/store/mc',campaign,dsname,filetype,reco,'*/*root'])

    filelist = glob.glob(path)
    nfiles = len(filelist)
    if nfiles==0:
        print "\033[93m WARNING: "+ds+" not found! Skip dataset. \033[0m"
        continue
    
    fnm = '_'.join(['flist',dsname,campaign,reco+'.txt'])
    f = open(subdir+'/'+fnm,"w")
    for file in filelist:
        file = string.replace(file,'/mnt/hadoop/cms','')
        f.write(file+'\n')
    f.close()
