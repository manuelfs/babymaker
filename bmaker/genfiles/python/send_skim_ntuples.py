#! /usr/bin/env python

from __future__ import print_function

import argparse
import glob
import os
import numpy
import itertools
import subprocess
import re

import utilities

def getSkimName(cut):
    cut = cut.replace(">=","GE")
    cut = cut.replace("<=","SE")
    cut = cut.replace("&","_");
    cut = cut.replace(">","G")
    cut = cut.replace("<","S")
    cut = cut.replace("=","");
    cut = cut.replace("(","")
    cut = cut.replace(")","")
    cut = cut.replace("+","");
    cut = cut.replace("[","")
    cut = cut.replace("]","")
    cut = cut.replace("|","_");
    cut = cut.replace("$","")
    cut = cut.replace(",","_")
    cut = cut.replace("!","NOT");
    cut = cut.replace(" ","")
    cut = cut.replace("@","")
    return cut

def splitJobs(files, num_jobs):
    return [ a.tolist() for a in numpy.array_split(numpy.array(files), num_jobs) if len(a.tolist()) > 0 ]

def sendSkimJob(in_files, out_files, cut, overwrite, cache, exe_name):
    python_dir = utilities.fullPath(os.path.dirname(__file__))
    run_dir = os.path.join(os.path.dirname(out_files[0]), "run")
    utilities.ensureDir(run_dir)
    run_file = os.path.join(run_dir, exe_name)

    with open(run_file, "w") as f:
        f.write('#! /usr/bin/env python\n')
        f.write('import sys\n')
        f.write('sys.path.append("'+python_dir+'")\n')
        f.write('import subprocess\n')
        f.write('import cache\n')
        for in_file, out_file in itertools.izip(in_files, out_files):
            if os.path.exists(out_file) and not overwrite:
                continue
            if cache:
                f.write('cache.cacheRun(["'+out_file+'","'+in_file+'"],["'
                        +os.path.join(python_dir,'skim_ntuple.py')
                        +'","'+cut+'","'+out_file+'","'+in_file
                        +'"],False,10000000000,0.5,False)\n')
            else:
                f.write('subprocess.call(["'+os.path.join(python_dir,'skim_ntuple.py')
                        +'","'+cut+'","'+out_file+'","'+in_file+'"])\n')
    os.chmod(run_file, 0755)

    subprocess.call(["JobSubmit.csh","run/wrapper.sh",run_file])

def sendSkims(in_dir, num_jobs, cut, out_parent, file_tag, overwrite, cache):
    in_dir = utilities.fullPath(in_dir)
    skim_name = getSkimName(cut)

    if out_parent == None:
        dir_pat = re.compile("(.*?/cms[0-9]+/cms[0-9]+r0/babymaker/babies/[0-9]{4}_[0-9]{2}_[0-9]{2}/.*?)/")
        match = dir_pat.search(in_dir+"/")
        out_parent = match.group(0)

    out_dir = os.path.join(out_parent,"skim_"+skim_name)
        
    in_files = [ f for f in glob.glob(utilities.fullPath(os.path.join(in_dir, "*"+file_tag+"*.root"))) ]
    out_files = [ f.replace(in_dir, out_dir).replace(".root","_"+skim_name+".root") for f in in_files ]

    in_files = splitJobs(in_files, num_jobs)
    out_files = splitJobs(out_files, num_jobs)

    total_jobs = 0
    for ijob in xrange(len(in_files)):
        total_jobs += 1
        sendSkimJob(in_files[ijob], out_files[ijob], cut, overwrite, cache,
                    skim_name+"_"+file_tag+"_"+str(ijob)+"_of_"+str(num_jobs)+".py")

    print("Submitted "+str(total_jobs)+" jobs.")
    print("Output sent to {}".format(out_dir))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submits jobs to skim non-SMS ntuples.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_dir", help="Directory from which to read pre-skim ntuples. E.g. /net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/data/unskimmed/alldata")
    parser.add_argument("cut", help="Skim cut to apply.")
    parser.add_argument("out_dir", default=None, nargs="?",
                        help="Parent directory in which to place skim_XYZ directory. If omitted, attempts to use the YYYY_MM_DD/data_or_mc directory corresponding to the input directory.")
    parser.add_argument("num_jobs", type=int, nargs="?", default=100,
                        help="Number of jobs over which to divide skimming.")
    parser.add_argument("file_tag", metavar="file_tag", default="", nargs="?",
                        help="Only skim files matching %(metavar)s. Matches all files if blank.")
    parser.add_argument("-o","--overwrite", action="store_true",
                        help="Remake skimmed output file even if it already exists.")
    parser.add_argument("--cache", action="store_true",
                        help="Enable use of file caching system")
    args = parser.parse_args()

    sendSkims(args.in_dir, args.num_jobs, args.cut, args.out_dir, args.file_tag, args.overwrite, args.cache)
