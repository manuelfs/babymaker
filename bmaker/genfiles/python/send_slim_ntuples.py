#! /usr/bin/env python

from __future__ import print_function

import argparse
import glob
import os
import subprocess 
import utilities
import shutil

def sendSlimJob(skim, slim, overwrite):
    mc_dir = os.path.dirname(skim)
    skim_name = os.path.basename(skim)
    slim_name = os.path.splitext(os.path.basename(slim))[0]
    out_dir = os.path.join(mc_dir, "merged_"+slim_name+"_"+skim_name).replace("_skim_","_")
    run_dir = os.path.join(out_dir, "run")
    python_dir = utilities.fullPath(os.path.dirname(__file__))

    utilities.ensureDir(run_dir)

    tags = utilities.findBaseSampleNames(skim)
    total_jobs = 0
    for tag in tags:
        in_files = os.path.join(skim,"*"+tag+"*.root")
        out_file = os.path.join(out_dir,"mergedbaby_"+tag+"_"+skim_name+"_"+slim_name+"_nfiles_"+str(len(glob.glob(in_files)))+".root")
        if os.path.exists(out_file) and not overwrite:
            print("Keeping pre-existing "+out_file)
            continue
        subprocess.call(['JobSubmit.csh','python/cache.py','-c',slim,out_file,'-e','python/slim_ntuple.py',slim,out_file,in_files])
        total_jobs += 1

    return total_jobs

def sendSlimJobs(input_dir, skims, slims, overwrite):
    input_dir = utilities.fullPath(input_dir)
    if skims == []: skims = ["*"]
    skims = [ skim for sublist in skims for skim in glob.glob(os.path.join(input_dir, "skim_"+sublist)) ]
    slims = [ utilities.fullPath(slim) for sublist in slims for slim in glob.glob(sublist)]

    total_jobs = 0
    for slim in slims:
        for skim in skims:
            total_jobs += sendSlimJob(skim, slim, overwrite)

    print("Submitted "+str(total_jobs)+" jobs.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submits jobs to slim and merge ntuples.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--input_dir", default="/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc",
                        help="Directory containing skim directories with ntuples.")
    parser.add_argument("-k","--skims", default=["*"], nargs="+",
                        help="List of skimmed subdirectories to slim. Prefix \"skim_\" is automatically added.")
    parser.add_argument("-l","--slims", default=["txt/slim_rules/*.txt"], nargs="+",
                        help="List of slims to generate.")
    parser.add_argument("-o","--overwrite", action="store_true",
                        help="Remake slimmed output file even if it already exists")
    args = parser.parse_args()

    sendSlimJobs(args.input_dir, args.skims, args.slims, args.overwrite)
