#! /usr/bin/env python

from __future__ import print_function

import argparse
import glob
import os
import subprocess 
import utilities
import shutil

def sendSlimJob(skim, slim, overwrite, out_dir):
    mc_dir = os.path.dirname(skim)
    skim_name = os.path.basename(skim)
    slim_name = os.path.splitext(os.path.basename(slim))[0]
    if out_dir == None:
        out_dir = os.path.join(mc_dir, "merged_"+slim_name+"_"+skim_name).replace("_skim_","_")
    else:
        out_dir = os.path.join(out_dir, "merged_"+slim_name+"_"+skim_name).replace("_skim_","_")
    run_dir = os.path.join(out_dir, "run")
    utilities.ensureDir(run_dir)

    tags = utilities.findBaseSampleNames(skim)
    total_jobs = 0
    for tag in tags:
        #if "TTJets_SingleLeptFromT_" not in tag: continue
        in_files = os.path.join(skim,"*"+tag+"*.root")
        out_name = "mergedbaby_"+tag+"_"+skim_name+"_"+slim_name+"_nfiles_"+str(len(glob.glob(in_files)))
        out_file = os.path.join(out_dir,out_name+".root")
        run_file = os.path.join(run_dir,out_name+".sh")
        
        if os.path.exists(out_file) and not overwrite:
            print("Keeping pre-existing "+out_file)
            continue
        with open(run_file, "wb") as f:
            f.write("#! /bin/bash\n\n")
            #f.write("python/cache.py -c "+slim+" "+out_file+" -e python/slim_ntuple.py "+slim+" "+out_file+" "+in_files+"\n")
            f.write("python/slim_ntuple.py "+slim+" "+out_file+" "+in_files+"\n")
            os.fchmod(f.fileno(),0755)
        subprocess.call(["JobSubmit.csh","run/wrapper.sh",run_file])
        total_jobs += 1

    return total_jobs

def sendSlimJobs(input_dir, skims, slims, overwrite, output_dir):
    input_dir = utilities.fullPath(input_dir)
    if skims == []: skims = ["*"]
    skims = [ utilities.fullPath(skim) for sublist in skims for skim in glob.glob(os.path.join(input_dir, "skim_"+sublist)) ]
    slims = [ utilities.fullPath(slim) for sublist in slims for slim in glob.glob(sublist)]

    total_jobs = 0
    for slim in slims:
        for skim in skims:
            total_jobs += sendSlimJob(skim, slim, overwrite, output_dir)

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
    parser.add_argument("--overwrite", action="store_true",
                        help="Remake slimmed output file even if it already exists")
    parser.add_argument("--output_dir", default=None,
                        help="Directory in which to put slimmed subdirectories. Uses input directory if omitted.")
    args = parser.parse_args()

    sendSlimJobs(args.input_dir, args.skims, args.slims, args.overwrite, args.output_dir)
