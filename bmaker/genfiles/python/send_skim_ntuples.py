#! /usr/bin/env python

from __future__ import print_function

import argparse
import glob
import os
import numpy
import itertools
import subprocess

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

def sendSkimJob(in_files, out_files, cut, overwrite, exe_name):
    python_dir = utilities.fullPath(os.path.dirname(__file__))
    run_dir = os.path.join(os.path.dirname(out_files[0]), "run")
    utilities.ensureDir(run_dir)
    run_file = os.path.join(run_dir, exe_name)

    with open(run_file, "w") as f:
        f.write("#! /usr/bin/env python\n\n")

        f.write("import tempfile\n")
        f.write("import getpass\n")
        f.write("import os\n")
        f.write("import shutil\n")
        f.write("import sys\n\n")

        f.write("sys.path.append(\""+python_dir+"\")\n")
        f.write("import utilities\n")
        f.write("import skim_ntuple\n\n")

        f.write("scratch = utilities.fullPath(os.path.join(\"/scratch\",getpass.getuser()))\n")
        for in_file, out_file in itertools.izip(in_files, out_files):
            if os.path.exists(out_file) and not overwrite:
                f.write("print \"Keeping pre-existing "+out_file+"\"\n")
            else:
                f.write("with tempfile.NamedTemporaryFile(dir=scratch,prefix=\"tmp_babymaker_skim_\",suffix=\".root\") as f:\n")
                f.write("    print \"Copying "+in_file+" to \"+str(f.name)\n")
                f.write("    shutil.copyfile(\""+in_file+"\",f.name)\n")
                f.write("    skim_ntuple.skimFiles([f.name],\""+out_file+"\",\""+cut+"\","+str(not overwrite)+")\n\n")
        os.fchmod(f.fileno(), 0755)
    subprocess.call(["JobSubmit.csh","run/wrapper.sh",run_file])

def sendSkims(in_dir, num_jobs, cut, file_tag, overwrite):
    in_dir = utilities.fullPath(in_dir)
    skim_name = getSkimName(cut)
    out_dir = os.path.join(os.path.dirname(in_dir),"skim_"+skim_name)

    in_files = [ f for f in glob.glob(utilities.fullPath(os.path.join(in_dir, "*"+file_tag+"*.root"))) ]
    out_files = [ f.replace(in_dir, out_dir).replace(".root","_"+skim_name+".root") for f in in_files ]

    in_files = splitJobs(in_files, num_jobs)
    out_files = splitJobs(out_files, num_jobs)

    total_jobs = 0
    for ijob in xrange(len(in_files)):
        total_jobs += 1
        sendSkimJob(in_files[ijob], out_files[ijob], cut, overwrite,
                    skim_name+"_"+file_tag+"_"+str(ijob)+"_of_"+str(num_jobs)+".py")

    print("Submitted "+str(total_jobs)+" jobs.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submits jobs to skim non-SMS ntuples.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_dir", help="Directory from which to read pre-skim ntuples.")
    parser.add_argument("num_jobs", type=int, help="Number of jobs over which to divide skimming.")
    parser.add_argument("cut", default="baseline", nargs="?", help="Skim cut to apply.")
    parser.add_argument("file_tag", metavar="file_tag", default="", nargs="?",
                        help="Only skim files matching %(metavar)s. Matches all files if blank.")
    parser.add_argument("-o","--overwrite", action="store_true",
                        help="Remake skimmed output file even if it already exists.")
    args = parser.parse_args()

    sendSkims(args.in_dir, args.num_jobs, args.cut, args.file_tag, args.overwrite)
