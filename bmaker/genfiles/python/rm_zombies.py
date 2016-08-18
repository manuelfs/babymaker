#! /usr/bin/env python

import argparse
import os
import glob

import ROOT

import utilities

def killZombies(in_dirs):
    in_dirs = [ utilities.fullPath(d) for sublist in in_dirs for d in glob.glob(sublist) ]
    ROOT.gErrorIgnoreLevel = 6000
    for d in in_dirs:
        for root, dirs, files in os.walk(d):
            print "In "+root
            for f in files:
                path = os.path.join(root, f)
                if os.path.splitext(f)[1] != ".root":
                    continue
                tfile = ROOT.TFile(path, "read")
                kill = tfile.IsZombie() or not tfile.IsOpen()
                tfile.Close()
                if kill:
                    print "Removing "+path
                    os.remove(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recursively removes zombie ROOT files from given directories.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_dirs", nargs="*", default=["."],
                        help="List of directories from which to purge zombie files.")
    args = parser.parse_args()

    killZombies(args.in_dirs)
