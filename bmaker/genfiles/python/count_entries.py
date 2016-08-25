#! /usr/bin/env python

import argparse
import os
import glob
import locale

import ROOT

from utilities import *

def countFile(cut, path, counts):
    counts[path] = (0,0)
    try:
        with ROOTFile(path, "read") as f:
            counts[path] = (0,1)
            tree = f.Get("tree")
            if not tree: return
            if not cut:
                counts[path] = (tree.GetEntries(), 1)
            else:
                counts[path] = (tree.GetEntries(cut), 1)
    except (NonROOTFileError, ROOTOpenError):
        pass

def countRecursive(cut, path, counts):
    if os.path.isfile(path):
        countFile(cut, path, counts)
        return

    for root, dirs, files in os.walk(path, topdown=False):
        num_entries = 0
        num_files = 0
        for f in files:
            p = fullPath(os.path.join(root, f))
            countFile(cut, p, counts)
            num_entries += counts[p][0]
            num_files += counts[p][1]
        for d in dirs:
            p = fullPath(os.path.join(root, d))
            num_entries += counts[p][0]
            num_files += counts[p][1]
        counts[fullPath(root)] = (num_entries, num_files)

def printCounts(path, counts, verbose):
    for root, dirs, files in os.walk(path):
        p = fullPath(root)
        level = root.replace(path, "").count(os.sep)
        indent = " "*2*level
        print(("{}"+Term.BOLD+Term.BLUE+"{}"+Term.END+Term.END+" [{} entries in {} ROOT files]")
              .format(indent, os.path.basename(root), 
                      locale.format("%d", counts[p][0], grouping=True),
                      locale.format("%d", counts[p][1], grouping=True)))
        if not verbose: continue
        subindent = " "*2*(level+1)
        for f in files:
            p = fullPath(os.path.join(root, f))
            print(("{}"+Term.BOLD+Term.GREEN+"{}"+Term.END+Term.END+" [{} entries]")
                  .format(subindent, f, 
                          locale.format("%d", counts[p][0], grouping=True)))

def countEntries(cut, file_dirs, verbose):
    locale.setlocale(locale.LC_ALL, "en_US")

    file_dirs = [ fullPath(f) for sublist in file_dirs for f in glob.glob(sublist) ]
    counts = dict()

    for path in file_dirs:
        countRecursive(cut, path, counts)
        printCounts(path, counts, verbose)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Counts the number of entries in ROOT files and directories",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c","--cut", default=None,
                        help="Count only arguments passing provided cut")
    parser.add_argument("-v","--verbose", action="store_true",
                        help="Print entries for individual files.")
    parser.add_argument("files", default=["."], nargs="*",
                        help="List of files and directories in which to count root files")
    args = parser.parse_args()

    countEntries(args.cut, args.files, args.verbose)
