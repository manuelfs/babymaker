#!/usr/bin/env python

###### Script to count the number of .root files in all subfolders
import os, sys, subprocess
import glob
import string
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-f", "--folder", help="Folder fo find root files in", default="./")
args = parser.parse_args()

def du(path):
    """disk usage in human readable format (e.g. '2,1GB')"""
    return subprocess.check_output(['du','-sh', path]).split()[0].decode('utf-8')
class bcolors:
    BOLD = '\033[1m'
    ENDC = '\033[0m'

print "\n==== Counting number of root files in subfolders of "+bcolors.BOLD+args.folder+ bcolors.ENDC+"\n"

subfolders = sorted([x[0] for x in os.walk(args.folder)])
for subfolder in subfolders:
    files = glob.glob(subfolder+'/*.root')
    if(len(files)>0): 
        sf_name = subfolder.split(args.folder)[1]
        print '{:>5}'.format(str(len(files)))+" .root files, size "+'{:>4}'.format(du(subfolder))+" in "+bcolors.BOLD+sf_name+ bcolors.ENDC

print 

sys.exit(0)
