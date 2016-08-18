#!/usr/bin/env python

from __future__ import print_function

import glob
import string
import os
import sys
import ROOT

## Finding basename for each dataset
def findBaseSampleNames(folder):
    infiles = set()
    for file in glob.glob(folder+'/*.root'):
        tag = file.split("RunII")[0]
        if ("TTJets_Tune" not in tag) and ("WJetsToLNu_Tune" not in tag) and ("DYJetsToLL_M-50_Tune" not in tag): 
            tag = tag.split("Tune")[0]
        tag = tag.split("13TeV")[0]
        tag = tag.split("-PromptReco")[0]
        tag = tag.split("_runs")[0]
        tag = tag.split("pythia")[0]
        tag = tag.split("baby_")[1]
        tag = tag.split("__")[0]
        if tag[0] != '_': tag = "_"+tag
        if tag[-1] != '_': tag = tag+"_"
        infiles.add(tag)
    sortedfiles = list()
    for file in infiles:
        sortedfiles.append(file)
    sortedfiles = sorted(sortedfiles)

    return sortedfiles

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def ePrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def flush():
    sys.stdout.flush()
    sys.stderr.flush()

class NonROOTFileError(Exception):
    def __init__(self, path):
        self.path = path
    def __str__(self):
        return self.path+" is not a ROOT file"

class ROOTOpenError(Exception):
    def __init__(self, path, mode):
        self.path = path
        self.mode = mode
    def __str__(self):
        return "Could not open "+self.path+" in "+self.mode+" mode"

class ROOTFile(object):
    def __init__(self, path, mode):
        if os.path.splitext(path)[1] != ".root":
            raise NonROOTFileError(path)
        self.path = path
        self.mode = mode
    def __enter__(self):
        self.file = ROOT.TFile(self.path, self.mode)
        if self.file.IsZombie() or not self.file.IsOpen():
            raise ROOTOpenError(self.path, self.mode)
        return self.file
    def __exit__(self, type, value, traceback):
        self.file.Close()

