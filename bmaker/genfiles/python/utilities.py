#!/usr/bin/env python

import glob
import string

## Finding basename for each dataset
def findBaseSampleNames(folder):
    infiles = set()
    for file in glob.glob(folder+'/*.root'):
        tag = file.split("RunII")[0]
        if ("TTJets_Tune" not in tag) and ("WJetsToLNu_Tune" not in tag) and ("DYJetsToLL_M-50_Tune" not in tag): 
            tag = tag.split("Tune")[0]
        tag = tag.split("13TeV")[0]
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
