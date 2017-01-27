#! /usr/bin/env python

from __future__ import print_function

import argparse
import glob
import os
import subprocess
import time

def FullPath(path):
    return os.path.abspath(os.path.expanduser(path))

def Resubmit(dir_list, interval):
    dir_list = [ FullPath(d) for sub_list in dir_list for d in glob.glob(sub_list) ]

    print("Resubmitting the following projects every {} seconds:".format(interval))
    for d in dir_list:
        print(d)
    print("")
    
    while(True):
        for d in dir_list:
            subprocess.call(["crab","resubmit","-d",d])
        time.sleep(interval)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Resubmits all failed jobs in specified directories repeatedly after specified interval.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("resub_dir", nargs="*", help="CRAB directories to resubmit")
    parser.add_argument("--interval", type=float, default=600, help="Seconds between resubmission attempts.")

    args = parser.parse_args()
    Resubmit(args.resub_dir, args.interval)
