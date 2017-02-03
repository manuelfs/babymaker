#! /usr/bin/env python

import argparse
import glob
import errno
import os
import pwd
import grp

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def lock(file_dir, unlock_dirs):
    ruid = 0
    rgid = 0
    try: ruid = pwd.getpwnam("root").pw_uid
    except KeyError: pass
    try: rgid = grp.getgrnam("root").gr_gid
    except KeyError: pass

    try:
        with open(file_dir) as f:
            try: os.fchmod(f.fileno(), 0444)
            except OSError as e:
                if e.errno != errno.EPERM: raise
            try: os.fchown(f.fileno(), ruid, -1)
            except OSError as e:
                if e.errno != errno.EPERM: raise
            try: os.fchown(f.fileno(), -1, rgid)
            except OSError as e:
                if e.errno != errno.EPERM: raise
    except IOError as e:
        if e.errno != errno.EISDIR: raise
        else:
            try: os.chmod(file_dir, 0777 if unlock_dirs else 0555)
            except OSError as e:
                if e.errno != errno.EPERM: raise
            try: os.chown(file_dir, ruid, -1)
            except OSError as e:
                if e.errno != errno.EPERM: raise
            try: os.chown(file_dir, -1, rgid)
            except OSError as e:
                if e.errno != errno.EPERM: raise

def lockFiles(file_dirs, unlock_dirs):
    file_dirs = [ fullPath(f) for sublist in file_dirs for f in glob.glob(sublist) ]

    for file_dir in file_dirs:
        for root, dirs, files in os.walk(file_dir):
            for f in files:
                lock(os.path.join(root,f), unlock_dirs)
            lock(root, unlock_dirs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sets directories to R+W and files to R permissions. Also tries to make the owner root.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file_or_dir", default = ["."], nargs = "*",
                        help="List of files and/or directories to lock.")
    parser.add_argument("-u","--unlock_dirs", action="store_true", help="Leave directories unlocked for addition/removal of files")
    args = parser.parse_args()

    lockFiles(args.file_or_dir, args.unlock_dirs)
