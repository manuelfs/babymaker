#! /usr/bin/env python

import argparse
import glob
import os
import tempfile
import shutil

import utilities

def singleFileDelete(path):
    dirname = os.path.dirname(path)
    prefix = os.path.splitext(os.path.basename(path))[0]+"_TMP_"
    with tempfile.NamedTemporaryFile(dir=dirname, prefix=prefix, suffix=".root") as tmp:
        try:
            with utilities.ROOTFile(path, "read") as orig:
                tree = orig.Get("tree")
                if not tree:
                    utilities.ePrint("Could not find tree in "+path+". Skipping.")
                    return
                with utilities.ROOTFile(tmp.name, "recreate") as outfile:
                    clone = tree.CloneTree(-1, "fast")
                    clone.Write()
            print "Deleting treeglobal from "+path
            shutil.copy(tmp.name, path)
        except (utilities.ROOTOpenError, utilities.NonROOTFileError) as r:
            utilities.ePrint("Could not open "+path+". Skipping.")
            return

def recursiveDelete(file_dir):
    try:
        singleFileDelete(file_dir)
        return
    except utilities.NonROOTFileError:
        pass
    except IOError as e:
        #Directory, not file
        if e.errno != 21:
            raise

    for root, dirs, files in os.walk(file_dir):
        for d in dirs:
            recursiveDelete(os.path.join(root, d))
        for f in files:
            try: singleFileDelete(os.path.join(root, f))
            except utilities.NonROOTFileError: pass

def deleteTreeglobal(in_files):
    in_files = [ utilities.fullPath(f) for sublist in in_files for f in glob.glob(sublist) ]
    for file_dir in in_files:
        recursiveDelete(file_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deletes treeglobal from provided list of files and directories, recursing through subdirectories.")
    parser.add_argument("files", nargs="+",
                        help="List of ROOT files from which to remove treeglobal. Directories will be recursively search for ROOT files.")
    args = parser.parse_args()

    deleteTreeglobal(args.files)
