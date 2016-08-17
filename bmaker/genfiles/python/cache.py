#! /usr/bin/env python

from __future__ import print_function

import argparse
import glob
import subprocess
import os
import tempfile
import shutil
import time
import sys

import utilities

def expand(files):
    expanded = []
    for f in files:
        globbed = glob.glob(f)
        if len(globbed) > 0:
            for g in globbed:
                expanded.append(utilities.fullPath(g))
        else:
            expanded.append(utilities.fullPath(f))
    return expanded

def isNetFile(path):
    return path.startswith("/net/")

def cachePath(path):
    cache_root = utilities.fullPath("/scratch/babymaker")
    return os.path.join(cache_root, path[5:])

def lastTime(path):
    return max(os.path.getctime(path), os.path.getmtime(path), os.path.getatime(path))

def touch(path, times):
    with open(path, "a"):
        os.utime(path, times)

def mapFiles(execute, file_map):
    #Replace executable arguments with cached equivalent

    expanded_args = []
    for arg in execute:
        globbed = glob.glob(arg)
        if len(globbed) > 0:
            #Argument represents file(s)
            for f in globbed:
                expanded_args.append(utilities.fullPath(f))
        else:
            expanded_args.append(arg)

    execute = []
    inv_file_map = dict((cached,net) for net,cached in file_map.iteritems())
    for arg in expanded_args:
        if ( arg in file_map and os.path.exists(file_map[arg])
             and (not os.path.exists(arg) or os.path.getmtime(file_map[arg])>=os.path.getmtime(arg)) ):
            #Check if generated cache for file
            execute.append(file_map[arg])
        elif isNetFile(arg):
            #Check if pre-existing cache
            cache_path = cachePath(arg)
            if ( os.path.exists(cache_path)
                 and (not os.path.exists(arg) or os.path.getmtime(cache_path)>=os.path.getmtime(arg)) ):
                execute.append(cache_path)
                inv_file_map[cache_path] = arg
            else:
                execute.append(arg)
        else:
            execute.append(arg)

    return execute, inv_file_map

def netCopy(src, dst):
    print("Copying "+src+" to "+dst+"\n")
    shutil.copy(src, dst)
    dst_time = os.path.getmtime(dst)
    src_time = os.path.getmtime(src)
    if src_time < dst_time:
        #Make sure cached copy is newer than /net copy
        touch(src, (dst_time+1., dst_time+1.))

def removeOldCache(file_map):
    #Deletes oldest cached file
    found_file = False
    oldest_mod_time = 0
    oldest_path = ""
    for root, dirs, files in os.walk(utilities.fullPath("/scratch/babymaker")):
        for f in files:
            path = os.path.join(root, f)
            if path in file_map.values(): continue
            mod_time = lastTime(path)
            if mod_time < oldest_mod_time or not found_file:
                found_file = True
                oldest_mod_time = mod_time
                oldest_path = path

    oldest_path = utilities.fullPath(oldest_path)
    if found_file:
        print("Deleting "+oldest_path+" from cache\n")
        os.remove(oldest_path)
        while oldest_path != "/" and oldest_path != "":
            try: os.rmdir(oldest_path)
            except OSError: pass
            finally: oldest_path = os.path.dirname(oldest_path)
        return True
    else:
        return False

def cacheCopy(src, dst, min_free, file_map, no_delete):
    #Cache a copy of src if possible, removing old files from cache if necessary

    src_size = os.stat(src).st_size * 2 #Safety factor of 2 to account for file growth if cached copy is modified

    du = os.statvfs(utilities.fullPath("/scratch/babymaker"))
    avail = du.f_bsize*du.f_bavail
    while avail-src_size < min_free:
        #Keep deleting until there's room
        if no_delete: return
        removed_file = removeOldCache(file_map)
        if not removed_file: return
        du = os.statvfs(utilities.fullPath("/scratch/babymaker"))
        avail = du.f_bsize*du.f_bavail
    print("Caching "+src+" to "+dst+"\n")
    shutil.copy(src, dst)
    os.chmod(dst, 0775)
    src_time = os.path.getmtime(src)
    dst_time = os.path.getmtime(dst)
    if dst_time < src_time:
        #Make sure cached copy is newer than /net copy
        touch(dst, (src_time+1., src_time+1.))

def execute(command, file_map, fragile):
    inv_file_map = dict()
    if not fragile:
        command, inv_file_map = mapFiles(command, file_map)
    else:
        inv_file_map = dict((cached,net) for net,cached in file_map.iteritems())

    if len(command) <= 0: return

    args = ["run/wrapper.sh"]
    for a in command:
        args.append(a.lstrip())
    command = args

    old_mod_times = dict()
    before_time = time.time()-1.
    # 1 second safety margin in case executable modifies within access
    # time resolution (typically 1 second)
    for f in inv_file_map:
        touch(f, (before_time, before_time))
        old_mod_times[f] = os.path.getmtime(f)

    print("Executing",command,"\n")
    utilities.flush()
    subprocess.call(command)
    utilities.flush()

    for f in inv_file_map:
        #Copy modified files back to /net
        if os.path.getmtime(f) > old_mod_times[f]:
            netCopy(f, inv_file_map[f])

def cacheRecurse(caches, file_map, command, fragile, min_free, no_delete):
    if len(caches)==0:
        #Caching done, run exectuable
        execute(command, file_map, fragile)
        return

    net_path = caches[0]
    if not isNetFile(net_path):
        utilities.ePrint("Cannot cache "+net_path+"\n")
        cacheRecurse(caches[1:], file_map, command, fragile, min_free, no_delete)
        return

    utilities.ensureDir(os.path.dirname(net_path))
    if not os.path.exists(net_path):
        #If /net file does not exist, create new empty file
        with open(net_path, "a"):
            pass
    cache_path = cachePath(net_path)
    utilities.ensureDir(os.path.dirname(cache_path))

    if not os.path.exists(cache_path) or os.path.getmtime(cache_path)<os.path.getmtime(net_path):
        #Cache doesn't exist or is outdated, so copy file from /net
        cacheCopy(net_path, cache_path, min_free, file_map, no_delete)

    if os.path.exists(cache_path) and os.path.getmtime(cache_path)>=os.path.getmtime(net_path):
        #Only use cached file if it was created and up-to-date
        file_map[net_path] = cache_path

    cacheRecurse(caches[1:], file_map, command, fragile, min_free, no_delete)

def cacheRun(caches, command, fragile, abs_limit, rel_limit, no_delete):
    if not os.path.isdir("/scratch/babymaker"):
        cacheRecurse([], dict(), command, True, 0, True)
        return
    caches = expand(caches)
    du = os.statvfs(utilities.fullPath("/scratch/babymaker"))
    min_free = max(abs_limit, du.f_bsize*du.f_blocks*rel_limit)
    cacheRecurse(caches, dict(), command, fragile, min_free, no_delete)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automatically creates and, if necessary, deletes local caches of files from /net/cmsX/cmsXr0/babymaker and remaps any files in the provided command to their cached version.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--cache", nargs="+", default=[],
                        help="Files that should be copied to cache or, if not already existing, created on the cache and moved to /net/cmsX/cmsXr0/babymaker upon completion.")
    parser.add_argument("-e", "--execute", nargs="+", default=[],
                        help="Command to execute. ./run/wrapper.sh is automatically prepended.")
    parser.add_argument("--fragile", action="store_true",
                        help="By default, wildcards are expanded and cached paths replaced in the arguments to the provided executable. Setting this flag runs the command \"as is.\" Files will still be cached, but the executable will not automatically use the cached version.")
    parser.add_argument("--abs_limit", default=10000000000, type=int,
                        help="Minimum number of bytes to leave available in cache.")
    parser.add_argument("--rel_limit", default=0.5, type=float,
                        help="Minimum fraction of cache to leave available.")
    parser.add_argument("--no_delete", action="store_true",
                        help="If cache is full, prevents deletion of old files to make room for new ones. Note that it is possible to delete a cached file currently being used by another script, so it is polite to use this flag if the cache is under heavy load.")
    args = parser.parse_args()

    cacheRun(args.cache, args.execute, args.fragile, args.abs_limit, args.rel_limit, args.no_delete)
