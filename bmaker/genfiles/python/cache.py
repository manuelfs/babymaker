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
import signal

import utilities

class SignalError(Exception):
    def __init__(self, signum, frame):
        self.signum = signum
        self.frame = frame
    def __str__(self):
        return "Caught signal number "+str(self.signum)

def signalHandler(signum, frame):
    raise SignalError(signum, frame)

def mkdirPermissions(path, mode):
    if not path or path == "/":
        return
    (head, tail) = os.path.split(path)
    mkdirPermissions(head, mode)
    try:
        os.mkdir(path)
        os.chmod(path, mode)
    except OSError:
        pass

def cacheUpToDate(cache_path, net_path):
    cache = None
    net = None
    try: cache = os.stat(cache_path)
    except OSError as e:
        if e.errno == 2: return False
        else: raise
    try: net = os.stat(net_path)
    except OSError as e:
        if e.errno == 2: return True
        else: raise
    return cache.st_mtime >= net.st_mtime and cache.st_size == net.st_size

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
    s = os.stat(path)
    return max(s.st_ctime, s.st_mtime, s.st_atime)

def mapFiles(command, file_map):
    #Replace executable arguments with cached equivalent

    expanded_args = []
    for arg in command:
        globbed = glob.glob(arg)
        if len(globbed) > 0:
            #Argument represents file(s)
            for f in globbed:
                expanded_args.append(utilities.fullPath(f))
        else:
            expanded_args.append(arg)

    command = []
    inv_file_map = dict((cached,net) for net,cached in file_map.iteritems())
    for arg in expanded_args:
        if arg in file_map and cacheUpToDate(file_map[arg], arg):
            #Check if generated cache for file
            command.append(file_map[arg])
        elif isNetFile(arg):
            #Check if pre-existing cache
            cache_path = cachePath(arg)
            if cacheUpToDate(cache_path, arg):
                command.append(cache_path)
                inv_file_map[cache_path] = arg
            else:
                command.append(arg)
        else:
            command.append(arg)

    return command, inv_file_map

def netCopy(src, dst):
    print("Copying "+src+" to "+dst+"\n")
    try:
        shutil.copy(src, dst)
        while not cacheUpToDate(src, dst):
            #Want cache to be newer so it's considered up to date
            now = time.time()
            os.utime(src, (now, now))
    except:
        try:
            os.remove(dst)
        finally:
            os.remove(src)
        utilities.ePrint("Failed to copy "+src+" to "+dst+"\n")
        raise

def removeOldCache(file_map):
    #Deletes oldest cached file
    found_file = False
    oldest_mod_time = 0
    oldest_path = ""
    for root, dirs, files in os.walk(utilities.fullPath("/scratch/babymaker")):
        for f in files:
            path = os.path.join(root, f)
            if path in file_map.itervalues(): continue
            mod_time = lastTime(path)
            if mod_time < oldest_mod_time or not found_file:
                found_file = True
                oldest_mod_time = mod_time
                oldest_path = path

    if time.time()-oldest_mod_time <= 86400.:
        #Don't delete files used in last 24 hours
        return False
    oldest_path = utilities.fullPath(oldest_path)
    if found_file:
        print("Deleting "+oldest_path+" from cache\n")
        try: os.remove(oldest_path)
        except: return False
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
    try:
        shutil.copy(src, dst)
        os.chmod(dst, 0775)
        while not cacheUpToDate(dst, src):
            now = time.time()
            os.utime(dst, (now, now))
    except:
        os.remove(dst)
        utilities.ePrint("Failed to cache "+src+" to "+dst+"\n")
        raise

def syncCache(net_path, cache_path):
    try:
        now = time.time()
        os.utime(cache_path, (now, now))
        cache_m_time = os.path.getmtime(cache_path)
        while not cacheUpToDate(cache_path, net_path):
            #Make sure cache is newer
            cache_m_time += 1.
            now = max(cache_m_time, time.time())
            os.utime(cache_path, (now, now))
    except:
        os.remove(cache_path)
        utilities.ePrint("Failed to sync cache times")
        raise

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

    try:
        old_mod_times = dict()
        before_time = round(time.time()-2.)
        # 2 second safety margin in case executable modifies within access
        # time resolution (typically 1 second)
        for f in inv_file_map.iterkeys():
            os.utime(f, (before_time, before_time))
            old_mod_times[f] = os.path.getmtime(f)

        exit_code = 0
        print("Executing",command,"\n")
        utilities.flush()
        try:
            exit_code = subprocess.call(command)
        except SignalError as e:
            if e.signum != signal.SIGCLD and e.signum != signal.SIGCHLD:
                raise e
        utilities.flush()

        if exit_code != 0:
            raise Exception("Executable returned non-zero exit code.")
    except:
        for f in inv_file_map.iterkeys():
            os.remove(f)
        utilities.ePrint("Failed to execute",command,"\n")
        raise
    else:
        for cache_path, net_path in inv_file_map.iteritems():
            if os.path.getmtime(cache_path) > old_mod_times[cache_path]:
                #Copy modified files back to /net
                netCopy(cache_path, net_path)
            else:
                syncCache(net_path, cache_path)

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

    mkdirPermissions(os.path.dirname(net_path), 0775)
    if not os.path.exists(net_path):
        #If /net file does not exist, create new empty file
        with open(net_path, "a"):
            pass
    cache_path = cachePath(net_path)
    mkdirPermissions(os.path.dirname(cache_path), 0775)

    if not cacheUpToDate(cache_path, net_path):
        #Cache doesn't exist or is outdated, so copy file from /net
        cacheCopy(net_path, cache_path, min_free, file_map, no_delete)

    if cacheUpToDate(cache_path, net_path):
        #Only use cached file if it was created and up-to-date
        file_map[net_path] = cache_path

    cacheRecurse(caches[1:], file_map, command, fragile, min_free, no_delete)

def cacheRun(caches, command, fragile, abs_limit, rel_limit, no_delete):
    for s in [sig for sig in dir(signal) if sig.startswith("SIG")
              and not sig.startswith("SIG_")
              and sig!="SIGKILL"
              and sig!="SIGSTOP"]:
        signum = getattr(signal, s)
        signal.signal(signum,signalHandler)

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
    parser.add_argument("--no_delete", action="store_false",
                        help="If cache is full, prevents deletion of old files to make room for new ones. Note that it is possible to delete a cached file currently being used by another script, so it is polite to use this flag if the cache is under heavy load.")
    args = parser.parse_args()

    cacheRun(args.cache, args.execute, args.fragile, args.abs_limit, args.rel_limit, args.no_delete)
