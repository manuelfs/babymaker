# `genfiles`

Code for processing existing babies without CMSSW dependencies.

## Code setup

The `babymaker/bmaker/genfiles` directory can be compiled as a
standalone set of scripts independent from the rest of the `babymaker`
repository. To do so, simply issue the following commands from this
directory:

    cmsenv
    ./compile.sh

Note that running scram and compiling in the root `babymaker`
directory are unnecessary.

## Post-processing of ntuple production

Post-processing relies on the UCSB tier 3, so all ntuples must first
be copied to a folder like `/net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/mc/to_renormalize/`. All commands
are issued from the `genfiles` directory for both MC and data.

### Data

1. Combine datasets removing duplicates using 

        python/send_combine_ntuples.py

    where you'll need to set the proper `infolder`, `outfolder`, `datasets`, and `jsonfile`. This script sends
    combination jobs for groups of `run_files` runs.

2. Skim the combined dataset. Each skim requires one execution of
`python/send_skim_ntuples.py`. Suppose, for example, that one want to
produce the standard and baseline skims for the "alldata" combined
dataset produced in the last step; then one would issue the commands

        python/send_skim_ntuples.py /net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/data/unskimmed/alldata standard
        python/send_skim_ntuples.py /net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/data/unskimmed/alldata baseline

    This will produce the subdirectories `skim_standard` and `skim_baseline` inside of `data`.

3. Slim and merge the skimmed ntuples. Multiple skims can be slimmed
with multiple slimming rules using `python/send_slim_ntuples.py`. The
script takes the outer product of the requested slims and skims, so if
one wanted the minimal and base_data slims for both the standard and
baseline skims produced in the last step, one would issue the command

        python/send_slim_ntuples.py --input_dir /net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/data/alldata --skims standard baseline --slims minimal base_data

    This will produce subdirectories `slim_minimal_standard`,
`slim_minimal_baseline`, `slim_base_data_standard`, and
`slim_base_data_baseline` inside `alldata`.

4. Validate! More extensive validation scripts are in progress, but
one should minimally check that the number of root files stays the
same after skimming and that the number of events stays the same after
slimming and merging. The number of files can be checked with the
command

        python/count_root_files.py -f /net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/data

### MC

1. If one is processing signal scans, separate the mass points. This
is done with the script below after setting the proper `infolder` and `outfolder` and `outname`

        ./run/send_split_scan.py 


2. Renormalize the weights so that the cross section is kept constant. Set up the correct `infolder` and `outfolder` and run

        python/send_change_weights.py

3. Validate `weight` by setting the proper `infolder`, `outfolder` in 

        ./python/validate_ntuples.py

    and running it. This script compares the sum of weights to an older production. The variable `weight`
    might have changed, so for instance when comparing Marmot and Capybara, the variables that would agree
    would be `weight/w_toppt/eff_trig` and `weight/w_isr/w_pu`

4. Follow steps 2 through 4 from [Data](#Data).

## Utility scripts and data caching

There are two useful tools for running on the batch system:
`run/wrapper.sh` and `python/cache.py`. The former should be inserted
between `JobSubmit.csh` and the issued command for any batch job to
ensure the environment variables are set correctly during execution on
the Tier 3. For example, one might send

    JobSubmit.csh run/wrapper.sh echo "Hello world"

The latter script, `python/cache.py`, is a bit more general. It has
two main options, `--cache` and `--execute.` The `--execute` option is
followed by a command to be executed, much like `run/wrapper.sh`
is. The `--cache` option is followed by an arbitrary number of file
paths (possibly containing wildcards) to be cached in
`/scratch/babymaker`. If the provided file path already exists, it is
assumed to be an input file and simply copied to `/scratch/babymaker`;
if it does not exists, a temporary file is created in
`/scratch/babymaker` and then copied to the provided path when the
command given to `--execute` is done. Note that if any file in the
`--cache` list or any file appearing in the command line arguments for
`--execute` is modified by the executed script, the modified file is
copied back from the cache to the original location after
completion. As a simple example, it can be used as a replacement for
`run/wrapper.sh`:

    JobSubmit.csh python/cache.py --execute echo "Hello world"

The caching utility is mainly intended for skimming (it is used by
`python/send_skim_ntuples.py`) and other processing which performs a
large number of small write operations. For example, one can skim a
single file with

    python/cache.py --cache output_file.root --execute python/skim_ntuple.py baseline output_file.root /net/cmsX/cmsXr0/babymaker/babies/YYYY_MM_DD/data/unskimmed/*.root " --keep-existing"

There are two things one should observe in this example. First, the
file `output_file.root` can be used as an argument passed to
`python/skim_ntuple.py`; the caching script knows to automatically
replace it with the temporary file path created in
`/scratch/babymaker`. If one does not want this behavior, simply add a
`--fragile` and file paths will be left intact. Second, the
`--keep-existing` argument for `python/skim_ntuple.py` is quoted and
has a space before the dashes. This is necessary to prevent
python/cache.py from thinking the option was intended for itself and
instead have it forward the option to `python/skim_ntuple.py`.

The `python/cache.py` script does some crude management of the
available disk space. To prevent `/scratch/babymaker` from becoming
full, it accepts options specifying an absolute (`--abs_limit`) and
relative (`--rel_limit`) amount of disk space to leave empty. If the
cache is becoming too full, it will delete the least recently used
files in the cache until it has space for the new files. In the
extreme scenario where the cache is nearly full with recently used
files and being access by multiple processes, this can result in one
process deleting a cached file that is being read by
another. Hopefully this occurs very rarely, but the option
`--no-delete` is available to prevent files from being deleted from
the cache.
