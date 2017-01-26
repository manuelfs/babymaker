#! /bin/bash

## Writing a tag because different code is need in 7.4.12 and before
if [ ! -e bmaker/interface/release.hh ]
then
    if [[ $CMSSW_BASE == *"CMSSW_7_4"* ]]
    then
	printf "#define PRE_8_0\n\n" > bmaker/interface/release.hh

	if [[ $CMSSW_BASE == *"CMSSW_7_4_6"* ]]
	then
	    printf "#define PRE_7_4_12\n\n" > bmaker/interface/release.hh
	else
	    printf "#define POST_7_4_12\n\n" > bmaker/interface/release.hh
	fi
    else
	printf "#define POST_7_4\n\n" > bmaker/interface/release.hh
    fi
fi

# ## Checking if the event filter list file exist, and if not, download it
# if [ ! -e data/csc_beamhalo_filter/csc2015_ee4sc_Jan13.txt ]
# then
#     echo "Donwloading filter list file"
#     cd data/csc_beamhalo_filter/
#     wget http://hep.ucsb.edu/people/manuelf/ra4/csc2015_ee4sc_Jan13.txt
#     cd ../..
# fi

## Compiling genfiles which generates baby_base and associated babies
cd bmaker/genfiles
exit_code=0;
if [ $# -ne 0 ] && [ "$1" == "clean" ]
then
    rm -rf run/*.exe bin/*.o bin/*.a bin/*.d *.exe *.out
    ./run/remove_backups.sh
    exit_code=$?
else
    tmp_file=mktemp
    
    make -j $(getconf _NPROCESSORS_ONLN) -k -r -R 2> >(tee $tmp_file >&2)
    exit_code=$?
    
    echo
    
    if [[ $exit_code != 0 ]] ; then
	echo "ERRORS AND WARNINGS:"
	cat $tmp_file >&2
    else
	echo "Generated baby header and source code files :o)"
	echo
    fi
    
    rm -rf $tmp_file
fi

cd ../../..  ## Going back to the base babymaker folder

## unpacking JECs
count=`ls -1 *.txt 2>/dev/null | wc -l`
if [ $count == 0 ]
then 
    cd babymaker/data/jec/
    for i in $( ls *.tar.gz); do
        tar --strip-components=2 -zxf $i
    done
    # the FastSim does not have nested folders
    tar -zxf Spring16_FastSimV1_MC.tar.gz    
    ls -1 *.txt | grep -v 'AK4PFchs.txt' | xargs rm
    cd -
fi 

## Scramming CMSSW code
if [ $# -ne 0 ] && [ "$1" == "clean" ]
then
    scram b clean
else
    scram b -j$(getconf _NPROCESSORS_ONLN)
fi

cd babymaker

exit $exit_code
