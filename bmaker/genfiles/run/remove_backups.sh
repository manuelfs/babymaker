#! /bin/bash

currentdir=$(pwd)

cd ../..
rm -f bmaker/interface/baby*hh
rm -f bmaker/src/baby*cc

shopt -s nullglob

clean_recur(){
    cd $1
    rm -f *~ *#
    for i in $(ls -A)
    do
	if [ -d $i ]
	then
	    if [ ! -L $i ]
	    then
		clean_recur $i
	    fi
	fi
    done
    cd ..
}

clean_recur $(pwd)

cd $currentdir
