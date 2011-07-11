#!/bin/bash
# mpp standard wrapper to call packaged python scripts
# will take the name of the script as the wrapper name
dir=`dirname $(readlink -f $0)`/Owlcat

usage() {
    echo "usage: `basename $0` <cmd> [args]"
    allowed=`ls ${dir}/*.py`
    echo "        where cmd is one of:"
    for cmd in ${allowed} ; do
        echo -n "        "
        echo `basename ${cmd} | cut -f1 -d.`
    done
}

if [ $# -lt 1 ]; then
   usage
   exit 1;
fi

script=${dir}/$1.py
if [ ! -x ${script} ]; then
   echo unknown command $1
   usage
   exit 1;
fi
shift
exec $script $*
