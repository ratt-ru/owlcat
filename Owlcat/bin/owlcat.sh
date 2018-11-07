#!/bin/bash
# mpp standard wrapper to call packaged python scripts
# will take the name of the script as the wrapper name
#
#% $Id$ 
#
#
# Copyright (C) 2002-2011
# The MeqTree Foundation & 
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

dir=`dirname $(readlink -f $0)`
cmdfile=$dir/commands.list

if [ ! -f $cmdfile ]; then
  echo "$cmdfile not found"
  exit 1
fi

listCmds() {
  cmds="`cut -d " " -f 1 $cmdfile`"
  echo $cmds
}

usage() {
    if [ -f $dir/version.sh ]; then
      source $dir/version.sh
      echo "This is Owlcat version" $OWLCAT_VERSION
    fi
    echo "Usage: `basename $0` <command> [args]"
    echo "    Available commands are:"
    cat $cmdfile
    
    echo "    To get more help on a specific command, use:"
    echo "`basename $0` <command> -h"
    echo "    Additional options:"
    echo "-h         show this help message and exit"
    echo "-l         list available commands and exit"
}


if [ $# -lt 1 ]; then
   usage
   exit 1
fi

if [ "$1" == "-h" ]; then
   usage
   exit 0
fi

commands="`listCmds`"

if [ "$1" == "-l" ]; then
   echo $commands
   exit 0
fi

# make sure only valid commands from the list are executed
for cmd in $commands; do
  if [ "$1" == "$cmd" ]; then
    script=${dir}/$1.py
    if [ ! -x ${script} ]; then
      script=${dir}/$1.sh
      if [ ! -x ${script} ]; then
        echo "Unexpected error: missing script for command $1. Please report this as a bug."
        exit 1
      fi
    fi
    shift
    exec $script $*
  fi
done

echo "Unknown command $1"
usage
exit 1
