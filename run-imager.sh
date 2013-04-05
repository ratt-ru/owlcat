#!/bin/bash
# -*- coding: utf-8 -*-
# if config doesn't exist, strange things happen!

# if dirty image is not made first, clean fails!
# Same message as when ddid_index is not set properly:
#               Maximum of approximate PSF for field 1 = 0 : renormalizing to unity
# But not sure what state persists between making a dirty image and a clean one

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

CONFFILE=imager.conf

# Unset all variables beginning with img_


# Default settings
img_lwimager=lwimager
img_image2fits=image2fits
img_data=DATA
img_ifrs=""
img_weight=natural
img_taper=""
img_robust=0
# wfov/wnpix is optional, only works with newer imagers, so only pass through if present in config files
#img_wfov=0rad
#img_wnpix=0
img_spwid=0
img_field=0
img_size=
img_npix=512
img_cellsize=1arcsec
img_mode=channel
img_stokes=I
img_padding=1.2
img_flux_scale=2
img_phasecenter=""
img_wprojplanes=0
img_select=""
img_uvmin=""
img_uvmax=""
img_mask=""

img_channels=""
img_img_channels=""

img_prefervelocity=False
img_remove_img=1

img_name_dirty=dirty
img_name_restored=restored
img_name_model=model
img_name_residual=residual
img_name_prefix=""
img_name_suffix=""

img_cachesize=512

img_oper=image

img_niter=1000
img_gain=.1
img_threshold=0Jy
img_fixed=0

img_export_column=0

# do we have a config file on the command line?
for arg in $*; do
  if [ "${arg%.conf}" != "$arg" -a "${arg#+}" == "$arg" -a -f $arg ]; then
    CONFFILE=$arg
    echo "Using imager config file $CONFFILE"
    break
  fi
done

# load conf file
if [ -f $CONFFILE ]; then
  source $CONFFILE
fi

# load extra config files
for arg in $*; do
  if [ "${arg%.conf}" != "$arg" -a "${arg#+}" != "$arg" -a -f ${arg#+} ]; then
    echo "Loading additional config file ${arg#+}"
    source ${arg#+}
  fi
done


trial=true
confirm=true

echo "run-imager.sh command line: $*"

# overwrite with command line
while [ "$1" != "" ]; do
  if [ "${1%=*}" != "$1" ]; then
    eval "img_$1"
  elif [ "$1" == "-h" -o "$1" == "--help" ]; then
    echo "`basename $0`: wrapper script for lwimager"
    echo "Imager parameters are read from $CONFFILE, and can also be overridden"
    echo "on the command line as parm=value."
    echo "See `dirname $0`/imager.conf.example for an example imager.conf"
    echo "Recognized parameters (and their current settings) are: "
    set | egrep \^img_ | cut -c 5-
    exit 0
  elif [ "$1" == "-nr" ]; then
    img_remove_img=0
  elif [ "$1" == "-trial" ]; then
    trial="echo 'Trial mode only, exiting'; exit 1";
  elif [ "$1" == "-i" ]; then
    confirm="echo 'Press Enter to continue with these settings, or Ctrl+C to stop'; read";
  fi
  shift
done

# expand vars in filenames
img_name_dirty=`eval echo $img_name_prefix$img_name_dirty$img_name_suffix`
img_name_restored=`eval echo $img_name_prefix$img_name_restored$img_name_suffix`
img_name_model=`eval echo $img_name_prefix$img_name_model$img_name_suffix`
img_name_residual=`eval echo $img_name_prefix$img_name_residual$img_name_suffix`

# MS must be set
if [ -z "$img_ms" ]; then
  echo "MS must be set with ms=name, or in $CONFFILE";
  exit 1
fi

# for weight=default, disable taper (otherwise the taper may be e.g. applied twice)
if [ "$img_weight" == "default" -a "$img_taper" != "" ]; then
  echo "WARNING: weight=default but taper=$img_taper is set. Disabling the taper!"
  img_taper=""
fi

# print final var settings
set | grep ^img_ | cut -c 5-

# expand size into arcmin and npix
if [ "$img_arcmin" != "" -a "$img_npix" != 0 ]; then
  img_cellsize=`python -c "print '%farcmin'%$img_arcmin/float($img_npix)"`
elif [ "$img_cellsize" == "" -o "$img_npix" == "" ]; then
  if [ "$img_size" == "" ]; then
    echo "ERROR: one of size, arcmin & npix, or cellsize & npix must be specified"
    exit 1
  fi
  local _arcmin=${img_size#*/}
  img_npix=${img_size%/*}
  img_cellsize=`python -c "print '%farcmin'%$_arcmin/float($img_npix)"`
fi

# print baseline settings
if [ "$img_uvmin" != "" ]; then
  if [ "${img_uvmin%,*}" == "$img_uvmin" ]; then
    img_select="${imgselect:+($img_select)&&}(SQRT(SUMSQUARE(UVW[1:2]))>=$img_uvmin)"
  else
    a=${img_uvmin%,*}
    b=${img_uvmin#*,}
    img_select="${imgselect:+($img_select)&&}((UVW[1]/$a)**2+(UVW[2]/$b)**2>=1)"
  fi
  echo "uvmin $img_uvmin, select string is $img_select"
fi
if [ "$img_uvmax" != "" ]; then
  img_select="${imgselect:+($img_select)&&}(SQRT(SUMSQUARE(UVW[1:2]))<=$img_uvmax)"
  echo "uvmax $img_uvmax, select string is $img_select"
fi

# print ifr settings
if [ "$img_ifrs" == "" ]; then
  echo "Using all interferometers"
else
  # find the ifr-subset-to-taql.sh utility
  dir=`dirname $(readlink -f $0)`
  subset2taql="$dir/ifr-subset-to-taql.sh"
  if [ ! -x $subset2taql ]; then
    subset2taql=`which ifr-subset-to-taql.sh`
    if [ ! -x $subset2taql ]; then
      echo "Can't find ifr-subset-to-taql.sh script"
      exit 1
    fi
  fi
  # run it to get the IFR selection
  ifr_select="`$subset2taql $img_ms ${img_ifrs// /,}`"
  echo "NB: IFR selection '$img_ifrs' applied, resulting in ${ifr_select%//*} ifrs"
  ifr_select="${ifr_select#*//}"
  if [ "$img_select" != "" ]; then
    img_select="($img_select)&&($ifr_select)"
  else
    img_select="$ifr_select"
  fi
  echo "select string is $img_select"
fi

eval $confirm

# kludge: export data before running
if [ "$img_export_column" == "$img_weight" ]; then
  export-ms-column.py $img_ms $img_data $img_name_dirty.expcol.gz
fi

# helper function: protects / and - in filenames with backslash,
# so image2fits does not interpret them as LEL arithmetic
sanitize ()
{
  name="$1"
  name="${name//\//\\/}"
  name="${name//-/\\-}"
  echo $name
}

# make dirty image
make_image ()
{
  # collect mandatory arguments
  cmd="$img_lwimager ms=$img_ms data=$img_data operation=$img_oper
      stokes=$img_stokes mode=$img_mode weight=$img_weight wprojplanes=$img_wprojplanes
      npix=$img_npix cellsize=$img_cellsize ${img_wfov:+wfov=$img_wfov} ${img_wnpix:+wnpix=$img_wnpix}
      spwid=$img_spwid field=$img_field padding=$img_padding cachesize=$img_cachesize
      prefervelocity=$img_prefervelocity
  "

  # add optional arguments
  if [ "$img_robust" != "" ]; then
    cmd="$cmd robust=$img_robust"
  fi
  if [ "$img_taper" != "" ]; then
    cmd="$cmd filter=${img_taper}arcsec,${img_taper}arcsec,0.000000deg"
  fi
  if [ "$img_channels" != "" ]; then
    declare -a chans=(${img_channels//,/ })
    cmd="$cmd chanmode=channel nchan=${chans[0]} chanstart=${chans[1]} chanstep=${chans[2]}";
  fi
  if [ "$img_img_channels" != "" ]; then
    declare -a chans=(${img_img_channels//,/ })
    cmd="$cmd img_nchan=${chans[0]} img_chanstart=${chans[1]} img_chanstep=${chans[2]}";
  fi
  if [ "$img_select" != "" ]; then
    cmd="$cmd select='$img_select'"
  fi
  if [ "$img_phasecenter" != "" ]; then
    cmd="$cmd phasecenter='$img_phasecenter'"
  fi
  if [ "$img_sigma" != "" ]; then
    cmd="$cmd sigma='$img_sigma'"
  fi
  if [ "$img_sigma" != "" ]; then
    cmd="$cmd targetflux='$img_targetflux'"
  fi

  # setup output files
  if [ "$img_oper" == "image" ]; then
    imgname="$img_name_dirty.img"
    imgname_fits="$img_name_dirty.fits"
    cmd="$cmd image=$imgname"
    rm -fr $imgname_fits $imgname &>/dev/null
    echo "Making dirty image: " $cmd
  else
    restored="$img_name_restored.img";
    model="$img_name_model.img";
    residual="$img_name_residual.img";
    restored_fits="$img_name_restored.fits";
    model_fits="$img_name_model.fits";
    residual_fits="$img_name_residual.fits";
    rm -fr $restored_fits $model_fits $residual_fits $restored $model $residual &>/dev/null
    cmd="$cmd model=$model restored=$restored residual=$residual
      niter=$img_niter gain=$img_gain threshold=$img_threshold fixed=$img_fixed
    "
    if [ "$img_maskblc" != "" -o "$img_masktrc" ]; then
      cmd="$cmd mask=${img_mask:-mask} maskblc=$img_maskblc masktrc=$img_masktrc maskvalue=1"
    fi
    echo "Making clean image: " $cmd
  fi

  eval $trial
  if eval time $cmd; then
    if [ "$img_remove_img" != "0" ]; then
      echo "Success, all *.img files will be removed after conversion to FITS"
      remove=rm
    else
      echo "Success"
      remove=true
    fi
    # call sanitize() (see above) to make sure / and - in filenames is not interpreted
    # as LEL expressions
    if [ "$img_oper" == "image" ]; then
      ${img_image2fits} in="`sanitize $imgname`"*$img_flux_scale out=$imgname_fits && $remove -fr $imgname
    else
      ${img_image2fits} in="`sanitize $model`"*$img_flux_scale out=$model_fits && $remove -fr $model
      ${img_image2fits} in="`sanitize $residual`"*$img_flux_scale out=$residual_fits && $remove -fr $residual
      ${img_image2fits} in="`sanitize $restored`"*$img_flux_scale out=$restored_fits && $remove -fr $restored
    fi
    return 0
  else
    echo "Imager failed";
    return 1
  fi
}

make_image
