# -*- coding: utf-8 -*-
# if config doesn't exist, strange things happen!

# if dirty image is not made first, clean fails! 
# Same message as when ddid_index is not set properly:
#               Maximum of approximate PSF for field 1 = 0 : renormalizing to unity
# But not sure what state persists between making a dirty image and a clean one

CONFFILE=imager.conf

# Unset all variables beginning with img_


# Default settings
img_lwimager=lwimager
img_data=DATA
img_ifrs="*"
img_weight=natural
img_spwid=0
img_field=0
img_size=512/60
img_mode=channel
img_stokes=I
img_padding=1.2
img_flux_scale=2
img_phasecenter=""

img_channels=""
img_img_channels=""

img_remove_img=1

img_name_dirty=dirty
img_name_restored=restored
img_name_model=model
img_name_residual=residual

img_cachesize=512

img_oper=image

img_niter=1000
img_gain=.1
img_threshold=0Jy

# load conf file
if [ -f $CONFFILE ]; then
  source $CONFFILE
fi

trial=true
confirm=true

# overwrite with command line
while [ "$1" != "" ]; do
  if [ "${1%=*}" != "$1" ]; then
    eval "img_$1"
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
img_name_dirty=`eval echo $img_name_dirty`
img_name_restored=`eval echo $img_name_restored`
img_name_model=`eval echo $img_name_model`
img_name_residual=`eval echo $img_name_residual`

# MS must be set
if [ -z "$img_ms" ]; then
  echo "MS must be set with MS=name, or in $CONFFILE";
fi

# print final var settings
set | grep ^img_ | cut -c 5-

# expand size into arcmin and npix
if [ "$img_arcmin" == "" -o "$img_npix" == "" ]; then
  img_arcmin=${img_size#*/}
  img_npix=${img_size%/*}
fi

# print ifr settings
if [ "$img_ifrs" == "" ]; then
  echo "Using all interferometers"
else
  ifr_select="`ifr-subset-to-taql.sh $img_ms ${img_ifrs// /,}`"
  echo "NB: IFR selection '$img_ifrs' applied, resulting in ${ifr_select%//*} ifrs"
  ifr_select="${ifr_select#*//}"
  if [ "$img_select" != "" ]; then
    img_select="($img_select)&&($ifr_select)"
  else
    img_select="$ifr_select"
  fi
  echo "select=$img_select"
fi

eval $confirm

CELLSIZE=`python -c "print $img_arcmin*60./float($img_npix)"`

# make dirty image
make_image ()
{
  # collect mandatory arguments
  cmd="$img_lwimager ms=$img_ms data=$img_data operation=$img_oper 
      stokes=$img_stokes mode=$img_mode weight=$img_weight
      npix=$img_npix cellsize=${CELLSIZE}arcsec 
      spwid=$img_spwid field=$img_field padding=$img_padding cachesize=$img_cachesize
  "

  # add optional arguments
  if [ "$img_weight" == "radial" -a "$img_taper" != "" ]; then
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
      niter=$img_niter gain=$img_gain threshold=$img_threshold
    "
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
    # the ridiculous-looking pattern below replaces "/" in imgname/model/residual/restored
    # with "\/", so that image2fits does not interpret them as a divide operation
    if [ "$img_oper" == "image" ]; then
      image2fits in=$img_flux_scale\*\\${imgname//\//\\\/} out=$imgname_fits && $remove -fr $imgname
    else
      image2fits in=$img_flux_scale\*\\${model//\//\\\/} out=$model_fits && $remove -fr $model
      image2fits in=$img_flux_scale\*\\${residual//\//\\\/} out=$residual_fits && $remove -fr $residual
      image2fits in=$img_flux_scale\*\\${restored//\//\\\/} out=$restored_fits && $remove -fr $restored
    fi
    return 0
  else
    echo "Imager failed";
    return 1
  fi
}

make_image
