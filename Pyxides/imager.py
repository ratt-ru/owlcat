"""Pyxis module for MS-related operations""";
from Pyxis.ModSupport import *

import pyrap.images
import os
import subprocess

import mqt,ms,std
 
# register ourselves with Pyxis and define the superglobals
register_pyxis_module(superglobals="MS LSM DESTDIR");

v.define("LSM","lsm.lsm.html",
  """current local sky model""");
  

# external tools  
lwimager = x.time.args("lwimager");
rm_fr = x.rm.args("-fr");
tigger_restore = x("tigger-restore");

# standard imaging options 
ifrs=""
npix=2048
cellsize="8arcsec"
mode="channel"
stokes="IQUV"
weight="briggs"
wprojplanes=0
cachesize=4096
niter=1000
gain=.1
threshold=0

# rescale images by factor
flux_rescale=1

# use velocity rather than frequency
velocity = False;

# known lwimager args -- these will be passed from keywords
_fileargs = set("image model restored residual".split(" ")); 
_lwimager_args = set(("ms spwid field prior image model restored residual data mode filter nscales weight noise robust wprojplanes padding "+
    "cachesize stokes nfacets npix cellsize phasecenter field spwid chanmode nchan chanstart chanstep img_nchan img_chanstart img_chanstep "+
    "select operation niter gain threshold targetflux sigma fixed constrainflux prefervelocity mask maskblc masktrc uservector maskvalue").split(" "));

def _run (**kw):
  # make dict of imager arguments that have been specified globally or locally
  args = dict([ (arg,globals()[arg]) for arg in _lwimager_args if arg in globals() and globals()[arg] is not None ]);
  args.update([ (arg,kw[arg]) for arg in _lwimager_args if arg in kw ]);
  # add ifrs, spwid and field arguments
  ms.IFRS is not None and args.setdefault('ifrs',ms.IFRS);
  ms.DDID is not None and args.setdefault('spwid',ms.DDID);
  ms.FIELD is not None and args.setdefault('field',ms.FIELD);
  # have an IFR subset? Parse that too
  msname,ifrs = kw['ms'],args.pop('ifrs',None);
  if ifrs and ifrs.lower() != "all":
    import Meow.IfrSet
    subset = Meow.IfrSet.from_ms(msname).subset(ifrs).taql_string();
    args['select'] = "(%s)&&(%s)"%(args['select'],subset) if 'select' in args else subset;
  # make image names
  fitsfiles = {};
  for x in _fileargs:
    if x in args:
      fitsfiles[x] = args[x];
      args[x] = args[x]+".img";
  # run the imager
  lwimager(**args);
  # convert to FITS
  fs = kw.get('flux_rescale') or flux_rescale;
  velo = kw.get('velocity') or velocity;
  for x in _fileargs:
    if x in args:
      im = pyrap.images.image(args[x]);
      if fs and fs != 1:
        im.putdata(fs*im.getdata());
      im.tofits(fitsfiles[x],overwrite=True,velocity=velo);  
      subprocess.call("rm -fr "+args[x],shell=True);

# filenames for images
define("DIRTY_IMAGE_Template", "${OUTFILE}.dirty.fits","output filename for dirty image");
define("RESTORED_IMAGE_Template", "${OUTFILE}.restored.fits","output filename for restored image");
define("RESIDUAL_IMAGE_Template", "${OUTFILE}.residual.fits","output filename for deconvolution residuals");
define("MODEL_IMAGE_Template", "${OUTFILE}.model.fits","output filename for deconvolution model");
define("FULLREST_IMAGE_Template", "${OUTFILE}.fullrest.fits","output filename for LSM-restored image");

# How to channelize the output image. 0 for average all, 1 to include all, 2 to average with a step of 2, etc.
# None means defer to 'imager' module options
define("IMAGE_CHANNELIZE",0,"image channels selection: 0 for all, 1 for per-channel cube")
# passed to tigger-restore when restoring models into images. Use e.g. "-b 45" for a 45" restoring beam.
define("RESTORING_OPTIONS","","extra options to tigger-restore for LSM-restoring")
# default clean algorithm
define("CLEAN_ALGORITHM","clark","CLEAN algorithm (clark, hogbom, csclean, etc.)")


def make_image (msname="$MS",column="CORRECTED_DATA",
                dirty=True,restore=False,restore_lsm=True,
                channelize=None,lsm="$LSM",config="",**kw0):
  """Makes image(s) from MS. Set dirty and restore to True or False to make the appropriate images. You can also
  set either to a dict of options to be passed to the imager. If restore=True and restore_lsm is True and 'lsm' is set, 
  it will also make a full-restored image (i.e. will restore the LSM into the image) with tigger-restore. Use this when 
  deconvolving residual images. Note that RESTORING_OPTIONS are passed to tigger-restore.
  
  'config' specifies a config file for run-imager. If empty, the default imager.conf is used.
  
  'channelize', if set, overrides the IMAGE_CHANNELIZE setting. If both are None, the options in the 'imager' module take effect. 
  
  Image names are determined by the globals DIRTY_IMAGE, RESTORED_IMAGE, RESIDUAL_IMAGE, MODEL_IMAGE and FULLREST_IMAGE"""
  msname,column,lsm = interpolate_locals("msname column lsm"); 
  makedir(DESTDIR);
  
  if restore and column != "CORRECTED_DATA":
    abort("Due to imager limitations, restored images can only be made from the CORRECTED_DATA column.");
  
  # setup imager options
  kw0.update(dict(chanstart=ms.CHANSTART,chanstep=ms.CHANSTEP,nchan=ms.NUMCHANS));
  if channelize is None:
    channelize = IMAGE_CHANNELIZE;
  if channelize == 0:
    kw0.update(img_nchan=1,img_chanstart=ms.CHANSTART,img_chanstep=ms.NUMCHANS);
  elif channelize > 0:
    kw0.update(img_nchan=ms.NUMCHANS//channelize,img_chanstart=ms.CHANSTART,img_chanstep=channelize);
    
  kw0.update(ms=msname,data=column);

  if dirty:
    info("Making dirty image DIRTY_IMAGE=$DIRTY_IMAGE");
    kw = kw0.copy();
    if type(dirty) is dict:
      kw.update(dirty);
    _run(operation="image",image=DIRTY_IMAGE,**kw);
  if restore:
    info("Making restored image RESTORED_IMAGE=$RESTORED_IMAGE");
    info("       (MODEL_IMAGE=$MODEL_IMAGE RESIDUAL_IMAGE=$RESIDUAL_IMAGE)");
    kw = kw0.copy();
    if type(restore) is dict:
      kw.update(restore);
    _run(operation=CLEAN_ALGORITHM,restored=RESTORED_IMAGE,model=MODEL_IMAGE,residual=RESIDUAL_IMAGE,**kw)
    if lsm and restore_lsm:
      info("Restoring LSM into FULLREST_IMAGE=$FULLREST_IMAGE");
      tigger_restore("$RESTORING_OPTIONS","-f",RESTORED_IMAGE,lsm,FULLREST_IMAGE);
      
document_globals(make_image,"*_IMAGE IMAGE_CHANNELIZE MS RESTORING_OPTIONS CLEAN_ALGORITHM ms.IFRS ms.DDID ms.FIELD ms.CHANRANGE");      

      