"""Pyxis module for MS-related operations""";
from Pyxis import *
import pyrap.images
import os
import subprocess

# Pyxides.mq is needed to make the Cattery available (see reference to Meow below)
import mqt

# register ourselves with Pyxis and define the superglobals
register_pyxis_module(superglobals="MS DDID FIELD");


# external tools  
lwimager = x.lwimager;
rm_fr = x.rm.args("-fr");

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
_lwimager_args = set(("spwid field prior image model restored residual data mode filter nscales weight noise robust wprojplanes padding "+
    "cachesize stokes nfacets npix cellsize phasecenter field spwid chanmode nchan chanstart chanstep img_nchan img_chanstart img_chanstep "+
    "select operation niter gain threshold targetflux sigma fixed constrainflux prefervelocity mask maskblc masktrc uservector maskvalue").split(" "));

def run (ms='$MS',ifrs=None,**kw):
  print ms,MS,v.MS
  # make dict of imager arguments that have been specified globally or locally
  args = dict([ (arg,globals()[arg]) for arg in _lwimager_args if arg in globals() and globals()[arg] is not None ]);
  args.update([ (arg,kw[arg]) for arg in _lwimager_args if arg in kw ]);
  # add spwid and field arguments
  args.setdefault('ms',ms);
  DDID is not None and args.setdefault('spwid',DDID);
  FIELD is not None and args.setdefault('field',FIELD);
  # have an IFR subset? Parse that too
  ifrs = ifrs or globals()['ifrs'];
  if ifrs and ifrs.lower() != "all":
    import Meow.IfrSet
    subset = Meow.IfrSet.from_ms(ms).subset(ifrs).taql_string();
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
