"""Pyxis module for MS-related operations""";
import numpy as np
import math
from pyrap.tables import table
import os.path
import pyfits

from Pyxis import *

# register ourselves with Pyxis, and define the superglobals
register_pyxis_module(superglobals="MS DDID FIELD IMAGE");


# external tools  
lwimager = x.lwimager;
tigger = xz.tigger

PLOTMS_ARGS = "";
plotms = x("plot-ms.py").args("$PLOTMS_ARGS ${-D <DDID} ${-F <FIELD}"); 

def msw (msname="$MS",subtable=None):
  """Opens the MS or a subtable read-write, returns table object"""
  return ms(msname,subtable,write=True);

def ms (msname="$MS",subtable=None,write=False):
  """Opens the MS or a subtable (read-only by default), returns table object"""
  msname = interpolate_locals("msname");
  if not msname:
    raise ValueError("'msname' or global MS variable must be set");
  if subtable:
    msname = os.path.join(msname,subtable);
  return table(msname,subtable,readonly=not write);

def _filename (base,newext):
  while base and base[-1] == "/":
    base = base[:-1];
  return os.path.splitext(base)[0]+"."+newext;
  
def uvcov (msname="$MS",save=None):
  """Makes uv-coverage plot"""
  msname,save = interpolate_locals("msname save");
  uv = ms(msname).getcol("UVW")[:,:2];
  import pylab
  pylab.plot(uv[:,0],uv[:,1],'.b');
  pylab.plot(-uv[:,0],-uv[:,1],'.r');
  pylab.savefig(save) if save else pylab.show();

def _strtodeg (arg):
  """Helper function, converts value+unit argument to degrees."""
  if isinstance(arg,str):
    for suffix,scale in ("'",60),("arcmin",60),('"',3600),("arcsec",3600),("rad",180/math.pi):
      if arg.endswith(suffix):
        return float(size[:-len(suffix)])*scale;
  return float(arg);

def dirtyimage (msname="$MS",output="${msname:BASE}.dirty.fits",show=True,size=None,**kw):
  """Makes a dirty image""";
  msname,output = interpolate_locals("msname output");
  if not msname:
    raise ValueError("'msname' or global MS variable must be set");
  # set output
  if not output:
    raise ValueError("'output' not set");
  if os.path.exists(output):
    os.remove(output);
  kw['fits'] = IMAGE = output;
  
  # set size, if specified
  if size is not None:
    # assume degrees, unless suffix is given
    degs = _strtodeg(size);
    # if npix is given, set cellsize
    if 'npix' in kw:
      kw['cellsize'] = cellsize = "%farcsec"%(degs*3600./kw['npix']);
      verbose(2,"specified size $size, set cellsize=$cellsize");
    elif 'cellsize' in kw:
      kw['npix'] = npix = round(degs/_strtodeg(kw['cellsize']));
      verbose(2,"specified size $size, set npix=$npix");
      
  lwimager(ms=msname,**kw);

  verbose(1,"Setting IMAGE=%s"%kw['fits']);
  if show:
    tigger(IMAGE);
  
def imgplot (image="$IMAGE",range=None,index=None,save=None):
  """Plots image""";
  image,save = interpolate_locals("image save");
  if not image:
    raise ValueError("'image' or global IMAGE variable must be set");
  # get FITS data
  dd = pyfits.open(image)[0].data.transpose();
  slicer = [slice(None),slice(None)] + (list(index) if index else [0]*(dd.ndim-2));
  dd = dd[tuple(slicer)]
  import pylab
  pylab.imshow(dd);
  if range:
    pylab.clim(*range);
  pylab.savefig(save) if save else pylab.show();

