import numpy as np
import math
from pyrap.tables import table
import os.path
import pyfits

MS = None;
IMAGE = None;

## template: strips off MS suffix to make MSBASE
#def MSBASE_Template ():
  #base = MS;
  #while base and base[-1] == "/":
    #base = base[:-1];
  #return os.path.splitext(base)[0];


# external tools  
lwimager = x.lwimager;
tigger = xz.tigger

def msw (msname="$MS"):
  """ms([msname]) Opens the MS (by argument, or global $MS variable) read-write, returns table object"""
  return ms(msname,write=True);

def ms (msname="$MS",write=False):
  """ms([msname]) Opens the MS (by argument, or global $MS variable) read-only, returns table object"""
  msname = interpolate_locals("msname");
  if not msname:
    raise ValueError("'msname' or global MS variable must be set");
  return table(msname,readonly=not write);

def _filename (base,newext):
  while base and base[-1] == "/":
    base = base[:-1];
  return os.path.splitext(base)[0]+"."+newext;
  
def ms_uvcov (msname="$MS",save=None):
  """Plots uv coverage of the specified MS, or the global 'MS'"""
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
  
def ms_dirtyimage (msname="$MS",output="$msname_BASE.dirty.fits",show=True,size=None,**kw):
  """Makes a dirty image from this MS""";
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
  
def img_plot (image="$IMAGE",range=None,index=None,save=None):
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
