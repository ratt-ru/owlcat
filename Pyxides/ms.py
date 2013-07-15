"""Pyxis module for MS-related operations""";
import math
import pyrap.tables
from pyrap.tables import table
import os.path
import pyfits

from Pyxis.ModSupport import *

import std

# register ourselves with Pyxis, and define the superglobals
register_pyxis_module();

v.define("MS","",
  """current measurement set""");
  
define("DDID",0,
  """current DATA_DESC_ID value""");
define("FIELD",0,
  """current FIELD value""");
define("IFRS","all",
  """interferometer subset""");
define("CHANRANGE",None,
  """channel range, as first,last[,step], or list of such tuples per DDID, or None for all""");
define("MS_TDL_Template",'ms_sel.msname=$MS ms_sel.ddid_index=$DDID ms_sel.field_index=$FIELD',
  """these options get passed to TDL scripts to specify an MS""");
  

lwimager = x.lwimager;
flagms          = x("flag-ms.py").args("$MS ${-I <IFRS} $CHAN_OWLCAT");
downweigh_redundant = x("downweigh-redundant-baselines.py").args("$MS ${-I <IFRS}");
aoflagger       = x.aoflagger.args("$MS");
addbitflagcol   = x.addbitflagcol.args("$MS");
wsrt_j2convert  = x.wsrt_j2convert.args("in=$MS");
mergems         = x("merge-ms.py");


define("PLOTVIS","CORRECTED_DATA:I",
  """passed to plot-ms to plot output visibilities. Set to None to skip plots.""");
define("PLOTMS_ARGS","",
  """extra plot-ms arguments""");
plotms = x("plot-ms.py").args("$MS $PLOTVIS $PLOTMS_ARGS ${-D <DDID} ${-F <FIELD} ${-I <IFRS} $CHAN_OWLCAT");


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


def prep (msname=None):
  if msname:
    v.MS = msname;
  info("adding imaging columns to $MS");
  pyrap.tables.addImagingColumns(v.MS);
  info("adding bitflag column");
  addbitflagcol();
  info("copying FLAG to bitflag 'legacy'");
  flagms("-Y +L -f legacy -c");
  
  
def copycol (msname="$MS",fromcol="DATA",tocol="CORRECTED_DATA"):
  """Copies data from one column of MS to another""";
  msname,fromcol,tocol = interpolate_locals("msname fromcol tocol");
  info("Copying $msname $fromcol to $tocol");
  tab = msw(msname)
  tab.putcol(tocol,tab.getcol(fromcol))
  tab.close()
  


def uvcov (msname="$MS",save=None):
  """Makes uv-coverage plot
  'msname' is superglobal MS by default.
  If 'save' is given, saves figure to file.
  """
  msname,save = interpolate_locals("msname save");
  uv = ms(msname).getcol("UVW")[:,:2];
  import pylab
  pylab.plot(uv[:,0],uv[:,1],'.b');
  pylab.plot(-uv[:,0],-uv[:,1],'.r');
  pylab.savefig(save) if save else pylab.show();

def rebin_freq (msname="$MS",output="$MSOUT",step=1):
  msname,output = interpolate_locals("msname output");
  std.runcasapy("""ms.open("$msname"); ms.split(outputms='$output',step=$step);""");

def _strtodeg (arg):
  """Helper function, converts value+unit argument to degrees."""
  if isinstance(arg,str):
    for suffix,scale in ("'",60),("arcmin",60),('"',3600),("arcsec",3600),("rad",180/math.pi):
      if arg.endswith(suffix):
        return float(size[:-len(suffix)])*scale;
  return float(arg);

#def dirtyimage (msname="$MS",output="${msname:BASE}.dirty.fits",show=True,size=None,**kw):
  #"""Makes a dirty image.
  #'msname' is superglobal MS by default.
  #'output' is output fits file, built from the MS name by default
  #If 'show' is True, invokes tigger on the file if available.
  
  #""";
  #msname,output = interpolate_locals("msname output");
  #if not msname:
    #raise ValueError("'msname' or global MS variable must be set");
  ## set output
  #if not output:
    #raise ValueError("'output' not set");
  #if os.path.exists(output):
    #os.remove(output);
  #kw['fits'] = IMAGE = output;
  
  ## set size, if specified
  #if size is not None:
    ## assume degrees, unless suffix is given
    #degs = _strtodeg(size);
    ## if npix is given, set cellsize
    #if 'npix' in kw:
      #kw['cellsize'] = cellsize = "%farcsec"%(degs*3600./kw['npix']);
      #verbose(2,"specified size $size, set cellsize=$cellsize");
    #elif 'cellsize' in kw:
      #kw['npix'] = npix = round(degs/_strtodeg(kw['cellsize']));
      #verbose(2,"specified size $size, set npix=$npix");
      
  #lwimager(ms=msname,**kw);

  #verbose(1,"Setting IMAGE=%s"%kw['fits']);
  #if show:
    #tigger(IMAGE);
  
#def imgplot (image="$IMAGE",range=None,index=None,save=None):
  #"""Plots image""";
  #image,save = interpolate_locals("image save");
  #if not image:
    #raise ValueError("'image' or global IMAGE variable must be set");
  ## get FITS data
  #dd = pyfits.open(image)[0].data.transpose();
  #slicer = [slice(None),slice(None)] + (list(index) if index else [0]*(dd.ndim-2));
  #dd = dd[tuple(slicer)]
  #import pylab
  #pylab.imshow(dd);
  #if range:
    #pylab.clim(*range);
  #pylab.savefig(save) if save else pylab.show();

## current spwid and number of channels. Note that these are set automatically from the MS by the _msddid_Template below
SPWID = 0
TOTAL_CHANNELS = 0

## whenever the channel range changes, setup strings for TDL & Owlcat channel selection (CHAN_TDL and CHAN_OWLCAT),
## and also CHANSTART,CHANSTEP,NUMCHANS
_chanspec = None;
def _chanspec_Template ():
  global CHAN_TDL,CHAN_OWLCAT,CHANSTART,CHANSTEP,NUMCHANS;
  chans = CHANRANGE;
  if isinstance(CHANRANGE,(list,tuple)) and type(CHANRANGE[0]) is not int:
    chans = CHANRANGE[DDID];
  # process channel specification 
  if chans is None:
    CHAN_OWLCAT = '';
    CHANSTART,CHANSTEP,NUMCHANS = 0,1,TOTAL_CHANNELS;
    CHAN_TDL = 'ms_sel.select_channels=0';
  else:
    if type(chans) is int:
      ch0,ch1,dch = chans,chans,1;
#      CHANSTART,CHANSTEP,NUMCHANS = chans,1,1;
    elif len(chans) == 1:
      ch0,ch1,dch = chans[0],chans[0],1;
#      CHANSTART,CHANSTEP,NUMCHANS = chans[0],1,1;
    elif len(chans) == 2:
      ch0,ch1,dch = chans[0],chans[1],1;
#      CHANSTART,CHANSTEP,NUMCHANS = chans[0],1,chans[1]-chans[0]+1;
    elif len(chans) == 3:
      ch0,ch1,dch = chans;
    CHANSTART,CHANSTEP,NUMCHANS = ch0,dch,((ch1-ch0)//dch+1);
    CHAN_OWLCAT = "-L %d~%d:%d"%(ch0,ch1,dch);
    CHAN_TDL = 'ms_sel.select_channels=1 ms_sel.ms_channel_start=%d ms_sel.ms_channel_end=%d ms_sel.ms_channel_step=%d'%\
               (ch0,ch1,dch);
  return CHANSTART,CHANSTEP,NUMCHANS;
