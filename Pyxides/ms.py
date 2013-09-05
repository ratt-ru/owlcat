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
_flagms  	= x("flag-ms.py");
flagms          = _flagms.args("$MS ${-I <IFRS} $CHAN_OWLCAT");
downweigh_redundant = x("downweigh-redundant-baselines.py").args("$MS ${-I <IFRS}");
_aoflagger      = x.aoflagger;
addbitflagcol   = x.addbitflagcol.args("$MS");
wsrt_j2convert  = x.wsrt_j2convert.args("in=$MS");
mergems         = x("merge-ms.py");
taql		= x("taql");

define("PLOTVIS","CORRECTED_DATA:I",
  """passed to plot-ms to plot output visibilities. Set to None to skip plots.""");
define("PLOTMS_ARGS","",
  """extra plot-ms arguments""");
plotms = x("plot-ms.py").args("$MS $PLOTVIS $PLOTMS_ARGS ${-D <DDID} ${-F <FIELD} ${-I <IFRS} $CHAN_OWLCAT");

class _TableGuard (object):
  def __init__ (self,tab):
    self.tab = tab;
  def __enter__ (self):
    return self.tab;
  def __exit__ (self):
    self.tab.close();

def msw (msname="$MS",subtable=None):
  """Opens the MS or a subtable read-write, returns table object."""
  return ms(msname,subtable,write=True);

def ms (msname="$MS",subtable=None,write=False):
  """Opens the MS or a subtable (read-only by default), returns table object."""
  msname = interpolate_locals("msname");
  if not msname:
    raise ValueError("'msname' or global MS variable must be set");
  if subtable:
    msname = table(msname).getkeyword(subtable);
  tab = table(msname,readonly=not write);
  return tab;

def _filename (base,newext):
  while base and base[-1] == "/":
    base = base[:-1];
  return os.path.splitext(base)[0]+"."+newext;


def prep (msname=None):
  """Prepares MS for use with MeqTrees: adds imaging columns, adds BITFLAG columns, copies current flags
  to 'legacy' flagset"""
  if msname:
    v.MS = msname;
  info("adding imaging columns to $MS");
  pyrap.tables.addImagingColumns(v.MS);
  info("adding bitflag column");
  addbitflagcol();
  info("copying FLAG to bitflag 'legacy'");
  flagms("-Y +L -f legacy -c");
  
  
  
def copycol (msname="$MS",fromcol="DATA",tocol="CORRECTED_DATA",rowchunk=500000):
  """Copies data from one column of MS to another.
  Copies 'rowchunk' rows at a time; decrease the default if you have low RAM.
  """;
  msname,fromcol,tocol = interpolate_locals("msname fromcol tocol");
  tab = msw(msname)
  nrows = tab.nrows();
  for row0 in range(0,nrows,rowchunk):
    nr = min(rowchunk,nrows-row0);
    info("Copying $msname $fromcol to $tocol (rows $row0 to %d)"%(row0+nr-1));
    tab.putcol(tocol,tab.getcol(fromcol,row0,nr),row0,nr)
  tab.close()

  
define('FIGURE_WIDTH',8,'width of plots, in inches');
define('FIGURE_HEIGHT',6,'height of plots, in inches');
define('FIGURE_DPI',100,'resolution of plots, in DPI');
  
def plot_uvcov (msname="$MS",width=None,height=None,dpi=None,save=None,select=None,limit=None,**kw):
  """Makes uv-coverage plot
  'msname' is superglobal MS by default.
  If 'save' is given, saves figure to file.
  Use width/height/dpi to override figure settings.
  Any additional keyword arguments are passed to plot(). Try e.g. ms=.1 to change the marker size.
  """
  msname,save = interpolate_locals("msname save");
  tab = ms(msname);
  if select:
    tab = tab.query(select);
  uv = tab.getcol("UVW")[:,:2];
  import matplotlib
  if save:
    matplotlib.use('agg');
#  else:
#    matplotlib.use('qt4agg');
  import pylab
  pylab.figure(figsize=(width or FIGURE_WIDTH,height or FIGURE_HEIGHT));
  pylab.plot(-uv[:,0],-uv[:,1],'.r',**kw);
  pylab.plot(uv[:,0],uv[:,1],'.b',**kw);
  if limit is not None:
    pylab.xlim(-limit,limit);
    pylab.ylim(-limit,limit);
  if save:
    pylab.savefig(save,dpi=(dpi or FIGURE_DPI));
    info("saved UV coverage plot to $save");
  else:
    pylab.show();
  
document_globals(plot_uvcov,"FIGURE_*");  

##
## ARCHIVE/UNARCHIVE FUNCTIONS
##
def_global("TARBALL_DIR",".","directory where tarballs of MSs are kept");

def load_tarball (msname="$MS"):
  """Unpacks fresh copy of MS from .tgz in ms.TARBALL_DIR"""; 
  msname = interpolate_locals("msname");
  if exists(msname):
    info("$msname exists, removing and unpacking fresh MS from tarball");
    x.sh("rm -fr $msname")
  else:
    info("$msname does not exist, unpacking fresh MS from tarball");
  x.sh("cd ${msname:DIR}; tar zxvf ${TARBALL_DIR}/${msname:FILE}.tgz");

def save_tarball (msname="$MS"):
  """Saves MS to .tgz in ms.TARBALL_DIR""";
  msname = interpolate_locals("msname");
  x.sh("cd ${msname:DIR}; tar zcvf ${TARBALL_DIR}/${msname:FILE}.tgz ${msname:FILE}");
                                
document_globals(load_tarball,"TARBALL_DIR");
document_globals(save_tarball,"TARBALL_DIR");


##
## MERGE/SPLIT/VIEW FUNCTIONS
##
def merge (output="merged.MS",options=""):
  """Merges the MSs given by MS_List into the output MS. Options are passed to merge-ms."""
  if not v.MS_List:
    abort("MS_List must be set before calling ms.merge()");
  output,options = interpolate_locals("output options");
  mergems("-f",options,output,*v.MS_List);
document_globals(merge,"MS_List");


def split_views (msname="$MS",output="${MS:DIR}/${MS:BASE}-%s.MS",column="OBSERVATION_ID",values=None):
  """Splits MS into views according to unique values of the column. If values is not specified,
  looks in the column for the set of unique values. Makes a set of reference MSs, using 'output' as
  a template.""";
  msname,output,column = interpolate_locals("msname output column");
  if not values:
    values = set(ms(msname).getcol(column));
  info("splitting $msname by $column (%s)"%" ".join(map(str,values)));
  for val in values:
    subset = output%val;
    taql("SELECT FROM $msname where $column == $val giving $subset",split_args=False);

    
def virtconcat (output="concat.MS",thorough=False,subtables=False):
  """Virtually concatenates the MSs given by MS_List into an output MS."""
  output,options = interpolate_locals("output options");
  auxfile = II("${output}.dat");
  if not v.MS_List:
    abort("MS_List must be set before calling ms.virtconcat()");
  if len(v.MS_List) != len(set([os.path.basename(x) for x in v.MS_List])):
    abort("ms.virtconcat: MSs to be concatenated need to have unique basenames. Please rename your MSs accordingly.");
  cmd = "";
  if thorough:
    cmd = """ms.open("%s", nomodify=False)\n"""%v.MS_List[0];
    for msl in v.MS_List[1:]:
      cmd += II("""ms.virtconcatenate("$msl","$auxfile",'1GHz','1arcsec')\n""");  
    cmd += II("""ms.close()\nos.remove('$auxfile')\n""");
  cmd += """if not ms.createmultims('$output',["%s"],
  [],
  True, # nomodify  
  False,# lock  
  $subtables): # copysubtables from first to all other members  
  os._exit(1);
ms.close();
"""%'","'.join(v.MS_List); 
  std.runcasapy(cmd);
  info("virtually concatenated %d inputs MSs into $output"%len(v.MS_List));

document_globals(virtconcat,"MS_List");


##
## CONVERSION
##
def from_uvfits (fitsfile,msname="$MS"):
  fitsfile,msname = interpolate_locals("fitsfile msname");
  if not msname:
    msname = fitsfile+".MS";
  std.runcasapy("""ms.fromfits(msfile='$msname',fitsfile='$fitsfile')""");

##
## RESAMPLING FUNCTIONS
##
def rebin_freq (msname="$MS",output="$MSOUT",step=1):
  """Resamples MS in frequency with the given step size""";
  msname,output = interpolate_locals("msname output");
  std.runcasapy("""ms.open("$msname"); ms.split(outputms='$output',step=$step);""");
  



##
## FLAGGING FUNCTIONS
##
def aoflagger (msname="$MS",strategy=None):
  """Runs AOFlagger with the specified strategy"""
  msname,strategy = interpolate_locals("msname strategy");
  if strategy:
    _aoflagger("-strategy $strategy $msname");
  else:
    _aoflagger(msname);


def flag_ifrs (msname="$MS",ifrs="",flagset="badifr"):
  """Flags specified baselines""";
  msname,ifrs,flagset = interpolate_locals("msname ifrs flagset");
  _flagms(msname,"-I $ifrs -f $flagset -c");


def_global("FLAG_CHANNELS_MULTIPLIER",1,"multiply channel numbers given to flag_channels() by N");
def_global("FLAG_TIMESLOTS_MULTIPLIER",1,"multiply timeslot numbers given to flag_timeslots() by N");

def flag_channels (msname="$MS",begin=0,end=0,ifrs="all",flagset="badchan"):
  """Flags specified channel range in the specified baselines""";
  msname,ifrs,flagset = interpolate_locals("msname ifrs flagset");
  begin *= FLAG_CHANNELS_MULTIPLIER;
  end *= FLAG_CHANNELS_MULTIPLIER;
  _flagms(msname,"-I $ifrs -L $begin~$end -f $flagset -c");
document_globals(flag_channels,"FLAG_CHANNELS_*");
  

def flag_timeslots (msname="$MS",begin=0,end=0,ifrs="all",flagset="badts"):
  """Flags specified channel range in the specified baselines""";
  msname,ifrs,flagset = interpolate_locals("msname ifrs flagset");
  begin *= FLAG_TIMESLOTS_MULTIPLIER;
  end *= FLAG_TIMESLOTS_MULTIPLIER;
  _flagms(msname,"-I $ifrs -T $begin~$end -f $flagset -c");
document_globals(flag_channels,"FLAG_TIMESLOTS_*");

  
  

###
### Various MS-related settings
###


## current spwid and number of channels. Note that these are set automatically from the MS by the _msddid_Template below
SPWID = 0
TOTAL_CHANNELS = 0

## whenever the MS or DDID changes, look up the corresponding info on channels and spectral windows 
_msddid = None;
def _msddid_Template ():
  global SPWID,TOTAL_CHANNELS,_ms_ddid;
  if II("$MS:$DDID") != _msddid and II("$MS") and DDID is not None:
    try:
      SPWID = ms(MS,"DATA_DESCRIPTION").getcol("SPECTRAL_WINDOW_ID",DDID,1)[0];
      TOTAL_CHANNELS = ms(MS,"SPECTRAL_WINDOW").getcol("NUM_CHAN",SPWID,1)[0];
      # make sure this is reevaluated
      _chanspec_Template();
      info("$MS ddid $DDID is spwid $SPWID with $TOTAL_CHANNELS channels"); 
    except:
      warn("Error accessing $MS");
#      traceback.print_exc();
      return None;
  return II("$MS:$DDID");

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

