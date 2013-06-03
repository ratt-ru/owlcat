"""Pyxis.Commands: implements commands, etc. for interactive use, as well as for Pyxides modules"""

import sys
import traceback
import inspect
import os
import os.path
import time
import itertools
import math

import Pyxis
import Pyxis.Internals
from Pyxis.Internals import _int_or_str,interpolate,loadconf,assign,assign_templates

import numpy as np

DEG = math.pi/180
ARCMIN = DEG/60
ARCSEC = ARCMIN/60

def _init (context):
  global x;
  global xo;
  global xz;
  global v;
  global E;
  # add some standard objects
  # the 'x' object is a shortcut for executing shell commands. E.g. x.ls('.')
  # the 'xo' object is a shortcut for executing shell commands that are allowed to fail. E.g. x.ls('.')
  x = Pyxis.Internals.ShellExecutorFactory(allow_fail=False);
  x.__name__ = 'x';

  xo = Pyxis.Internals.ShellExecutorFactory(allow_fail=True);
  xo.__name__ = 'xo';
  
  xz = Pyxis.Internals.ShellExecutorFactory(allow_fail=True,bg=True);
  xz.__name__ = 'xz';

  v = Pyxis.Internals.GlobalVariableSpace(context);
  object.__setattr__(v,'__name__','v');
  
  E = Pyxis.Internals.ShellVariableSpace();
  object.__setattr__(E,'__name__','E');


def _I (string,level=1):
  """_I(string): interpolates (i.e. replaces by value) $VAR and ${VAR} occurrences of local and 
  global variables. Local variables are interpreted as those in the context of the caller (i.e. 1 level away).
  
  _i(string,level): change the number of caller levels for lookup of locals, e.g. 2 means caller of caller
  _i(string,-1): interpolate across all callers
  """;
  frame = inspect.currentframe();
  for i in range(level):
    frame = frame and frame.f_back;
  return interpolate(string,frame);
  
def _II (*strings):
  """_II(string): interpolates multiple strings. Does an implicit assign_templates call beforehand. Useful as
  the opening line of a function, as e.g. arg1,arg2 = _II(arg1,arg2);
  """;
  Pyxis.Internals.assign_templates();
  frame = inspect.currentframe().f_back;
  ret = [ x and interpolate(x,frame) for x in strings ];
  return ret[0] if len(strings)<2 else ret;

II = _II;  
  
_subprocess_id = None;  
  
def _timestamp ():
  ts = time.strftime("%Y/%m/%d %H:%M:%S");
  if _subprocess_id is not None:
    ts += " [%d]"%_subprocess_id;
  return ts;

def _message (*msg,**kw):
  output = " ".join(map(str,msg));
  if sys.stdout is not sys.__stdout__:
    print output;
    if ( kw.get('critical') or not Pyxis.Context.get('QUIET') ):
      sys.__stdout__.write(output+"\n");
  else:
    if not Pyxis.Context.get('QUIET'):
      print output;
  
def _debug (*msg):
  """Prints debug message(s) without interpolation""";
  _message(_timestamp(),"DEBUG:",*msg);

def _verbose (level,*msg):
  """Prints verbosity message(s), if verbosity level is >= Context.pyxis_verbosity""";
  try:
    verb = int(Pyxis.Context['VERBOSE']);
  except:
    verb = 1;
  if level <= verb:
    _message("PYXIS%s:"%((" [%d]"%_subprocess_id) if _subprocess_id is not None else ""),*msg);

def _info (*msg):
  """Prints info message(s) without interpolation""";
  _message(_timestamp(),"INFO:",*msg);

def _warn (*msg):
  """Prints warning message(s) without interpolation""";
  _message(_timestamp(),"WARNING:",*msg);

def _error (*msg):
  """Prints error message(s) without interpolation""";
  _message(_timestamp(),"ERROR:",*msg);

def _abort (*msg):
  """Prints error message(s) without interpolation and aborts""";
  _message(_timestamp(),"ABORT:",critical=True,*msg);
  Pyxis.Internals.flush_log();
  sys.exit(1);
  
def debug (*msg):
  """Prints debug message(s)""";
  _debug(*[ _I(x,2) for x in msg ]);

def verbose (level,*msg):
  """Prints verbosity message(s), if verbosity level is >= Context.pyxis_verbosity""";
  _verbose(level,*[ _I(x,2) for x in msg ]);

def info (*msg):
  """Prints info message(s)""";
  _info(*[ _I(x,2) for x in msg ]);

def warn (*msg):
  """Prints warning message(s)""";
  _warn(*[ _I(x,2) for x in msg ]);

def error (*msg):
  """Prints error message(s)""";
  _error(*[ _I(x,2) for x in msg ]);

def abort (*msg):
  """Prints error message(s) and aborts""";
  _abort(*[ _I(x,2) for x in msg ]);
  
def printvars ():
  """Prints the current variable settings in Context""";
  # make sorted list of globals (excepting those starting with underscore)
  globs = [ (name,value) for name,value in sorted(Pyxis.Context.iteritems()) 
            if name[0]!="_" and name not in Pyxis._predefined_names ];
  # from this, extract the non-callables
  varlist = [ (name,value) for name,value in globs if not callable(value)  ];
  print "Ordinary variables:";
  for var,val in varlist:
    if not var.endswith("_List") and not var.endswith("_Template") and isinstance(val,(str,int)):
      print "  %s=%s"%(var,val);
  print "Lists:";
  for var,val in varlist:
    if var.endswith("_List"):
      print "  %s=%s"%(var,val);
  print "Templates:";
  for var,val in globs:
    if var.endswith("_Template"):
      print "  %s=%s"%(var,val);
  # now deal with the callables
  print "Functions:";
  print " ",", ".join([ name for name,value in globs 
    if callable(value) and not isinstance(value,Pyxis.Internals.ShellExecutor) and not name.endswith("_Template") and not name in globals() ]);
  print "External tools:";
  for name,value in globs:
    if isinstance(value,Pyxis.Internals.ShellExecutor):
      print "  %s=%s"%(name,value.path);
  print "Pyxis built-ins:";
  print " ",", ".join([ name for name,value in globs 
    if callable(value) and not isinstance(value,Pyxis.Internals.ShellExecutor) and not name.endswith("_Template") and name in globals() ]);

def exists (filename):
  """Returns True if filename exists, interpolating the filename""";
  return os.path.exists(_I(filename,2));

def _per (varname,*commands):
  saveval = Pyxis.Context.get(varname,None);
  varlist = Pyxis.Context.get(varname+"_List",None);
  cmdlist = ",".join([ x if isinstance(x,str) else getattr(x,"__name__","?") for x in commands ]);
  if varlist is None:
    _verbose(1,"per(%s,%s): %s_List is empty"%(varname,cmdlist,varname));
    return;
  try:
    if type(varlist) is str:
      varlist = map(_int_or_str,varlist.split(","));
    elif not isinstance(varlist,(list,tuple)):
      _abort("PYXIS: per(%s,%s): %s_List has invalid type %s"%(varname,cmdlist,str(type(varlist))));
    nforks = Pyxis.Context.get("JOBS",0);
    stagger = Pyxis.Context.get("JOB_STAGGER",0);
    # unforked case
    _verbose(1,"per(%s,%s): iterating over %s=%s"%(varname,cmdlist,varname," ".join(map(str,varlist))));
    global _subprocess_id;
    if nforks < 2 or len(varlist) < 2 or _subprocess_id is not None:
      # do the actual iteration
      for value in varlist:
        assign(varname,value,namespace=Pyxis.Context,interpolate=False);
        Pyxis.Internals.run(*commands);
    else:
      # else split varlist into forked subprocesses
      per_fork = max(len(varlist)//nforks,1);
      _verbose(1,"splitting into %d jobs, %d %s's per job, staggered by %ds"%(min(nforks,len(varlist)),per_fork,varname,stagger));
      forked_pids = set();
      for i in range(0,len(varlist),per_fork):
        if i and stagger:
          time.sleep(stagger);
        pid = os.fork();
        if not pid:
          # child fork: run commands
          _subprocess_id = i//per_fork;
          try:
            for value in varlist[i:i+per_fork]:
              assign(varname,value,namespace=Pyxis.Context,interpolate=False);
              Pyxis.Internals.run(*commands);
          except:
            traceback.print_exc();
            sys.exit(1);
          sys.exit(0);
        else: # parent pid: append to list
          forked_pids.add(pid);
      _verbose(1,"%d forked subprocesses launched, waiting for finish"%len(forked_pids));
      while forked_pids:
        pid,status = os.waitpid(-1,0);
        if pid in forked_pids:
          forked_pids.remove(pid);
          status <<= 8;
          if status:
            _abort("subprocess %d exited with error status %d, aborting"%(pid,status));
          else:
            _verbose(1,"subprocess %d finished, waiting for %d more"%(pid,len(forked_pids)));
  finally:
    # assign old value
    if saveval is not None:
      assign(varname,saveval,namespace=Pyxis.Context,interpolate=False);
      _verbose(1,"restoring %s=%s"%(varname,saveval));
      Pyxis.Internals.assign_templates();

def per (varname,*commands):
  """Iterates over variable 'varname', and executes commands. That is, for every value
  in varname_List, sets varname to that value, then calls the commands."""
  _per(varname,*commands);

def per_ms (*commands):
  """Iterates over variable 'MS', and executes commands. That is, for every value
  in MS_List, sets MS to that value, then calls the commands."""
  _per("MS",*commands);

def per_ddid (*commands):
  """Iterates over variable 'DDID', and executes commands. That is, for every value
  in DDID_List, sets MS to that value, then calls the commands."""
  _per("DDID",*commands);

def is_true (arg):
  """Returns True if argument evaluates to boolean truth, possibly as a string.
  is_true('False') == is_true('0') == is_true(0) == False
  is_true('True') == is_true('1') == is_true(1) == True
  """;
  if isinstance(arg,str):
    try:
      return bool(eval(arg));
    except:
      raise TypeError,"is_true('%s'): invalid string argument"%arg;
  try:
    return bool(arg);
  except:
    raise TypeError,"is_true(%s): invalid argument %s"%(str(arg),str(type(arg)));

import tempfile   
import cPickle
import fcntl
    
class Safelist (object):
  """A Safelist provides a list of objects that is shared among parallel processes. Multiple jobs launched
  by pyxis can append to a Safelist in a multiprocess-safe manner: write locks are enforced""";
  def __init__ (self,filename=None):
    """Creates a safelist, associated with the given filename. If not supplied, a filename is
    chosen randomly"""
    if filename:
      filename = _I(filename,2);
      if os.path.exists(filename):
        os.remove(filename);
      self.filename = filename;
    else:
      self.filename = tempfile.NamedTemporaryFile(delete=False).name;

  def reset (self):
    """Resets safelist to empty (by deleting the associated file)""";
    if os.path.exists(self.filename):
      os.remove(self.filename);
      
  def add (self,obj):
    """Adds an object to the safelist, in an MP-safe manner""";
    if isinstance(obj,str):
      obj = _I(obj,2);
    ff = file(self.filename,"ab");
    fcntl.flock(ff,fcntl.LOCK_EX);
    try:
      cPickle.dump(obj,ff);
    finally:
      fcntl.flock(ff,fcntl.LOCK_UN);
    
  def read (self):
    """Reads all objects accumulated in the safelist, in an MP-safe manner""";
    ret = [];
    if os.path.exists(self.filename):
      ff = file(self.filename);
      fcntl.flock(ff,fcntl.LOCK_EX);
      try:
        while True:
          ret.append(cPickle.load(ff));
      except EOFError:
        pass;
      finally:
        fcntl.flock(ff,fcntl.LOCK_UN);
    return ret;
    