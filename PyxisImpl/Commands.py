import sys
import traceback
import inspect
import os
import os.path
import time
import itertools

import PyxisImpl
import PyxisImpl.Internals
from PyxisImpl.Internals import _int_or_str,interpolate,run,loadconf,assign_templates,register_pyxis_module

def _init (context):
  global x;
  global xo;
  global xz;
  global v;
  global E;
  # add some standard objects
  # the 'x' object is a shortcut for executing shell commands. E.g. x.ls('.')
  # the 'xo' object is a shortcut for executing shell commands that are allowed to fail. E.g. x.ls('.')
  x = PyxisImpl.Internals.ShellExecutorFactory(allow_fail=False);
  x.__name__ = 'x';

  xo = PyxisImpl.Internals.ShellExecutorFactory(allow_fail=True);
  xo.__name__ = 'xo';
  
  xz = PyxisImpl.Internals.ShellExecutorFactory(allow_fail=True,bg=True);
  xz.__name__ = 'xz';

  v = PyxisImpl.Internals.OptionalVariable(context,_assign_global);
  object.__setattr__(v,'__name__','v');
  
  E = PyxisImpl.Internals.ShellVariable();
  object.__setattr__(E,'__name__','E');


def _I (string,level=1):
  """_I(string): interpolates (i.e. replaces by value) $VAR and ${VAR} occurrences of local and 
  global variables. Local variables are interpreted as those in the context of the caller (i.e. 1 level away).
  
  _i(string,level): change the number of caller levels for lookup of locals, e.g. 2 means caller or caller
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
  PyxisImpl.Internals.assign_templates();
  frame = inspect.currentframe().f_back;
  ret = [ x and interpolate(x,frame) for x in strings ];
  return ret[0] if len(strings)<2 else ret;

II = _II;  
  
def interpolate_locals (*varnames):
  """interpolates the variable names (from the local context) given by its argument(s).
  Returns new values in the order given. Useful as the opening line of a function, for example:
  
  def f(a="$A",b="$a.1"):
    a,b = interpolate_locals("a","b")   # or "a", "b"
    
  will assign the global variable A to a, and the value of a (in this case, the value of A) plus ".1" to b.
  """;
  PyxisImpl.Internals.assign_templates();
  # interpolate the whole locals() dict
  frame = inspect.currentframe().f_back;
  locs = interpolate(frame.f_locals,frame,depth=2);
  # return variables in the order listed
  ret = [ locs.get(name) for name in itertools.chain(*[ v.split(" ") for v in varnames ]) ];
  return ret if len(ret) != 1 else ret[0];
    
def _timestamp ():
  return time.strftime("%Y/%m/%d %H:%M:%S");
  
def _debug (*msg):
  """Prints debug message(s) without interpolation""";
  print _timestamp(),"DEBUG:",' '.join(map(str,msg));

def _verbose (level,*msg):
  """Prints verbosity message(s), if verbosity level is >= Context.pyxis_verbosity""";
  try:
    verb = int(PyxisImpl.Context['VERBOSE']);
  except:
    verb = 1;
  if level <= verb:
    print "PYXIS:",' '.join(map(str,msg));

def _info (*msg):
  """Prints info message(s) without interpolation""";
  print _timestamp(),"INFO:",' '.join(map(str,msg));

def _warn (*msg):
  """Prints warning message(s) without interpolation""";
  print _timestamp(),"WARNING:",' '.join(map(str,msg));

def _abort (*msg):
  """Prints error message(s) without interpolation and aborts""";
  print _timestamp(),"ERROR:",' '.join(map(str,msg));
  PyxisImpl.Internals.flush_log();
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

def abort (*msg):
  """Prints error message(s) and aborts""";
  _abort(*[ _I(x,2) for x in msg ]);
  
def printvars ():
  """Prints the current variable settings in Context""";
  # make sorted list of globals (excepting those starting with underscore)
  globs = [ (name,value) for name,value in sorted(PyxisImpl.Context.iteritems()) 
            if name[0]!="_" and name not in PyxisImpl._predefined_names ];
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
    if callable(value) and not isinstance(value,PyxisImpl.Internals.ShellExecutor) and not name.endswith("_Template") and not name in globals() ]);
  print "External tools:";
  for name,value in globs:
    if isinstance(value,PyxisImpl.Internals.ShellExecutor):
      print "  %s=%s"%(name,value.path);
  print "Pyxis built-ins:";
  print " ",", ".join([ name for name,value in globs 
    if callable(value) and not isinstance(value,PyxisImpl.Internals.ShellExecutor) and not name.endswith("_Template") and name in globals() ]);

def exists (filename):
  return os.path.exists(_I(filename,2));
  
def _assign_global (name,value,interpolate=True,frame=None,verbose_level=2):
  assign(name,value,namespace=PyxisImpl.Context,interpolate=interpolate,frame=frame,verbose_level=verbose_level);
  
def assign (name,value,namespace=None,interpolate=True,frame=None,verbose_level=2):
  frame = frame or inspect.currentframe().f_back;
  # find namespace
  if not namespace and '.' in name:
    nsname,name = name.rsplit(".",1);
    namespace = PyxisImpl._namespaces.get(nsname);
    if namespace is None:
      raise ValueError,"invalid namespace %s"%nsname;
  elif not namespace:
    namespace = frame.f_globals;
  # templates are not interpolated, all others are
  if name.endswith("_Template"):
    _verbose(verbose_level,"setting %s=%s"%(name,value));
    namespace[name] = value;
  else:
    value1 = PyxisImpl.Internals.interpolate(value,frame) if interpolate else value;
    namespace[name] = value1;
    _verbose(verbose_level,"setting %s=%s%s"%(name,value1,(" (%s)"%value) if str(value) != str(value1) else ""));
  PyxisImpl.Internals.assign_templates();

def _per (varname,*commands):
  saveval = PyxisImpl.Context.get(varname,None);
  varlist = PyxisImpl.Context.get(varname+"_List",None);
  if varlist is None:
    _verbose(1,"per(%s,%s): %s_List is empty"%(varname,",".join(map(str,commands)),varname));
    return;
  if type(varlist) is str:
    varlist = map(_int_or_str,varlist.split(","));
  elif not isinstance(varlist,(list,tuple)):
    _abort("PYXIS: per(%s,%s): %s_List has invalid type %s"%(varname,",".join(map(str,commands)),varname,str(type(varlist))));
  _verbose(1,"per(%s,%s): iterating over %s=%s"%(varname,",".join(map(str,commands)),varname," ".join(map(str,varlist))));
  # do the actual iteration
  for value in varlist:
    assign(varname,_I(value,3),interpolate=False);
    PyxisImpl.Internals.run(*commands);
  # assign old value
  if saveval is not None:
    globals()[varname] = saveval;
    _verbose(1,"restoring %s=%s"%(varname,saveval));
    PyxisImpl.Internals.assign_templates();

def per (varname,*commands):
  _per(varname,*commands);

def per_ms (*commands):
  _per("MS",*commands);

