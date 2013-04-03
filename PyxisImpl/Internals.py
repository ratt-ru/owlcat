import glob
import traceback
import subprocess
import re
import string
import os
import os.path
import inspect
import sys
import itertools

import PyxisImpl

def init (context):
  """init internals, attach to the given context""";
  global _debug
  global _info
  global _abort
  global _verbose
  global _warn
  from PyxisImpl.Commands import _debug,_info,_abort,_verbose,_warn
  # set default verbosity to 1
  preset_verbosity = context.get("VERBOSE",None);
  context.setdefault("VERBOSE",1);
  # set context and init stuff
  PyxisImpl.Context = context;
  PyxisImpl._predefined_names = set(context.iterkeys());
  PyxisImpl.Commands._init(context);
  # report verbosity
  if preset_verbosity is None and context['VERBOSE'] == 1:
    _verbose(1,"VERBOSE=1 by default");
  else:
    _verbose(1,"VERBOSE=%d"%context['VERBOSE']);

def _int_or_str (x):
  """helper function: converts argument to int if possible, else returns string""";
  try:
    return int(x);
  except:
    return str(x);
    
def _interpolate_args (args,kws,localdict): 
  """Helper function to interpolate argument list and keywords using the local dictionary, plus PyxisImpl.Context globals""";
  return [ interpolate(arg,localdict) for arg in args ], \
         dict([ (kw,interpolate(arg,localdict)) for kw,arg in kws.iteritems() ]);
  
class ShellExecutor (object):    
  """This is a ShellExecutor object, which can be used to execute a particular shell command. 
The command (and any fixed arguments) are specified in the constructor of the object. ShellExecutors are
typically created via the Pyxis x, xo or xz built-ins, e.g. as 
      ls = x.ls 
      lsl = x.ls("-l")
      ls()             # runs "ls"
      dir="./test"
      lsl("$dir")      # runs "ls -l ./test"
      ls(x=1)          # runs "ls x=1"
  """;
  
  def __init__ (self,name,path,allow_fail=False,bg=False,localdict={},*args,**kws):
    self.name,self.path = name,path;
    self.allow_fail = allow_fail;
    self.bg  = bg;
    self._add_args,self._add_kws = _interpolate_args(args,kws,localdict);
    
  def args (self,*args,**kws):
    """Creates instance of executor with additional args. Local variables of caller are interpolated."""
    localdict = inspect.currentframe().f_back.f_locals;
    args,kws = _interpolate_args(args,kws,localdict);
    kws0 = dict(self._add_kws);
    kws0.update(kws);
    return ShellExecutor(self.name,self.path,self.allow_fail,{},*(self._add_args+args),**kws0);
    
  def __str__ (self):
    return " ".join([self.path]+self._add_args+["%s=%s"%(a,b) for a,b in self._add_kws.iteritems()]);
    
  def __repr__ (self):
    return "ShellExecutor: %s"%str(self);
    
  def __call__ (self,*args,**kws):
    """Runs the associated shell command, with additional supplied arguments. Normal arguments are simply converted
    to strings. Keywords are converted to key=value arguments. Local variables of the caller are interpolated."""
    if self.path is None:
      if self.allow_fail:
        _abort("PYXIS: shell command '%s' not found"%self.name);
      _warn("PYXIS: shell command '%s' not found"%self.name);
    else:
      localdict = inspect.currentframe().f_back.f_locals;
      args,kws = _interpolate_args(args,kws,localdict);
      kws0 = dict(self._add_kws);
      kws0.update(kws);
      return _call_exec(self.path,allow_fail=self.allow_fail,bg=self.bg,*(self._add_args+args),**kws0);

class ShellExecutorFactory (object):
  """The Pyxis "x", "xo" and "xz" built-ins can be used to create proxies for shell commands called
ShellExecutors. For example:
    ls = x.ls     # the ls object is now a ShellExecutor for the ls command
    ls()          # executes ls
    ls("-l")      # executes ls -l
    ls = x("run-imager.sh")  # alternative syntax for creating a ShellExecutor, in this case for run-imager.sh
Executors created with 'x' are mandatory, while those created with 'xo' are optional. A mandatory executor
will terminate the Pyxis script if its command fails. An optional executor will report an error and continue.
Executors created via 'xz' are optional, and run commands in the background.
"""  
  def __init__ (self,allow_fail=False,bg=False):
    self.allow_fail = allow_fail;
    self.bg = bg;
    
  def __getattr__ (self,command,default=None):
    """Creates a ShellExecutor for a given shell command. For example,
    ls = x.ls     # the ls object is now a ShellExecutor for the ls command
    ls()          # executes ls
    ls("-l")      # executes ls -l
 """
    if command.find('/') >= 0:
      path = command if os.access(command,os.X_OK) else None;
    else:
      path = find_exec(command);
    return ShellExecutor(command,path,self.allow_fail,self.bg,inspect.currentframe().f_back.f_locals);
    
  def __call__ (self,*args):
    """An alternative way to make ShellExecutors, e.g. as x("command arg1 arg2").
    Useful when the command contains e.g. dots or slashes, thus making the x.command syntax unsuitable."""
    if len(args) == 1:
      args = args[0].split(" ");
    return ShellExecutor(args[0],args[0],self.allow_fail,self.bg,inspect.currentframe().f_back.f_locals,*args[1:]);
    
  def sh (self,*command):
    """Directly invokes the shell with a command and arguments"""
    commands = [ str(interpolate(command,inspect.currentframe().f_back.f_locals)) for command in commands ];
    # run command
    _verbose(1,"executing '%s':"%(" ".join(commands)));
    flush_log();
    retcode = subprocess.call(*commands,shell=True,stdout=sys.stdout,stderr=sys.stderr);
    if retcode:
      if self.allow_fail:
        _warn("PYXIS: '%s' returns error code %d"%(commands[0],retcode));
        return;
      else:
        _abort("PYXIS: '%s' returns error code %d"%(commands[0],retcode));
    else:
      _verbose(2,"'%s' succeeded"%commands[0]);
      
  def __repr__ (self):
    name = self.__name__;
    return "Pyxis built-in %s: access to shell commands. Use help(%s) for details."%(name,name);
    

class OptionalVariable (object):
  """This object provides "smart" access to global Pyxis variables. Note that most Pyxis code runs with
globals directly accessible as Python variables anyway, but "v" provides a number of extra features:

  v.VARNAME evaluates to the variable VARNAME, or to the empty string, if VARNAME is not defined
  v('VARNAME',default) evaluates to the variable VARNAME, or to default if VARNAME is not defined.
     If default is a string, local variables of the caller are interpolated.
  v.VARNAME=value assigns to a global variable, and causes templates to be re-evaluated, and other 
    implicit variable-related actions to be taken. In particular, v.LOG="logfile" will set a new
    log destination.
  """;
  
  def __init__ (self,namespace,assigner=None):
    object.__setattr__(self,'namespace',namespace);
    object.__setattr__(self,'_assign_func',assigner);
    
  def __call__ (self,name,default=""):
    if isinstance(default,str):
      default = interpolate(default,inspect.currentframe().f_back.f_locals);
    return object.__getattribute__(self,'namespace').get(name,default);
    
  def __getattr__ (self,attr,default=""):
    return self(attr,default);
    
  def __setattr__ (self,attr,value):
    assigner = object.__getattribute__(self,'_assign_func') or object.__setattr__;
    return assigner(attr,value); 
    
  def __repr__ (self):
    name = object.__getattribute__(self,'__name__');
    return "Pyxis built-in %s: smart access to global variables. Try %s.VARNAME, or help(%s)."%(name,name,name);
    
class ShellVariable (OptionalVariable):
  """This object provides quick access to environment (i.e. shell) variables. Use e.g.
  E.HOME or E("HOME",default) to access a shell variable. If default is a string, local 
  variables of the caller are interpolated.
  """;
  def __init__ (self):
    OptionalVariable.__init__(self,os.environ);
    
  def __repr__ (self):
    name = object.__getattribute__(self,'__name__');
    return "Pyxis built-in %s: quick access to shell variables. Try %s.VARNAME, or help(%s)."%(name,name,name);

def interpolate (arg,localdicts=[],ignore=set(),skip=set()):
  """Interpolates strings: substitutes $var and ${var} with the corresponding variable value from PyxisImpl.Context.
  Also does the old-style %{var}s interpolation.
  
  If arg is a string, does interpolation and returns new string.
  
  If arg is a dict, does interpolation on every string-type key in the dict (except for those in 'skip'), using the 
    dict itself as a source of symbols (plus the global variables). Returns the dict.
  
  If set, 'localdicts' is local dict where additional symbols will be looked up, or a list of dicts 
  (with symbols in the head of the list overriding the tail).
  
  If set, 'ignore' is a container of symbols which will interpolate to an empty string.
  """;
  if isinstance(arg,dict):
    arg = arg.copy();
    # interpolate until things stop changing, but quit after 20 loops
    for count in range(20):
      updates = {};
      for key,value in arg.iteritems():
        # interpolate string variables
        if key not in skip and isinstance(value,str):
          newvalue = interpolate(value,[arg],ignore=[value]);
#          print "%s: %s->%s"%(key,value,newvalue);
          if newvalue != value:
            updates[key] = newvalue;
      # apply updates, unless things stop
#      print updates;
      if not updates:
        break;
      arg.update(updates);
    return arg;
  # strings are interpolated
  elif isinstance(arg,str):
    syms = DictProxy(PyxisImpl.Context);
    for ld in ([localdicts] if isinstance(localdicts,dict) else localdicts[-1::-1]):
      syms.update(ld);
    syms.remove(ignore);
    return string.Template(str(arg)%syms).safe_substitute(syms);
  # all other types returned as-is
  else:
    return arg;

class DictProxy (object):
  def __init__ (self,dict0):
    self.dict0 = dict0;
    
  def update (self,dict1):
    self.dict0.update(dict1);
    
  def remove (self,what):
    for key in what:
      self.dict0.pop(key,None);
    
  def __getitem__ (self,item):
    if item in self.dict0:
      return self.dict0.get(item);
    elif item.endswith("_BASE") and item[:-5] in self.dict0:
      base = self.dict0.get(item[:-5]);
      if not isinstance(base,str):
        return "";
      while base and base[-1] == "/":
        base = base[:-1];
      return os.path.splitext(base)[0];
    else:
      return "";
    
  def __contains__ (self,item):
    return True;
    
def assign_templates ():
  """For every variable in PyxisImpl.Context that ends with "_Template", assigns value to it by interpolating the template.""";
  updated = True;
  count = 1000;
  while updated:
    updated = False;
    count -= 1;
    if not count:
      _abort("Too many template assignment steps. This can be caused by variable templates that cross-reference each other");
    newvalues = {};
    # interpolate new values for each variable that has a _Template equivalent
    for var,value in PyxisImpl.Context.iteritems():
      if var.endswith("_Template"):
        # string templates are interpolated, callable ones are called
        if isinstance(value,str):
          newvalues[var[:-len("_Template")]] = interpolate(value);
        elif callable(value):
          newvalues[var[:-len("_Template")]] = value();
    # update dict
    for var,value in newvalues.iteritems():
      oldval = PyxisImpl.Context.get(var,None);
      if oldval != value:
        updated = True;
        PyxisImpl.Context[var] = value;
        _verbose(2,"%s templated value %s=%s"%("initialized" if oldval is None else "updated",var,value));
  set_logfile(PyxisImpl.Context.get('LOG',None));

_current_logfile = None;
_current_logobj = None;

def flush_log ():
  _current_logobj and _current_logobj.flush();

def set_logfile (filename):
  """Starts logging to the specified file""";
  global _current_logfile;
  global _current_logobj;
  if filename and not isinstance(filename,str):
    _warn("invalid LOG variable of type %s, ignoring"%str(type(filename)));
    return;
  if filename == "-" or not filename:
    filename = None;
  if filename != _current_logfile:
    _info("redirecting log output to %s"%(filename or "console"));
    if filename is None:
      sys.stdout,sys.stderr = sys.__stdout__,sys.__stderr__;
      _current_logobj = None;
    else:
      mode = "w";
      if filename[0] == '+':
        filename = filename[1:];
        mode = "a";
      _current_logobj = sys.stdout = sys.stderr = open(filename,"w");
    if _current_logfile:
      _info("log continued from %s"%_current_logfile);
    else:
      _info("log started");
    _current_logfile = filename;
        
def initconf (*files):
  """Loads configuration from specified files, or from default file""";
  if not files:
    files = glob.glob("pyxis*.conf")+glob.glob("pyxis*.py");
  # load config files -- all variable assignments go into the PyxisImpl.Context scope
  if files:
    _verbose(1,"auto-loading config files and scripts from 'pyxis*.{conf,py}'. To disable this, set pyxis_skip_config=True before importing Pyxis");
  for filename in files:
    PyxisImpl.Commands.loadconf(filename);

def loadconf (filename):
  """Loads config file""";
  filename = interpolate(filename,inspect.currentframe().f_back.f_locals);
  verbose(1,"loading %s"%filename);
  load_package("config",filename);

def load_package (name,filename,report=True):
  oldstuff = PyxisImpl.Context.copy();
  try:
    exec(file(filename),PyxisImpl.Context);
  except:
    traceback.print_exc();
    _abort("PYXIS: error parsing %s, see output above for details"%filename);
  newnames =  [ (name,obj) for name,obj in PyxisImpl.Context.iteritems() 
                  if not name.startswith("_") and not name in oldstuff ];
  newvars = sorted([name for name,obj in newnames if not callable(obj) and not inspect.ismodule(obj) 
                                                              and not name.endswith("_Template") ]);
  newfuncs = sorted([name for name,obj in newnames if callable(obj) and not name.endswith("_Template") ]);
  newtemps = sorted([name[:-9] for name,obj in newnames if name.endswith("_Template") ]);
  if newfuncs:
    _verbose(1,"  provides functions:"," ".join(newfuncs));
  if newvars:
    _verbose(1,"  provides variables:"," ".join(newvars));
  if newtemps:
    _verbose(1,"  provides templates for:"," ".join(newtemps));
    
def find_exec (cmd):
  """Finds shell executable in PATH"""
  for path in os.environ["PATH"].split(":"):
    filename = os.path.join(path,cmd);
    if os.access(filename,os.X_OK):
      return filename;
  return None;
  
_bg_processes = [];  
  
def _call_exec (path,*args,**kws):
  """Helper function: calls external program with the given arguments and keywords
  (each kw dict element is turned into a name=value argument)""";
  allow_fail = kws.pop('allow_fail',False);
  bg = kws.pop('bg',False);
  # default is to split each argument at whitespace, but split_args=False passes them as-is
  split = kws.pop('split_args',True);
  args1 = [path];
  if split:
    for arg in args:
      args1 += arg.split(" ") if isinstance(arg,str) else [ str(arg) ];
  else:
    args1 += map(str,args);
  # run command
  args = args1+["%s=%s"%(a,b) for a,b in kws.iteritems()];
  flush_log();
  if bg:
    global _bg_processes;
    po = subprocess.Popen(args);
    _bg_processes.append(po);
    _verbose(1,"executing '%s' in background: pid %d"%(" ".join(args),po.pid));
  else:
    _verbose(1,"executing '%s':"%(" ".join(args)));
    if type(sys.stdout) is not file or type(sys.stderr) is not file:
      retcode = subprocess.call(args);
    else:
      retcode = subprocess.call(args,stdout=sys.stdout,stderr=sys.stderr);
    if retcode:
      if allow_fail:
        _warn("PYXIS: '%s' returns error code %d"%(path,retcode));
        return;
      else:
        _abort("PYXIS: '%s' returns error code %d"%(path,retcode));
    else:
      _verbose(2,"'%s' succeeded"%path);

def find_command (comname):
  """Locates command by name. If command is present (as a callable) in PyxisImpl.Context, returns that.
  Otherwise checks the path for a binary by that name, and returns a callable to call that.
  Else aborts.""";
  comname = interpolate(comname);
  # look for a predefined command
  comcall = PyxisImpl.Context.get(comname);
  if callable(comcall):
    return comcall;
  # else look for a shell command
  if comname[0] == "?":
    comname = comname[1:];
    allow_fail = True;
  else:
    allow_fail = False;
  path = find_exec(comname);
  if path is None:
    _abort("PYXIS: undefined command '%s'"%comname);
  # make callable for this shell command
  return lambda *args:_call_exec(path,allow_fail=allow_fail,*args);

def run (*commands):
  """Runs list of commands""";
  # _debug("running",commands);
  for step,command in enumerate(commands):
    # set step counter
    # interpolate the command
    command = command.strip();
    _info("PYXIS: executing command %s"%command);
    # syntax 1: VAR=VALUE or VAR:=VALUE
    match = re.match("^(\w+)(=)(.*)$",command);
    if match:
      name,op,value = match.groups();
      if name == "LOG":
        set_logfile(value);
      # assign variable -- note that templates are not interpolated
      PyxisImpl.Context.assign(name,value,uselocals=False,verbose=False);
      continue;
    # syntax 2: command(args) or command[args]. command can have a "?" prefix
    match = re.match("^(\\??\w+)\\[(.*)\\]$",command) or re.match("^(\\??\w+)\\((.*)\\)$",command);
    if match:
      comname,comargs = match.groups();
      comcall = find_command(comname);
      # split up arguments
      args = [];
      kws = {};
      for arg in re.split(",," if comargs.find(",,") >=0 else ",",comargs):
        arg = interpolate(arg).strip();
        match = re.match("^(\w+)=(.*)$",arg);
        if match:
          kws[match.group(1)] = _int_or_str(match.group(2));
        else:
          args.append(arg);
      # exec command
      comcall(*args,**kws);
      assign_templates();
      continue;
    # syntax 3: standalone command. This better be found!
    comcall = find_command(command);
    comcall();
    assign_templates();
  
