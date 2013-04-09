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
  # loaded modules
  global _namespaces,_superglobals,_modules;
  _namespaces = dict();
  _superglobals = dict();
  _modules = set();
  # The "v" and "" namespaces correspond to the global context
  # All variables of that context are superglobals.
  _namespaces[''] = _namespaces['v'] = context;
  _superglobals[id(context)] = context;
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
    
def _interpolate_args (args,kws,frame): 
  """Helper function to interpolate argument list and keywords using the local dictionary, plus PyxisImpl.Context globals""";
  return [ interpolate(arg,frame) for arg in args ], \
         dict([ (kw,interpolate(arg,frame)) for kw,arg in kws.iteritems() ]);
  
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
  
  def __init__ (self,name,path,frame,allow_fail=False,bg=False,*args,**kws):
    self.name,self.path = name,path;
    self.allow_fail = allow_fail;
    self.bg  = bg;
    self.argframe = frame;
    self._add_args,self._add_kws = list(args),kws;
    
  def args (self,*args,**kws):
    """Creates instance of executor with additional args. Local variables of caller are interpolated."""
    args0,kws0 = _interpolate_args(self._add_args,self._add_kws,self.argframe);
    kws0.update(kws);
    return ShellExecutor(self.name,self.path,inspect.currentframe().f_back,self.allow_fail,self.bg,*(args0+list(args)),**kws0);
    
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
      args0,kws0 = _interpolate_args(self._add_args,self._add_kws,self.argframe);
      args,kws = _interpolate_args(args,kws,inspect.currentframe().f_back);
      kws0.update(kws);
      return _call_exec(self.path,allow_fail=self.allow_fail,bg=self.bg,*(args0+args),**kws0);

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
    command = interpolate(command,inspect.currentframe().f_back);
    if command.find('/') >= 0:
      path = command if os.access(command,os.X_OK) else None;
    else:
      path = find_exec(command);
    return ShellExecutor(command,path,None,self.allow_fail,self.bg);
    
  def __call__ (self,*args,**kws):
    """An alternative way to make ShellExecutors, e.g. as x("command arg1 arg2").
    Useful when the command contains e.g. dots or slashes, thus making the x.command syntax unsuitable."""
    args,kws = _interpolate_args(args,kws,inspect.currentframe().f_back);
    if len(args) == 1:
      args = args[0].split(" ");
    return ShellExecutor(args[0],args[0],None,allow_fail=self.allow_fail,bg=self.bg,*args[1:],**kws);
    
  def sh (self,*args,**kws):
    """Directly invokes the shell with a command and arguments"""
    commands,kws = _interpolate_args(args,kws,inspect.currentframe().f_back);
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

class _DictAccessor (object):
  """Helper class that maps dicts to attributes""";
  def __init__ (self,namespace):
    object.__setattr__(self,'namespace',namespace);
    
  def __call__ (self,name,default=""):
    if isinstance(default,str):
      default = interpolate(default,inspect.currentframe().f_back);
    return object.__getattribute__(self,'namespace').get(name,default);
    
  def __getattr__ (self,name,default=""):
    if isinstance(default,str):
      default = interpolate(default,inspect.currentframe().f_back);
    return object.__getattribute__(self,'namespace').get(name,default);
    
  def __setattr__ (self,name,value):
    object.__getattribute__(self,'namespace')[name] = value;

class GlobalVariableSpace (_DictAccessor):
  """The "v" object provides access to global Pyxis variables. Global variables assigned with
  v.VARNAME=value are propagated across all imported Pyxis modules. They can then be referenced
  directly as VARNAME. The "v" object itself provides a number of extra features:

  v.VARNAME evaluates to the variable VARNAME, or to the empty string, if VARNAME is not defined
  v('VARNAME',default) evaluates to the variable VARNAME, or to default if VARNAME is not defined.
     If default is a string, local variables of the caller are interpolated.
  v.VARNAME=value assigns to a global variable, and also causes templates to be re-evaluated, and 
    other implicit variable-related actions to be taken. In particular, v.LOG="logfile" will set a 
    new log destination.
  """;
  
  def __init__ (self,namespace):
    _DictAccessor.__init__(self,namespace);
    
  def __setattr__ (self,attr,value):
    assign(attr,value,namespace=self.namespace,frame=inspect.currentframe().f_back); 
    
  def __repr__ (self):
    name = object.__getattribute__(self,'__name__');
    return "Pyxis built-in %s: global variable space. Try %s.VARNAME, or help(%s)."%(name,name,name);
    
class ShellVariableSpace (_DictAccessor):
  """This object provides quick access to environment (i.e. shell) variables. Use e.g.
  E.HOME or E("HOME",default) to access a shell variable. If default is a string, local 
  variables of the caller are interpolated.
  """;
  def __init__ (self):
    _DictAccessor.__init__(self,os.environ);
    
  def __repr__ (self):
    name = object.__getattribute__(self,'__name__');
    return "Pyxis built-in %s: quick access to shell variables. Try %s.VARNAME, or help(%s)."%(name,name,name);

def interpolate (arg,frame,depth=1,ignore=set(),skip=set()):
  """Interpolates strings: substitutes $var and ${var} with the corresponding variable value from 
  (in order of lookup):
  
  * the locals and globals of the given frame (must be a frame object: see inspect module for details)
  * the locals and globals of outer frames to a depth of 'depth' (if >1)
  * the global Pyxis context. 
  
  Alternatively, if frame is a dict, then lookup happens in frame, then the global Pyxis context.
  
  If arg is a string, does interpolation and returns new string.
  
  If arg is a dict, does interpolation on every string-type key in the dict (except for those in 'skip'), using the 
    dict itself as a source of symbols (plus the global variables). Returns the dict.
    
  If set, 'ignore' is a container of symbols which will interpolate to an empty string.
  """;
  # setup lookup dictionaries based on frame and depth
  if isinstance(frame,dict):
    lookups = [ frame,PyxisImpl.Context ];
  else:
    lookups = [];
    while depth>=0 and frame:
      lookups += [ frame.f_locals,frame.f_globals ];
      frame = frame.f_back
      depth -= 1
    lookups.append(PyxisImpl.Context);
  # interpolate either a single string, or a dict recursively
  if isinstance(arg,dict):
    arg,arg0 = arg.copy(),arg;
    if arg0 is not lookups[0]:
      lookups = [arg] + lookups;
    defdict = DictProxy(lookups,ignore);
  # interpolate until things stop changing, but quit after 20 loops
    for count in range(20):
      updates = {};
      for key,value in arg.iteritems():
        # interpolate string variables
        if key not in skip and isinstance(value,str):
          defdict.ignores = [value];
          newvalue = SmartTemplate(value).safe_substitute(defdict);
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
    defdict = DictProxy(lookups,ignore);
    return SmartTemplate(str(arg)%defdict).safe_substitute(defdict);
  # all other types returned as-is
  else:
    return arg;

# RE pattern matching the [PREFIX<][NAMESPACES.]NAME[?DEFAULT][:BASE][>SUFFIX] syntax
_substpattern = \
  "(?i)((?P<prefix>[^{}]+)<)?(?P<name>[._a-z][._a-z0-9]*)(\\?(?P<defval>[^}\\$]*?))?(:(?P<command>BASE|DIR))?(>(?P<suffix>[^{}]+))?"
    
class SmartTemplate (string.Template):
  pattern = "(?P<escaped>\\$\\$)|(\\$(?P<named>[_a-z][_a-z0-9]*))|(\\${(?P<braced>%s)})|(?P<invalid>\\$)"%_substpattern;

class DictProxy (object):
  itempattern = re.compile(_substpattern+"$");
  
  def __init__ (self,dicts,ignores):
    self.dicts = dicts;
    self.ignores = frozenset(ignores);
    
  def __getitem__ (self,item):
    # parse the item as a [PREFIX<]NAME[?DEFAULT][:BASE][>suffix] combo
    match = DictProxy.itempattern.match(item);
    if not match:
      return ""#,item;
    prefix,name,defval,command,suffix = match.group("prefix","name","defval","command","suffix");
    if name in self.ignores:
      return "";
    # is there an explicit namespace? otherwise use the default "merged" one
    if '.' in name:
      namespace,name = name.rsplit(".",1);
      namespace = _namespaces.get(namespace);
      if namespace is None:
        return ""; #,item;
      value = namespace.get(name,defval);
    else:
      for dd in self.dicts:
        value = dd.get(name);
        if value is not None:
          break;
    if value is None:
      value = defval;
    # check for commands
    if isinstance(value,str) and command:
      if command.upper() == "BASE":
        while value and value[-1] == "/":
          value = value[:-1];
        value = os.path.splitext(value)[0];
      elif command.upper() == "DIR":
        value = os.path.dirname(value) or ".";
    return (prefix or "")+str(value)+(suffix or "") if value not in ('',None) else "";
    
  def __contains__ (self,item):
    return True;
    
def assign (name,value,namespace=None,interpolate=True,frame=None,verbose_level=2):
  frame = frame or inspect.currentframe().f_back;
  # find namespace
  if not namespace:
    if '.' in name:
      nsname,name = name.rsplit(".",1);
      namespace = _namespaces.get(nsname);
      if namespace is None:
        raise ValueError,"invalid namespace %s"%nsname;
    else:
      namespace = frame.f_globals;
  # get superglobals associated with this namespace
  superglobs = _superglobals.get(id(namespace),[]);
  # templates are not interpolated, all others are
  if name.endswith("_Template"):
    _verbose(verbose_level,"setting %s=%s"%(name,value));
    namespace[name] = value;
  else:
    value1 = PyxisImpl.Internals.interpolate(value,frame) if interpolate else value;
    namespace[name] = value1;
    _verbose(verbose_level,"setting %s=%s%s"%(name,value1,(" (%s)"%value) if str(value) != str(value1) else ""));
  # if variable is superglobal, propagate across all namespaces defining that superglobal
  if name in superglobs:
    for ns in _namespaces.itervalues():
      if name in _superglobals.get(id(ns)):
        ns[name] = value;
  # reprocess templates
  assign_templates();

_in_assign_templates = False;

def assign_templates ():
  """For every variable in PyxisImpl.Context that ends with "_Template", assigns value to it by interpolating the template.""";
  ## non-reentrant, otherwise any call to Pyxis methods from within a template is liable to cause recursion
  global _in_assign_templates;
  if _in_assign_templates:
    return;
  _in_assign_templates = True;
  for count in range(100):
    updated = False;
    for modname,context in list(_namespaces.iteritems())+[ ("",PyxisImpl.Context) ]:
      superglobs = _superglobals[id(context)];
      newvalues = {};
      # interpolate new values for each variable that has a _Template equivalent
      for var,value in list(context.iteritems()):
        if var.endswith("_Template"):
          # string templates are interpolated, callable ones are called
          if isinstance(value,str):
            newvalues[var[:-len("_Template")]] = interpolate(value,context);
          elif callable(value):
            newvalues[var[:-len("_Template")]] = value();
      # update dict
      for var,value in newvalues.iteritems():
        oldval = context.get(var,None);
        if oldval != value:
          updated = True;
          context[var] = value;
          # if variable is superglobal, propagate across all namespaces defining that superglobal
          if var in superglobs:
            for ns in _namespaces.itervalues():
              if var in _superglobals.get(id(ns)):
                ns[var] = value;
          _verbose(3,"%s templated value %s.%s=%s"%("initialized" if oldval is None else "updated",modname,var,value));
    if not updated:
      break;
  else:
    _abort("Too many template assignment steps. This can be caused by templates that cross-reference each other");
  set_logfile(PyxisImpl.Context.get('LOG',None));
  _in_assign_templates = False;

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

_initconf_done = False;        
def initconf (force=False,*files):
  """Loads configuration from specified files, or from default file""";
  if not force and not PyxisImpl.Context.get("PYXIS_LOAD_CONFIG",True):
    return;
  global _initconf_done;
  if not _initconf_done:
    _initconf_done = True;
  if not files:
    files = glob.glob("pyxis*.py") + glob.glob("pyxis*.conf");
  # load config files -- all variable assignments go into the PyxisImpl.Context scope
  if files:
    if force:
      _verbose(1,"auto-loading config files and scripts from 'pyxis*.{py,conf}'");
    else:
      _verbose(1,"auto-loading config files and scripts from 'pyxis*.{py,conf}'. Preset PYXIS_LOAD_CONFIG=False to disable.");
  for filename in files:
    PyxisImpl.Commands.loadconf(filename,inspect.currentframe().f_back);
  assign_templates();
  _report_symbols("global",[],
      [ (name,obj) for name,obj in PyxisImpl.Context.iteritems() if not name.startswith("_") and name not in PyxisImpl.Commands.__dict__ ]);
  # make all auto-imported Pyxides modules available to global context
  toplevel = [ m for m in _modules if not '.' in m and (m in sys.modules or "Pyxides."+m in sys.modules) ];
  if PyxisImpl.Context.get("PYXIS_AUTO_IMPORT_MODULES",True) and toplevel:
    _verbose(1,"importing top-level modules (%s) for you. Preset PYXIS_AUTO_IMPORT_MODULES=False to disable."%", ".join(toplevel));
    for mod in toplevel:
      PyxisImpl.Context[mod] = sys.modules.get(mod,sys.modules.get("Pyxides."+m));
  
def loadconf (filename,frame=None):
  """Loads config file""";
  filename = interpolate(filename,frame or inspect.currentframe().f_back);
  _verbose(1,"loading %s"%filename);
  load_package(os.path.splitext(os.path.basename(filename))[0],filename);

def load_package (pkgname,filename,report=True):
  """Loads 'package' file into the Context namespace and reports on new global symbols"""
#  oldstuff = PyxisImpl.Context.copy();
  try:
    exec(file(filename),PyxisImpl.Context);
  except:
    traceback.print_exc();
    _abort("PYXIS: error parsing %s, see output above for details"%filename);
#  newnames =  [ (name,obj) for name,obj in PyxisImpl.Context.iteritems() 
#                 if not name.startswith("_") and not name in oldstuff ];
#  _report_symbols(pkgname,newnames);
  
def register_pyxis_module (superglobals=""):
  """Registers a module (the callee) as part of the Pyxis environment""";
  import PyxisImpl.Commands
  frame = inspect.currentframe().f_back;
  globs = frame.f_globals;
  # check for double registration
  if id(globs) in _superglobals:
    raise RuntimeError,"module '%s' is already registered"%modname;
  modname = globs['__name__'];
  _modules.add(modname);
  if modname.startswith("Pyxides."):
    modname = modname.split(".",1)[-1];
  # build list of superglobals
  if isinstance(superglobals,str):
    superglobs = superglobals.split();
  else:
    superglobs = itertools.chain(*[ x.split() for x in superglobals ]);
  superglobs = frozenset(superglobs);
  _verbose(1,"registered module '%s'"%modname);
  _namespaces[modname] = globs;
  _superglobals[id(globs)] = superglobs;
  # add superglobals
  for sym in superglobs:
    # if superglobal is already defined, copy its value to the new module
    # if superglobal was not yet defined, get its value form the module (or use None),
    # and propagate this value super-globally via assign
    if sym in PyxisImpl.Context:
      globs[sym] = PyxisImpl.Context[sym];
    else:
      assign(sym,globs.get(sym,None),namespace=PyxisImpl.Context,frame=frame)
  # report 
  _report_symbols(modname,superglobs,
      [ (name,obj) for name,obj in globs.iteritems() if not name.startswith("_") and name not in PyxisImpl.Commands.__dict__ and name not in superglobs ]);
  
def _report_symbols (pkgname,superglobs,syms):
  if PyxisImpl.Context['VERBOSE'] >= 2:
    varibs = sorted([name for name,obj in syms if not callable(obj) and not inspect.ismodule(obj) 
                                                              and not name.endswith("_Template") ]);
    funcs = sorted([name for name,obj in syms if callable(obj) and not name.endswith("_Template") and not isinstance(obj,ShellExecutor) ]);
    shtools = sorted([name for name,obj in syms if callable(obj) and not name.endswith("_Template") and isinstance(obj,ShellExecutor) ]);
    temps = sorted([name[:-9] for name,obj in syms if name.endswith("_Template") ]);
    if funcs:
      _verbose(2,"%s functions:"%pkgname," ".join(funcs));
    if shtools:
      _verbose(2,"%s external tools:"%pkgname," ".join(shtools));
    if superglobs:
      _verbose(2,"%s superglobals:"%pkgname," ".join(superglobs));
    if varibs:
      _verbose(2,"%s variables:"%pkgname," ".join(varibs));
    if temps:
      _verbose(2,"%s templates for:"%pkgname," ".join(temps));
  
  
    
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
  stdin =  kws.pop('stdin',None) or sys.stdin;
  stdout = kws.pop('stdout',None) or sys.stdout;
  stderr = kws.pop('stderr',None) or sys.stderr;
  # default is to split each argument at whitespace, but split_args=False passes them as-is
  split = kws.pop('split_args',True);
  # build list of arguments
  args1 = [path];
  if split:
    for arg in args:
      args1 += arg.split(" ") if isinstance(arg,str) else [ str(arg) ];
    # eliminate empty strings
    args1 = [ x for x in args1 if x ];
  else:
    args1 += map(str,args);
  # eliminate empty strings when splitting
  args = args1+["%s=%s"%(a,b) for a,b in kws.iteritems()];
  # run command
  args = [ x for x in args if x ];
  flush_log();
  if bg:
    global _bg_processes;
    po = subprocess.Popen(args,preexec_fn=_on_parent_exit('SIGTERM'));
    _bg_processes.append(po);
    _verbose(1,"executing '%s' in background: pid %d"%(" ".join(args),po.pid));
  else:
    _verbose(1,"executing '%s':"%(" ".join(args)));
    type(sys.stdout) is file and sys.stdout.flush();
    type(sys.stderr) is file and sys.stderr.flush();
    if type(stdout) is not file or type(stderr) is not file:
      po = subprocess.Popen(args,preexec_fn=_on_parent_exit('SIGTERM'));
      #retcode = subprocess.call(args);
    else:
      po = subprocess.Popen(args,preexec_fn=_on_parent_exit('SIGTERM'),stdin=stdin,stdout=stdout,stderr=stderr);
      #retcode = subprocess.call(args,stdin=stdin,stdout=stdout,stderr=stderr);
    po.wait();
    if po.returncode:
      if allow_fail:
        _warn("PYXIS: '%s' returns error code %d"%(path,po.returncode));
        return;
      else:
        _abort("PYXIS: '%s' returns error code %d"%(path,po.returncode));
    else:
      _verbose(2,"'%s' succeeded"%path);


## Function to ensure that child processes are killed when Pyxis is killed, see subprocess.Popen calls above
## source: http://www.evans.io/posts/killing-child-processes-on-parent-exit-prctl/
import signal
from ctypes import cdll

def _on_parent_exit(signame):
    """Return a function to be run in a child process which will trigger
    SIGNAME to be sent when the parent process dies.""";
    signum = getattr(signal, signame)
    # Constant taken from http://linux.die.net/include/linux/prctl.h
    _PR_SET_PDEATHSIG = 1
    def set_parent_exit_signal():
        # http://linux.die.net/man/2/prctl
        result = cdll['libc.so.6'].prctl(_PR_SET_PDEATHSIG, signum)
        if result != 0:
            raise RuntimeError('prctl failed with error code %s' % result)
    return set_parent_exit_signal
    

def find_command (comname,frame=None):
  """Locates command by name. If command is present (as a callable) in PyxisImpl.Context, returns that.
  Otherwise checks the path for a binary by that name, and returns a callable to call that.
  Else aborts.""";
  comname = interpolate(comname,frame or inspect.currentframe().f_back);
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

_re_assign = re.compile("^([\w.]+)(=)(.*)$");
_re_command1 = re.compile("^(\\??\w+)\\[(.*)\\]$");
_re_command2 = re.compile("^(\\??\w+)\\((.*)\\)$");
  
def _parse_cmdline_value (value):
  """Helper function. Parses value as a python expression. If not successful, uses a string directly"""
  try:
    exec("a=%s"%value);
    return a;
  except:
    return value;
  
def run (*commands):
  """Runs list of commands""";
  import PyxisImpl.Commands
  assign_templates();
  # _debug("running",commands);
  frame = inspect.currentframe().f_back;
  for step,command in enumerate(commands):
    # set step counter
    # interpolate the command
    command = command.strip();
    _verbose(1,"executing command %s"%command);
    # syntax 1: VAR=VALUE or VAR:=VALUE
    match = _re_assign.match(command);
    if match:
      name,op,value = match.groups();
      # assign variable -- note that templates are not interpolated
      PyxisImpl.Commands.assign(name,_parse_cmdline_value(value),frame=frame);
      continue;
    # syntax 2: command(args) or command[args]. command can have a "?" prefix
    match = _re_command1.match(command) or _re_command2.match(command);
    if match:
      comname,comargs = match.groups();
      # split up arguments
      args = [];
      kws = {};
      for arg in re.split(",," if comargs.find(",,") >=0 else ",",comargs):
        arg = interpolate(arg,frame).strip();
        match = re.match("^(\w+)=(.*)$",arg);
        if match:
          kws[match.group(1)] = _parse_cmdline_value(match.group(2));
        else:
          args.append(_parse_cmdline_value(arg));
      # exec command
      _initconf_done or initconf(force=True);  # make sure config is loaded
      comcall = find_command(comname,inspect.currentframe().f_back);
      comcall(*args,**kws);
      assign_templates();
      continue;
    # syntax 3: standalone command. This better be found!
    _initconf_done or initconf();  # make sure config is loaded
    comcall = find_command(command,inspect.currentframe().f_back);
    comcall();
    assign_templates();
  
