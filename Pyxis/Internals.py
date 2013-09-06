"""Pyxis.Internals: various internal machinery of Pyxis"""

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
import time
import fnmatch
import shutil

import Pyxis

_preset_doc = """The following global variables control Pyxis behaviour. These may be set from the command line 
as VAR=VALUE, or specified in Pyxis recipes or config files or recipes (as v.VAR = VALUE).

LOG: all output will be logged to this file (except for level-0 status messages, which are duplicated
to the console). Can also be specified as a template (e.g. LOG_Template={$MS:BASE}.log) to make the log
dependent on other variables.

VERBOSE: Pyxis verbosity level, default is 0. Higher is more verbose.

OUTDIR: output directory, for use in recipes and config files (note that Pyxis itself in no way enforces
any output to this directory -- it is up to recipes and configs to use $OUTDIR consistently in their naming
schemes filenames).

PYXIS_LOAD_CONFIG: pre-load recipes (pyxis-*.py) and config (pyxis-*.conf) files from current directory when
Pyxis is initialized. Default is true.

PYXIS_AUTO_IMPORT_MODULES: automatically import (into the global namespace) all top-level Pyxides modules 
directly or indirectly invoked by the pre-loaded recipes. Useful for interactive sessions. Default is True.

JOBS: split out up to this many subprocesses to work in parallel, when executing per() commands. 
Default is 1.

JOB_STAGGER: stagger launch of subprocesses by this many seconds. Can be useful to e.g. de-syncronize
disk access in subprocesses. Default is 10.
"""

def init (context):
  """init internals, attach to the given context""";
  global _debug
  global _info
  global _abort
  global _verbose
  global _warn
  from Pyxis.Commands import _debug,_info,_abort,_verbose,_warn
  # set default output dir
  context.setdefault("OUTDIR",".");
  context.setdefault("JOBS",0);
  context.setdefault("JOB_STAGGER",10);
  context.setdefault("PYXIS_LOAD_CONFIG",True);
  context.setdefault("PYXIS_AUTO_IMPORT_MODULES",True);
  # set default verbosity to 1
  preset_verbosity = context.get("VERBOSE",None);
  context.setdefault("VERBOSE",1);
  # set context and init stuff
  Pyxis.Context = context;
  Pyxis._predefined_names = set(context.iterkeys());
  Pyxis.Commands._init(context);
  # loaded modules
  global _namespaces,_superglobals,_modules;
  _namespaces = dict();
  _superglobals = dict();
  _modules = dict();
  # The "v" and "" namespaces correspond to the global context
  # All variables of that context are superglobals.
  _namespaces['v'] = context;
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
    
def interpolate_args (args,kws,frame,convert_lists=False): 
  """Helper function to interpolate argument list and keywords using the local dictionary, plus Pyxis.Context globals""";
  return [ interpolate(arg,frame,convert_lists=convert_lists) for arg in args ], \
         dict([ (kw,interpolate(arg,frame,convert_lists=convert_lists)) for kw,arg in kws.iteritems() ]);
  
class ShellExecutor (object):    
  """This is a ShellExecutor object, which is associated with a particular shell command, and can be
used to run that command. The command (and any fixed arguments) are specified in the constructor of 
the object. ShellExecutors are typically created via the Pyxis x, xo or xz built-ins, e.g. as 
      ls = x.ls 
      lsl = x.ls("-l")
      ls()             # runs "ls"
      dir="./test"
      lsl("$dir")      # runs "ls -l ./test"
      ls(x=1)          # runs "ls x=1"
  """;
  
  def __init__ (self,name,path,frame,allow_fail=False,get_output=False,bg=False,verbose=1,*args,**kws):
    self.name,self.path = name,path;
    self.allow_fail = allow_fail;
    self.bg  = bg;
    self.argframe = frame;
    self.get_output = get_output;
    self.verbose = verbose;
    doc = kws.pop('doc',None);
    self._add_args,self._add_kws = list(args),kws;
    self.__doc__ = "(Shell command: %s)"%str(self);
    if doc is not None:
      self.doc(doc);
      
  def doc (self,doc,append=True):
    """Sets doc string""";
    if append:
      self.__doc__ += "\n%s"%doc;
    else:
      self.__doc__ = doc;
    
  def args (self,*args,**kws):
    """Creates instance of executor with additional args. Local variables of caller are interpolated."""
    args0,kws0 = interpolate_args(self._add_args,self._add_kws,self.argframe);
    kws0.update(kws);
    return ShellExecutor(self.name,self.path,inspect.currentframe().f_back,self.allow_fail,
        self.get_output,self.bg,self.verbose,*(args0+list(args)),**kws0);
    
  def __str__ (self):
    return " ".join([self.path or ""]+self._add_args+["%s=%s"%(a,b) for a,b in self._add_kws.iteritems()]);
    
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
      args0,kws0 = interpolate_args(self._add_args,self._add_kws,self.argframe,convert_lists=True);
      args,kws = interpolate_args(args,kws,inspect.currentframe().f_back,convert_lists=True);
      kws0.update(kws);
      return _call_exec(self.path,get_output=self.get_output,allow_fail=self.allow_fail,bg=self.bg,verbose=self.verbose,*(args0+args),**kws0);

class ShellExecutorFactory (object):
  """This object can be used to create proxies for shell commands called ShellExecutors."""
  def __init__ (self,allow_fail=False,bg=False,get_output=False,verbose=0):
    self.allow_fail = allow_fail;
    self.bg = bg;
    self.get_output = get_output;
    self.verbose = verbose;
    self.doc_proto = """The '%(name)s' built-in provides an interface for shell commands. Invoking
%(name)s.sh("command") runs a command directly; invoking %(name)s.command returns a ShellExecutor for  
"command", i.e. a Python object that can used to run the command later. For example:
    ls = %(name)s.ls                   # creates a proxy for the ls command
    ls()                        # executes ls
    ls("-l")                    # executes ls -l
    ls = %(name)s("ls")                # alternative syntax for creating a ShellExecutor.\n""";
    
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
    return ShellExecutor(command,path,None,self.allow_fail,self.get_output,self.bg,self.verbose);
    
  def __call__ (self,*args,**kws):
    """An alternative way to make ShellExecutors, e.g. as x("command arg1 arg2").
    Useful when the command contains e.g. dots or slashes, thus making the x.command syntax unsuitable."""
#    args,kws = interpolate_args(args,kws,inspect.currentframe().f_back);
    if len(args) == 1:
      args = args[0].split(" ");
    return ShellExecutor(args[0],args[0],None,allow_fail=self.allow_fail,bg=self.bg,
                verbose=self.verbose,get_output=self.get_output,*args[1:],**kws);
    
  def sh (self,*args,**kws):
    """Directly invokes the shell with a command and arguments"""
    commands,kws = interpolate_args(args,kws,inspect.currentframe().f_back,convert_lists=True);
    # run command
    _verbose(self.verbose,"executing '%s':"%(" ".join(commands)));
    flush_log();
    po = subprocess.Popen(commands,preexec_fn=_on_parent_exit('SIGTERM'),shell=True,stdout=subprocess.PIPE if self.get_output else sys.stdout,stderr=sys.stderr);
    if self.get_output:
      (output,err_output) = po.communicate();
    else:
      po.wait();
      output = None;
    if po.returncode:
      if self.allow_fail:
        _warn("PYXIS: '%s' returns error code %d"%(commands[0],po.returncode));
        return;
      else:
        _abort("PYXIS: '%s' returns error code %d"%(commands[0],po.returncode));
    else:
      _verbose(self.verbose+1,"'%s' succeeded"%commands[0]);
    return output;
      
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
    
  def __contains__ (self,name):
    return name in object.__getattribute__(self,'namespace');

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
    
  v.define('VARNAME',value,doctext): assigns to variable VARNAME, and sets a documentation string
  v.doc('VARNAME'): returns the documentation string for a variable, or '' if none
  v.doc('VARNAME',doctext): sets the documentation string for a variable
  """;
  
  def __init__ (self,namespace):
#    object.__setattr__(self,'_docs',{});
    _DictAccessor.__init__(self,namespace);
    
  def __setattr__ (self,attr,value):
    assign(attr,value,namespace=self.namespace,frame=inspect.currentframe().f_back); 
    
  def __repr__ (self):
    name = object.__getattribute__(self,'__name__');
    return "Pyxis built-in %s: global variable space. Try %s.VARNAME, or help(%s)."%(name,name,name);
  
  def define (self,name,defvalue,doc=None):
    """Defines a variable, with an optional documentation string.""";
    frame = inspect.currentframe().f_back;
    register_superglobal(frame,name if not name.endswith("_Template") else name[:-len("_Template")]);
    ns = object.__getattribute__(self,'namespace');
    if name not in ns:
      assign(name,defvalue,namespace=ns,frame=frame); 
    if doc:
      ns.setdefault('_symdocs',{})[name] = doc;
  
  def doc (self,name,text=None):
    """doc(name) gets the documentation string for a variable, or '' if none is set.
    doc(name,text) sets the documentation string for the variable.""";
    ns = object.__getattribute__(self,'namespace');
    docs = ns.get('_symdocs',{});
    if text is None:
      return docs.get(name,'');
    docs[name] = text;
    return text;
    
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

def interpolate (arg,frame,depth=1,ignore=set(),skip=set(),convert_lists=False):
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
    lookups = [ frame,Pyxis.Context ];
  else:
    lookups = [];
    while depth>=0 and frame:
      lookups += [ frame.f_locals,frame.f_globals ];
      frame = frame.f_back
      depth -= 1
    lookups.append(Pyxis.Context);
  # convert lists to strings
  if isinstance(arg,(list,tuple)) and convert_lists:
    arg = ",".join(map(str,arg));
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

# RE pattern matching the [PREFIX<][NAMESPACES.]NAME[?DEFAULT][:BASE|DIR][>SUFFIX] syntax
_substpattern = \
  "(?i)((?P<prefix>[^{}]+)<)?(?P<name>[._a-z][._a-z0-9]*)(\\?(?P<defval>[^}\\$]*?))?(:(?P<command>BASE|DIR|FILE))?(>(?P<suffix>[^{}]+))?"
    
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
    if value is None or value == '':
      value = defval;
    # check for commands
    if isinstance(value,str) and command:
      if command.upper() == "BASE":
        value = value and os.path.basename(value);
        value = value and os.path.splitext(value)[0];
      elif command.upper() == "DIR":
        value = (value and os.path.dirname(value)) or ".";
      elif command.upper() == "FILE":
        value = value and os.path.basename(value);
    return (prefix or "")+str(value)+(suffix or "") if value not in ('',None) else "";
    
  def __contains__ (self,item):
    return True;

def _resolve_namespace (name,frame,autoimport=False):
  if '.' in name:
    nsname,name = name.rsplit(".",1);
    if autoimport:
      _autoimport(nsname);
    namespace = _namespaces.get(nsname);
    if namespace is None:
      raise ValueError,"invalid namespace %s"%nsname;
  else:
    namespace = frame.f_globals;
  return namespace,name;
    
    
def assign (name,value,namespace=None,interpolate=True,frame=None,append=False,autoimport=False,verbose_level=2):
  """Assigns value to variable, then reevaluates templates etc."""
  frame = frame or inspect.currentframe().f_back;
  # find namespace
  if not namespace:
    namespace,name = _resolve_namespace(name,frame,autoimport=autoimport);
  modname = namespace.get('__name__',"???") if namespace is not Pyxis.Context else "v";
  # interpolate if asked to, unless this is a template, which are never interpolated
  if interpolate and not name.endswith("_Template"):
    value1 = Pyxis.Internals.interpolate(value,frame);
  else:
    value1 = value;
  # append if asked to
  if append:
    value0 = str(namespace.get(name,""));
    if value0:
      value1 = "%s %s"%(value0,value1);
  # get list of namespaces in which to set variable
  namespaces = [ namespace ];
#  namespace[name] = value1;
#  _verbose(verbose_level,"setting %s.%s=%s"%(modname,name,value1));
  # get superglobals associated with this namespace: if the variable is one of them, then propagate the new value
  # across all other namespaces using that superglobal
  superglobs = _superglobals.get(id(namespace),[]);
  if name in superglobs:
    namespaces += [ ns for ns in _namespaces.itervalues() if name in _superglobals.get(id(ns)) and ns is not namespace ];
  # now assign
  for ns in namespaces:
    _verbose(verbose_level,"setting %s.%s=%s"%(ns['__name__'] if ns is not Pyxis.Context else "v",name,value1));
    ns[name] = value1;
    # if assigning a template, make sure the template is re-enabled
    if name.endswith("_Template"):
      ns.get('__pyxis_template_ids',{}).pop(name[:-len("_Template")],None);
  # reprocess templates
  assign_templates();

def unset (name,namespace=None,frame=None,verbose_level=2):
  """Unsets variable. Reactivates any templates associated with variable."""
  frame = frame or inspect.currentframe().f_back;
  # find namespace
  if not namespace:
    namespace,name = _resolve_namespace(name,frame);
  # get list of namespaces from which to unset
  namespaces = [ namespace ];
  # get superglobals associated with this namespace: if the variable is one of them, then propagate the un-setting
  # across all other namespaces using that superglobal
  superglobs = _superglobals.get(id(namespace),[]);
  if name in superglobs:
    namespaces += [ ns for ns in _namespaces.itervalues() if name in _superglobals.get(id(ns)) and ns is not namespace ];
  # now loop over all namespaces in which to unset
  for ns in namespaces:
    modname = ns.get('__name__',"???") if ns is not Pyxis.Context else "v";
    if name in ns:
      _verbose(verbose_level,"unsetting %s.%s"%(modname,name,value1));
      ns.pop(name);
    tmpldict = ns.get('__pyxis_template_ids',{});
    if name in tmpldict:
      tmpldict.pop(name);
      _verbose(verbose_level,"  and re-enabling associated template %s.%s_Template"%(modname,name));
  # reprocess templates
  assign_templates();


class DependsHandler (object):
  _handlers = {};
  
  def __init__ (self,varname,vardepend,frame=None):
    # resolve namespaces
    frame = frame or inspect.currentframe().f_back;
    self.namespace,self.name = _resolve_namespace(varname,frame);
    self.depnamespace,self.depname = _resolve_namespace(vardepend,frame);
    # add to list
    self.patterns = [];
    DependsHandler._handlers.setdefault((id(self.depnamespace),self.depname),[]).append(self);
    
  def __call__ (self,pattern,value):
    self.patterns.append((pattern,value));
    return self;
    
  def default (self,value):
    self.defval = value;
    return self;
    
  @staticmethod
  def update (namespace,name,value):
    result = False;
    for hdl in DependsHandler._handlers.get((id(namespace),name),[]):
      result |= hdl.check(value);
    return result;
    
_in_assign_templates = False;

def assign_templates ():
  """For every variable in a Pyxis module (or the global context) that ends with "_Template", assigns value 
  to it by interpolating the template.""";
  ## non-reentrant, otherwise any call to Pyxis methods from within a template is liable to cause recursion
  global _in_assign_templates;
  if _in_assign_templates:
    return;
  _in_assign_templates = True;
  for count in range(100):
    updated = False;
    for modname,context in list(_namespaces.iteritems()):
      superglobs = _superglobals[id(context)];
      newvalues = {};
      templdict = context.setdefault("__pyxis_template_ids",{});
      # interpolate new values for each variable that has a _Template equivalent
      for var,value in list(context.iteritems()):
        if var.endswith("_Template"):
          # get old value of variable
          varname = var[:-len("_Template")];
          oldvalue = varvalue = context.get(varname);
          # check if template is defined in the wrong place, superglobal templates must be defined
          # in the superglobal context
          if varname in superglobs and context is not Pyxis.Context:
            _abort("%s.%s defined for superglobal v.%s. Fix your scripts please: use v.%s instead"%
                    (modname,var,varname,var));
          # check if template is still active
          if varname in templdict:
            idvar = templdict[varname];
            # check if template has been disabled by setting __pyxis_template_ids[var] = None
            if idvar is None:
              _verbose(3,"template for %s.%s has been disabled, ignoring"%(modname,var));
              continue;
            # check if variable has been explicitly assigned to since the last time the template
            # got evaluated (i.e. if the id has changed). If so, disable template
            if idvar != id(oldvalue):
              _verbose(3,"value for %s.%s has been explicitly set [%x, was %x], disabling template"%(modname,var,id(oldvalue),idvar));
              templdict[varname] = None;
              continue;
          # catch all template errors below
          try:
            # string templates are simply interpolated 
            if isinstance(value,str):
              varvalue = interpolate(value,context);
            # templates of the form
            # SELECT_EXPR,{pattern:value,pattern:value,...}  [,ELSEVALUE]
            # or SELECT_EXPR,[ (pattern,value),(pattern,value),... ]  [,ELSEVALUE]
            elif isinstance(value,tuple):
              if len(value) not in (2,3) or not isinstance(value[0],str) or not isinstance(value[1],(list,tuple,dict)):
                raise TypeError,"invalid select clause";
              select_expr = interpolate(value[0],context);
              # loop through patterns, if one matches the select expression, return that value
              for pair in ( value[1] if not isinstance(value[1],dict) else value[1].iteritems() ):
                if len(pair) != 2 or not isinstance(pair[0],str):
                  raise TypeError,"invalid element in select clause";
                if fnmatch.fnmatch(select_expr,interpolate(pair[0],context)):
                  varvalue = pair[1];
                  break;
              # no patterns match? Look for ELSEVALUE, if defined
              else:
                if len(value) == 3:
                  varvalue = value[2];
                elif isinstance(value[1],dict) and 'default' in value[1]:
                  varvalue = value[1]['default'];
            # list templates are interpolated per-element
            elif isinstance(value,list):
              varvalue = [ interpolate(val) for val in value ];
            # callable templates are called directly
            elif callable(value):
              varvalue = value();
          except:
            traceback.print_exc();
            _warn("PYXIS: error evaluating template %s"%var);
          if varvalue is not oldvalue:
            newvalues[varname] = varvalue;
      # update dict
      for var,value in newvalues.iteritems():
        oldval = context.get(var);
        if oldval != value:
          updated = True;
          # set value, and set id in templdict for later comparison
          context[var] = value;
          templdict[var] = id(value); 
          # if variable is superglobal, propagate it to global context
          if var in superglobs:
            Pyxis.Context[var] = value;
          _verbose(3,"%s templated value %s.%s=%s [%x]"%("initialized" if oldval is None else "updated",modname,var,value,id(value)));
    if not updated:
      break;
  else:
    _abort("Too many template assignment steps. This can be caused by templates that cross-reference each other");
  # propagate superglobals
  for ns in _namespaces.itervalues():
    if ns is not Pyxis.Context:
      for var in _superglobals.get(id(ns)):
        ns[var] = Pyxis.Context[var];
  # set logger, in case LOG value has changed
  set_logfile(Pyxis.Context.get('LOG',None));
  _in_assign_templates = False;

_current_logfile = None;
_current_logobj = None;
_warned_nolog = None;
_visited_logfiles = set();
_disable_logfiles = False;

def flush_log ():
  global _current_logobj,_current_logfile;
  _current_logobj and _current_logobj.flush();

def get_logfile ():
  global _current_logobj,_current_logfile;
  return _current_logobj,_current_logfile;

def disable_logfiles ():
  global _disable_logfiles;
  _disable_logfiles = True;

def set_logfile (filename,quiet=False):
  """Starts logging to the specified file""";
  global _current_logfile,_current_logobj,_warned_nolog;
  if _disable_logfiles:
    if not _warned_nolog:
      _warn("logfiles disabled, forcing Pyxis output to console and ignoring LOG assignments");
      _warned_nolog = True;
    return;
  import Pyxis.ModSupport
  if filename is not None:
    filename = str(filename);
  if filename == "-" or not filename:
    filename = None;
  if filename != _current_logfile:
    if Pyxis.Context.get('get_ipython'):
      if not _warned_nolog:
        _warn("running inside ipython, forcing Pyxis output to console and ignoring LOG assignments");
        _warned_nolog = True;
      return;
    if not quiet:
      _info("redirecting log output to %s"%(filename or "console"),console=True);
    if filename is None:
      sys.stdout,sys.stderr = sys.__stdout__,sys.__stderr__;
      _current_logobj = None;
    else:
      mode = "w";
      # append to file if name starts with +, or if file has already been used as a log this session
      if filename[0] == '+':
        filename = filename[1:];
        mode = "a";
      if filename in _visited_logfiles:
        mode = "a";
      Pyxis.ModSupport.makedir(os.path.dirname(filename),no_interpolate=True);
      _current_logobj = sys.stdout = sys.stderr = open(filename,mode);
      _visited_logfiles.add(filename);
#    if _current_logfile:
#      pass;
#      _info("log continued from %s"%_current_logfile);
#    else:
#      _info("log started");
    _current_logfile = filename;

_initconf_done = False;  
_config_files = [];

def initconf (force=False,*files):
  """Loads configuration from specified files, or from default file""";
#  print "initconf",force,Pyxis.Context.get("PYXIS_LOAD_CONFIG",True);
  if not force and not Pyxis.Context.get("PYXIS_LOAD_CONFIG",True):
    return;
  global _initconf_done;
  if not _initconf_done:
    _initconf_done = True;
  if not files:
    files = glob.glob("pyxis*.py") + glob.glob("pyxis*.conf");
  global _config_files;
  _config_files = files;
  # remember current set of globals
  oldsyms = frozenset(Pyxis.Context.iterkeys());
  # load config files -- all variable assignments go into the Pyxis.Context scope
  if files:
    if force:
      _verbose(1,"auto-loading config files and scripts from 'pyxis*.{py,conf}'");
    else:
      _verbose(1,"auto-loading config files and scripts from 'pyxis*.{py,conf}'. Preset PYXIS_LOAD_CONFIG=False to disable.");
  for filename in files:
    Pyxis.Commands.loadconf(filename,inspect.currentframe().f_back);
  assign_templates();
  # report on global symbols
  report_symbols("global",[],
      [ (name,obj) for name,obj in Pyxis.Context.iteritems() 
        if name not in oldsyms and not name.startswith("_") and name not in ("In","Out") and name not in Pyxis.Commands.__dict__ ]);
  # make all auto-imported Pyxides modules available to global context
  toplevel = [ m for m in _modules.iterkeys() if not '.' in m and (m in sys.modules or "Pyxides."+m in sys.modules) ];
  if Pyxis.Context.get("PYXIS_AUTO_IMPORT_MODULES",True) and toplevel:
    _verbose(1,"importing top-level modules (%s) for you. Preset PYXIS_AUTO_IMPORT_MODULES=False to disable."%", ".join(toplevel));
    for mod in toplevel:
      Pyxis.Context[mod] = sys.modules.get(mod,sys.modules.get("Pyxides."+m));
  
def loadconf (filename,frame=None):
  """Loads config file""";
  filename = interpolate(filename,frame or inspect.currentframe().f_back);
  _verbose(1,"loading %s"%filename);
  load_package(os.path.splitext(os.path.basename(filename))[0],filename);
  
  
def saveconf ():
  """Saves config files to OUTDIR""";
  OUTDIR = Pyxis.Context['OUTDIR'] or '.';
  # make set of all config files
  configs = set(_config_files);
  for m,globs in _namespaces.iteritems():
    for fvar in globs.get('_config_files',[]):
      if fvar in globs:
        configs.add(globs[fvar]);
  # now back them up
  for ff in configs:
    dest = os.path.join(OUTDIR,os.path.basename(ff));
    if os.path.exists(ff):
      if not os.path.exists(dest) or os.path.getmtime(dest) < os.path.getmtime(ff):
        _info("copying %s to %s"%(os.path.basename(ff),OUTDIR));
        shutil.copyfile(ff,dest);
  

def load_package (pkgname,filename,report=True):
  """Loads 'package' file into the Context namespace and reports on new global symbols"""
#  oldstuff = Pyxis.Context.copy();
  try:
    exec(file(filename),Pyxis.Context);
  except SystemExit:
    raise;
  except KeyboardInterrupt:
    raise;
  except:
    traceback.print_exc();
    _abort("PYXIS: error parsing %s, see output above and/or log for details"%filename);
#  newnames =  [ (name,obj) for name,obj in Pyxis.Context.iteritems() 
#                 if not name.startswith("_") and not name in oldstuff ];
#  report_symbols(pkgname,newnames);
  
def register_superglobal (frame,sym):
  """Helper function. Registers the given symbol as a superglobal for the module given by 'frame'.
  If this is not a registered Pyxis module, does nothing.""";
  globs = frame.f_globals;
  modname = globs['__name__'];
  _sgs = _superglobals.get(id(globs));
  if _sgs is not None:
    if sym in _sgs:
      _verbose(3,"%s already registered as a superglobal for module '%s'"%(sym,modname));
      return;
    _verbose(2,"defining superglobal %s in module '%s'"%(sym,modname));
    _sgs.add(sym);
  else:
    _verbose(3,"'%s' is not a registered module, ignoring superglobal '%s'"%(modname,sym));
  
def is_superglobal (globs,sym):
  """Returns True if symbol is registered as a superglobal in the module whose globals are glob"""
  return sym in _superglobals.get(id(globs),{});
  
def report_symbols (pkgname,superglobs,syms):
  if Pyxis.Context['VERBOSE'] >= 2:
    # remove modules from symbols
    syms = [ (name,obj) for name,obj in syms if not inspect.ismodule(obj) ];
    varibs = sorted([name for name,obj in syms if not callable(obj) and not name.endswith("_Template") ]);
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
  verbose = kws.pop('verbose',1);
  if kws.pop('get_output',None):
    stdout = subprocess.PIPE;
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
  # if command is 'time', then actual command is second argument
  cmdname = args[1] if args[0] == "time" or args[0] == "/usr/bin/time" else args[0];
  # run command
  args = [ x for x in args if x ];
  flush_log();
  if bg:
    global _bg_processes;
    po = subprocess.Popen(args,preexec_fn=_on_parent_exit('SIGTERM'));
    _bg_processes.append(po);
    _verbose(verbose,"executing '%s' in background: pid %d"%(" ".join(args),po.pid));
  else:
    _verbose(verbose,"executing '%s':"%(" ".join(args)));
    type(sys.stdout) is file and sys.stdout.flush();
    type(sys.stderr) is file and sys.stderr.flush();
#    if type(stdout) is not file or type(stderr) is not file:
#      po = subprocess.Popen(args,preexec_fn=_on_parent_exit('SIGTERM'));
#      #retcode = subprocess.call(args);
#    else:
    po = subprocess.Popen(args,preexec_fn=_on_parent_exit('SIGTERM'),stdin=stdin,stdout=stdout,stderr=stderr);
    if stdout is subprocess.PIPE:
      (output,err_output) = po.communicate();
    else:
      po.wait();
      output = None;
    if po.returncode:
      if allow_fail:
        _warn("%s returned error code %d"%(cmdname,po.returncode));
        return;
      else:
        _abort("%s returned error code %d"%(cmdname,po.returncode));
    else:
      _verbose(verbose+1,"%s succeeded"%cmdname);
    return output;


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

def _autoimport (modname):
  if modname not in _namespaces:
    _verbose(1,"auto-importing module %s"%modname);
    Pyxis.Context[modname] = __import__(modname,Pyxis.Context);
    
    
_re_mod_command = re.compile("^(\w[\w.]+)\\.(\w+)$");
    
def find_command (comname,frame=None,autoimport=True):
  """Locates command by name. 
  If command is of form "module.command", attempts to resolve it to a callable.
  Will import module if it is not found and autoimport=True.
  Else if command is present (as a callable) in Pyxis.Context, returns that.
  Otherwise checks the path for a binary by that name, and returns a callable to call that.
  Else aborts.""";
  comname = interpolate(comname,frame or inspect.currentframe().f_back);
  # try auto-import
  m = _re_mod_command.match(comname);
  if m and autoimport:
    _autoimport(m.group(1));
  # first, try to evaluate to a callable function
  try:
    comcall = eval(comname,Pyxis.Context);
    if callable(comcall):
      return comcall;
  except:
    pass;
  # failing that, look for a shell command
  if comname[0] == "?":
    comname = comname[1:];
    allow_fail = True;
  else:
    allow_fail = False;
  path = find_exec(comname);
  if path is None:
    _abort("undefined command '%s'"%comname);
  # make callable for this shell command
  return lambda *args:_call_exec(path,allow_fail=allow_fail,*args);

_re_assign = re.compile("^(\w[\w.]*)(\\+?=)(.*)$");
_re_command1 = re.compile("^(\\??\w[\w.]*)\\[(.*)\\]$");
_re_command2 = re.compile("^(\\??\w[\w.]*)\\((.*)\\)$");
  
def _parse_cmdline_value (value):
  """Helper function to parse the VALUE portion of VAR=VALUE commands.
  Evaluates the value string as a Python expression, using eval() with superglobals.
  If this fails, then returns the value directly as a string.
  
  This implies that command-line arguments will work as follows:
  
  Argument:              Corresponding Python code:
  ---------              --------------------------
  VAR=x                  VAR = "x" if superglobal x is undefined, else VAR=x
  "VAR='x'"              VAR = "x"    note that shell will swallow the outer quotes
  VAR=1                  VAR = 1
  VAR=complex(1)         VAR = 1+0j
  "VAR=dict(x='y')"      VAR = dict(x='y')
  VAR=x=1                VAR = "x=1"
  """
  try:
    return eval(value,Pyxis.Context);
  except:
    return value;

def run (*commands):
  """Runs list of commands""";
  import Pyxis.Commands
  assign_templates();
  # _debug("running",commands);
  frame = inspect.currentframe().f_back;
  for step,command in enumerate(commands):
    # if command is callable, call directly
    if not callable(command):
      # interpolate the command
      command = command.strip();
      _verbose(1,"running command %s"%command);
      # syntax 1: VAR=VALUE or VAR+=VALUE
      match = _re_assign.match(command);
      if match:
        name,op,value = match.groups();
        # assign variable -- note that templates are not interpolated
        Pyxis.Commands.assign(name,_parse_cmdline_value(value),frame=frame,append=(op=="+="),autoimport=True);
        continue;
      # syntax 2: command(args) or command[args]. command can have a "?" prefix to make it optional
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
        # if command is 'help', disable logging
        if comname == "help" and _current_logobj:
          logfile = _current_logfile;
          set_logfile(None);
        else:
          logfile = None;
        # exec command
        _initconf_done or initconf(force=True);  # make sure config is loaded
        comcall = find_command(comname,inspect.currentframe().f_back,autoimport=True);
        comcall(*args,**kws);
        assign_templates();
        # reset logging, if disabled for 'help'
        if logfile:
          set_logfile(logfile);
        continue;
      # syntax 3: standalone command. This better be found!
      _initconf_done or initconf(force=True);  # make sure config is loaded
      comcall = find_command(command,inspect.currentframe().f_back);
    # fall through here if command is callable
    else:
      comcall = command;
      _verbose(1,"running command %s"%command.__name__);
    comcall();
    assign_templates();
  
