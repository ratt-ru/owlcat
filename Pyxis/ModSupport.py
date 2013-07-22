"""Pyxis.ModSupport: functions for programming Pyxides modules"""

import fnmatch    

import Pyxis
from Pyxis import *

from Pyxis.Commands import _verbose,_warn,_abort,makedir
from Pyxis.Internals import _superglobals,_namespaces,_modules

  
def register_pyxis_module (superglobals=""):
  """Registers a module (the callee) as part of the Pyxis environment.
  'superglobals' can be a list of super-global variables defined by the module.
  Superglobals are propagated across all modules that register for them.
  It is also possible to generate them as v.define() instead.""";
  frame = inspect.currentframe().f_back;
  globs = frame.f_globals;
  modname = globs['__name__'];
  module = sys.modules[modname];
  if modname.startswith("Pyxides."):
    modname = modname.split(".",1)[-1];
  # check for double registration
  if id(globs) in _superglobals:
    if _modules[modname] is not module:
      raise RuntimeError,"a different Pyxis module named '%s' is already registered"%modname;
  _modules[modname] = module;
  # build list of superglobals
  if isinstance(superglobals,str):
    superglobs = superglobals.split();
  else:
    superglobs = itertools.chain(*[ x.split() for x in superglobals ]);
  superglobs = set(superglobs);
  _verbose(1,"registered module '%s'"%modname);
  _namespaces[modname] = globs;
  _superglobals[id(globs)] = superglobs;
  # add superglobals
  for sym in superglobs:
    # if superglobal is already defined, copy its value to the new module
    # if superglobal was not yet defined, get its value form the module (or use None),
    # and propagate this value super-globally via assign
    if sym in Pyxis.Context:
      globs[sym] = Pyxis.Context[sym];
    else:
      assign(sym,globs.get(sym,None),namespace=Pyxis.Context,frame=frame)
  # report 
  Pyxis.Internals.report_symbols(modname,superglobs,
      [ (name,obj) for name,obj in globs.iteritems() if not name.startswith("_") and name not in Pyxis.Commands.__dict__ and name not in superglobs ]);
      

def def_global (name,default,doc=None,config=False):
  """Defines a module global with the given name, default value and documentation string.
  Mainly useful in Pyxides modules, to provide documentation on their globals.
  At module level, calling def_global('NAME',value,doc) is equivalent to
      NAME = value
      _doc_NAME = doc
  Also, if config=True, then the global variable refers to a config file that will be automatically
  backed up to OUTDIR.
  """;
  globs = inspect.currentframe().f_back.f_globals;
  globs[name] = default;
  if name.endswith("_Template"):
    name = name[:-9];
  globs.setdefault("_symdocs",{})[name] = doc;
  if config:
    globs.setdefault("_config_files",[]).append(name);
  
define = def_global;


import itertools
  
def document_globals (obj,*patterns):
  """Updates an object's documentation string with a list of globals that match the specified patterns.
  Can take multiple pattern arguments; each argument will be split at whitespace.
  Mainly useful in Pyxides modules, to form documentation strings for functions, e.g.:
  
    def_global(FOO_A,"1","option A for function foo");
    def_global(FOO_B,"2","option B for function foo");
    
    def foo ():
      # uses FOO_A and FOO_B
    
    document_globals(foo,"FOO_*");
    
  This results in foo.__doc__ being extended with documentation for FOO_A and FOO_B
  """
  patterns = list(itertools.chain(*[ patt.split() for patt in patterns ]));
  globs0 = inspect.currentframe().f_back.f_globals;
  modname0 = globs0['__name__'];
  globdocs = Pyxis.Context.get('_symdocs',{});
  # for each module, make a set of all global symbols, plus symbols that may not yet be defined, but already have documentation
  # (i.e. have entries in _symdocs)
  # each entry in allsyms is a pair of set_of_symbols,dict_of_docs
  allsyms = {}
  for modname,globs in _namespaces.iteritems():
    docs = globs.get('_symdocs',{});
    allsyms[modname] = set(itertools.chain(globs.iterkeys(),docs.iterkeys())),docs;
  if modname0 not in allsyms:
    docs = globs0.get('_symdocs',{});
    allsyms[modname0] = set(itertools.chain(globs0.iterkeys(),docs.iterkeys())),docs;
  # keep track of what's been documented, to avoid duplicates
  documented = set();
  doclist = [];
  for patt in patterns:
    # if patt matches a superglobal defined in calling module, add to list as a superglobal
    if Pyxis.Internals.is_superglobal(globs0,patt):
      doclist += [ ("v."+patt,globdocs.get(patt,'')) ];
    # else look for matching module-level globals
    else:
      if '.' in patt:
        modname,patt = patt.split('.',1);
        if modname not in allsyms:
          raise ValueError,"document_globals(\"%s.%s\"): module '%s' not registered"%(modname,patt,modname);
      else:
        modname = modname0;
      modsyms,moddocs = allsyms[modname];
      # make sorted list of matching module globals for the pattern, and loop over them
      for sym in sorted([ sym for sym in modsyms if fnmatch.fnmatch(sym,patt) and not sym.endswith("_Template") ]):
        fqsym = "%s.%s"%(modname,sym);
        # add to doclist, unless already documented
        if fqsym not in documented:
          doclist.append((fqsym,moddocs.get(sym,'')));
        documented.add(fqsym);
  # now generate documentation string
  if doclist:
    text = "\nThe following variables also apply:\n\n"
    for sym,doc in doclist:
      text += "  %-20s %s\n"%(sym,doc);
    obj.__doc__ = (obj.__doc__ or "")+ text;
  
def interpolate_locals (*varnames):
  """interpolates the variable names (from the local context) given by its argument(s).
  Returns new values in the order given. Useful as the opening line of a function, for example:
  
  def f(a="$A",b="$a.1"):
    a,b = interpolate_locals("a b")   # can also be specified as "a","b"
    
  will assign the global variable A to a, and the value of a (in this case, the value of A) plus ".1" to b.
  """;
  ## NB: the rationale for implementing it like this, as opposed to directly manipulating f_locals
  ## of the caller frame, is because f_locals of the caller can be read-only depending on Python implementation.
  Pyxis.Internals.assign_templates();
  # interpolate the whole locals() dict
  frame = inspect.currentframe().f_back;
  locs = interpolate(frame.f_locals,frame,depth=2);
  # return variables in the order listed
  ret = [ locs.get(name) for name in itertools.chain(*[ v.split(" ") for v in varnames ]) ];
  return ret if len(ret) != 1 else ret[0];
    
