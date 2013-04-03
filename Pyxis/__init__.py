import inspect
import os.path

import PyxisImpl
import PyxisImpl.Commands
import PyxisImpl.Internals 

def init (context):
  """Initializes Pyxis with a global context dict.
  """;
  PyxisImpl.Internals.init(context);  
  from PyxisImpl.Commands import verbose,warn;
  
  verbose(1,"===[ Pyxis: Python eXtensions for Inteferometry Scripting (C) 2013 by Oleg Smirnov <oms@ska.ac.za> ]===");

  # import basic Pyxis commands into the context
  verbose(1,"loading Pyxis into context '%s'"%context.get('__name__'));
  exec('from PyxisImpl.Commands import *',context);
  exec('from PyxisImpl.Commands import _I,_II',context);

  # import standard modules, unless a specific other set is given
  if not context.get("pyxis_preload"):
    import StandardModules
    for mod in StandardModules.pyxis_preload:
      filename = os.path.join(os.path.dirname(StandardModules.__file__),mod)+".py";
      if os.path.exists(filename):
        verbose(1,"loading standard package '%s' from %s"%(mod,filename));
        PyxisImpl.Internals.load_package(mod,filename);        
      else:
        warn("can't find standard package %s"%filename);
  
  # init config, unless disabled
  if not context.get("pyxis_skip_config"):
    PyxisImpl.Internals.initconf();
  
# initialize Pyxis using the globals of whoever imported the module  
init(inspect.currentframe().f_back.f_globals);

