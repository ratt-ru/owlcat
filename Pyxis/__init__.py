import inspect
import os.path
import math
import numpy
np = numpy

# some useful constants
DEG = math.pi/180
ARCMIN = DEG/60
ARCSEC = DEG/3600

import PyxisImpl
#from PyxisImpl.Commands import *
import PyxisImpl.Internals 

_initialized = False;

if not _initialized:
  """Initializes Pyxis with a global context dict.""";
  # initialize Pyxis using the globals of whoever imported the module  
  _context = inspect.currentframe().f_back.f_globals;
  
  PyxisImpl.Internals.init(_context);  
  from PyxisImpl.Commands import *
 
  verbose(1,"===[ Pyxis: Python eXtensions for Inteferometry Scripting (C) 2013 by Oleg Smirnov <oms@ska.ac.za> ]===");

  # import basic Pyxis commands into the context
  verbose(1,"loading Pyxis into context '%s'"%_context.get('__name__'));
  exec('from PyxisImpl.Commands import *',_context);
  exec('from PyxisImpl.Commands import _I,_II',_context);
  
  import Pyxides
  # add Pyxides path to module includes (so we can do stuff like "import ms" instead of "from Pyxides import ms"
  sys.path.append(os.path.dirname(Pyxides.__file__));

  ## import standard modules, unless a specific other set is given
  #if not context.get("pyxis_preload"):
    #import StandardModules
    #for mod in StandardModules.pyxis_preload:
      #filename = os.path.join(os.path.dirname(StandardModules.__file__),mod)+".py";
      #if os.path.exists(filename):
        #verbose(1,"loading standard package '%s' from %s"%(mod,filename));
        #PyxisImpl.Internals.load_package(mod,filename);        
      #else:
        #warn("can't find standard package %s"%filename);
  
  # init config, unless disabled
  if not _context.get("pyxis_skip_config"):
    PyxisImpl.Internals.initconf();
  
