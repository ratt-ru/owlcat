import inspect
import os.path
import math
import numpy
np = numpy

# some useful constants
DEG = math.pi/180
ARCMIN = DEG/60
ARCSEC = DEG/3600

#from Pyxis.Commands import *
import Pyxis.Internals 

_initialized = False;

if not _initialized:
  """Initializes Pyxis with a global context dict.""";
  # initialize Pyxis using the globals of whoever imported the module  
  _context = inspect.currentframe().f_back.f_globals;
  
  Pyxis.Internals.init(_context);  
  from Pyxis.Commands import *
 
  verbose(1,"===[ Pyxis: Python eXtensions for Inteferometry Scripting (C) 2013 by Oleg Smirnov <oms@ska.ac.za> ]===");

  # import basic Pyxis commands into the context
  verbose(1,"loading Pyxis into context '%s'"%_context.get('__name__'));
  exec('from Pyxis.Commands import *',_context);
  exec('from Pyxis.Commands import _I,_II',_context);
  
  import Pyxides
  # add Pyxides path to module includes (so we can do stuff like "import ms" instead of "from Pyxides import ms"
  if _context.get("ADD_PYXIDES_PATH",True):
    verbose(2,"adding Pyxides to import path. Set ADD_PYXIDES_PATH=False to disable");
    sys.path.append(os.path.dirname(Pyxides.__file__));

  ## import standard modules, unless a specific other set is given
  #if not context.get("pyxis_preload"):
    #import StandardModules
    #for mod in StandardModules.pyxis_preload:
      #filename = os.path.join(os.path.dirname(StandardModules.__file__),mod)+".py";
      #if os.path.exists(filename):
        #verbose(1,"loading standard package '%s' from %s"%(mod,filename));
        #Pyxis.Internals.load_package(mod,filename);        
      #else:
        #warn("can't find standard package %s"%filename);
  
  Pyxis.Internals.initconf();
  
