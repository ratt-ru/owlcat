from Pyxis import *
import os.path

## find the Cattery
import Timba
_cattery_path = Timba.packages()['Cattery'][0]
sys.path.append(_cattery_path);

## default multithread setting
MULTITHREAD = 2

## default TDL config file
TDLCONFIG = "tdlconf.profiles"

## extra TDL options applied to all scripts
EXTRA_TDLOPTS = ""

## pipeliner tool
pipeliner = x.time.args("meqtree-pipeliner.py");

def run (script,job,config="$TDLCONFIG",section=None,args=[],options={}):
  """compiles the specified script and runs the specified job""";
  script,config,section = interpolate_locals("script config section");
  section = section or os.path.splitext(os.path.basename(script))[0];
  # run pipeliner
  pipeliner(*(
    [ "--mt $MULTITHREAD" if MULTITHREAD > 1 else "" ] +
    [ "-c $config [$section]" ] + list(args) + [ "%s=%s"%(a,b) for a,b in options.iteritems() ]+
    [ "$EXTRA_TDLOPTS $script =$job" ]
  ));

  
# register ourselves with Pyxis
register_pyxis_module();