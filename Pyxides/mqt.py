import os.path

# register ourselves with Pyxis
from Pyxis.ModSupport import *
register_pyxis_module();

## find the Cattery
import Timba
_cattery_path = Timba.packages()['Cattery'][0]
sys.path.append(_cattery_path);

## default multithread setting
MULTITHREAD = 2

## default TDL config file

## extra TDL options applied to all scripts
EXTRA_TDLOPTS = ""

## pipeliner tool
pipeliner = x.time.args("meqtree-pipeliner.py");

SCRIPT = ""
JOB = ""
SECTION = ""
TDLCONFIG = "tdlconf.profiles"

def run (script="$SCRIPT",job="$JOB",config="$TDLCONFIG",section="$SECTION",args=[],options={}):
  """Uses meqtree-pipeliner to compile the specified MeqTrees 'script', using 'config' file and config 'section',
  then runs the specified 'job'.
  Use a list of 'args' to pass extra arguments to meqtree-pipeliner. Use a dict of 'options' to
  pass extra arguments as key=value.""";
  script,job,config,section = interpolate_locals("script job config section");
  section = section or os.path.splitext(os.path.basename(script))[0];
  # run pipeliner
  pipeliner(*(
    [ "--mt $MULTITHREAD" if MULTITHREAD > 1 else "" ] +
    [ "-c $config [$section]" ] + list(args) + [ "%s=%s"%(a,b) for a,b in options.iteritems() ]+
    [ "$EXTRA_TDLOPTS $script =$job" ]
  ));

  
