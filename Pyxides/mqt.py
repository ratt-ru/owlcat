import os.path

# register ourselves with Pyxis
from Pyxis.ModSupport import *
register_pyxis_module();

## find the Cattery
import Timba
import os.path
_cattery_path = Timba.packages()['Cattery'][0]
sys.path.append(_cattery_path);
if v("ADD_PYXIDES_PATH",True):
  path = os.path.join(_cattery_path,"Pyxides");
  verbose(2,"adding %s to import path. Set ADD_PYXIDES_PATH=False to disable"%path);
  sys.path.append(path);

def_global('CATTERY',_cattery_path,"default path to Cattery scripts");

## default multithread setting
def_global('MULTITHREAD',2,"max number of meqserver threads");

## default TDL config file

## extra TDL options applied to all scripts
def_global('EXTRA_TDLOPTS',"","extra options passed to all TDL scripts");

## pipeliner tool
pipeliner = x.time.args("meqtree-pipeliner.py");

def_global("SCRIPT","default TDL script");
def_global("JOB","default TDL job to run");
def_global("SECTION","default section to use in TDL config file");
def_global("TDLCONFIG","tdlconf.profiles","default TDL config file",config=True);

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

document_globals(run,"MULTITHREAD EXTRA_TDLOPTS SCRIPT JOB SECTION TDLCONFIG");
