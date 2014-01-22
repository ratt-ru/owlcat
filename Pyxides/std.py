from Pyxis.ModSupport import *

import tempfile

# register ourselves with Pyxis, and define the superglobals
register_pyxis_module();

v.define("OUTDIR","",
  """base output directory""");
v.define("SUFFIX_Template",'${spw<ms.DDID}',
  """suffix added to filenames, default is "-spwX""");
v.define("DESTDIR_Template",'${OUTDIR>/}plots${MS:BASE>-}${spw<ms.DDID}',
  """destination directory for plots, images, etc.""");
v.define("OUTFILE_Template",'${DESTDIR>/}${MS:BASE>-}${SUFFIX>-}${s<STEP>-}${LABEL}',
  """base output filename for plots, images, etc.""");
v.define("STEP",1,
  """step counter, automatically incremented. Useful for decorating filenames.""")
v.define("LABEL","",
  """decorative label, mainly used for decorating filenames.""")
  
remove          = xo.rm.args("-fr");
copy            = x.cp.args("-a");
plotparms       = x("plot-parms.py").args("$PLOTPARMS_ARGS");
fitstool        = x("fitstool.py");

casapy = x.casapy.args("--nologger --log2term -c");

def runcasapy (command):
  command = interpolate_locals("command");
  # write command to script file
  tf = tempfile.NamedTemporaryFile(suffix=".py");
  tf.write(command+"\nexit\n");
  tf.flush();
  tfname = tf.name;
  # run casapy
  info("Running casapy $tfname. Content:\n$command\n");
  retcode = casapy(tfname);
  tf.close();
  if retcode:
    abort("casapy failed with return code %d. Check the logs for errors.");

def printpaths ():
  info("OUTDIR=$OUTDIR, DESTDIR=$DESTDIR, OUTFILE=$OUTFILE");
