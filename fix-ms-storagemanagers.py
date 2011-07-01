#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os.path
import re
import time
import math
import Owlcat

ROWCHUNK = 50000;

time_start = time.time();

import pyrap.tables
from pyrap.tables import table,tablecopy,tableexists,tabledelete,tablecreatearraycoldesc

def timestamp (format="%H:%M:%S:"):
  return time.strftime(format,time.gmtime(time.time()-time_start));

def progress (message,newline=True):
  sys.stdout.write("%s %-60s%s"%(timestamp(),message,'\n' if newline else '\r'));
  sys.stdout.flush();

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser
  parser = OptionParser(usage="""%prog: [options] MS""",
      description=\
      """Sets up tiled storage managers for the DATA, FLAG and BITFLAG columns of an MS,
      and sets a fixed shape on the columns. This can speed up access to some MSs, e.g.
      such as those produced by the uvfits2ms converter."""
  );

  (options,msnames) = parser.parse_args();

  if len(msnames) != 1:
    parser.error("MS name not supplied.");

  msname = msnames[0];
  ms = table(msname,readonly=False);

  # read in DATA and FLAG cols
  progress("Reading DATA column");
  datacol = ms.getcol("DATA");
  progress("Read DATA column of shape %s"%(datacol.shape,));
  progress("Reading FLAG column");
  flagcol = ms.getcol("FLAG");

  if datacol.shape != flagcol.shape:
    print "Error: FLAG shape (%s) does not match DATA shape"%(flagcol.shape,);

  # kill columns
  cols = [ col for col in ["DATA","FLAG","BITFLAG","MODEL_DATA","CORRECTED_DATA","IMAGING_WEIGHT"]
           if col in ms.colnames() ];
  if cols:
    progress("Removing columns %s"%", ".join(cols));
    ms.removecols(cols);

  # now reinsert DATA and FLAG columns
  progress("Inserting tiled DATA, FLAG and BITFLAG columns");
  ms.close();
  ncorr = datacol.shape[-1];
  nfreq = datacol.shape[-2];
  tilerow = 512;
  for colname in ("DATA","FLAG","BITFLAG"):
    progress("Adding tiled %s column"%colname);
    rc = os.system("addtiledmscol %s %s complex %d %d %d %d %d 4"%
                 (msname,colname,ncorr,nfreq,min(ncorr,4),min(nfreq,8),tilerow));
    if rc:
      print "addtiledmscol process killed by signal %d"%(rc&0xFF);
      rc>>=8;
      sys.exit(rc);

  # reopen table
  ms = table(msname,readonly=False);
  progress("Rewriting DATA column");
  ms.putcol('DATA',datacol);
  progress("Rewriting FLAG column");
  ms.putcol('FLAG',flagcol);
  progress("Closing MS");
  ms.close();
  progress("Adding imaging columns");
  pyrap.tables.addImagingColumns(msname);
  progress("Done");
