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
      description="""Sets up proper tiled storage managers for the DATA and FLAG columns of an MS.""");
  #parser.add_option("-o","--output",dest="output",type="string",
                    #help="name of output FITS file");
  #parser.add_option("-z","--zoom",dest="output_quad",type="string",
                    #help="name of zoomed output FITS file");
  #parser.add_option("-t","--tolerance",dest="tolerance",type="float",
                    #help="How close (in meters) two baselines need to be to each other to be considered redundant (default .1)");
  #parser.add_option("-I","--ifrs",dest="ifrs",type="string",
                    #help="subset of interferometers to use.");
  #parser.add_option("-s","--select",dest="select",action="store",
                    #help="additional TaQL selection string. Note that redundant baselines are counted only within the subset "
                         #"given by the --ifrs and --select options.");
#  parser.set_defaults(tolerance=.1,select="",ifrs="");

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
  os.system("addtiledmscol %s DATA complex %d %d %d %d %d 4"%
              (msname,ncorr,nfreq,min(ncorr,4),min(nfreq,8),tilerow));
  os.system("addtiledmscol %s FLAG bool %d %d %d %d %d 4"%
              (msname,ncorr,nfreq,min(ncorr,4),min(nfreq,8),tilerow));
  os.system("addtiledmscol %s BITFLAG int %d %d %d %d %d 4"%
              (msname,ncorr,nfreq,min(ncorr,4),min(nfreq,8),tilerow));

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
