#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#% $Id$ 
#
#
# Copyright (C) 2002-2011
# The MeqTree Foundation & 
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

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

def call_system (command,args):
  rc = os.system("%s %s"%(command,args));
  if rc:
    stat = rc>>8;
    sig = rc&0xFF;
    if sig:
      print "%s process killed by signal %d"%(command,sig);
    else:
      print "%s call failed with exit status",(command,stat);
      if stat == 127:
        print "Please check your MeqTrees installation for the %s utility!"%command;
  return rc;
  


if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser
  parser = OptionParser(usage="""%prog: [options] MS""",
      description=\
      "Sets up tiled storage managers for the DATA, FLAG and BITFLAG columns of an MS, "+
      "and sets a fixed shape on the columns. This can speed up access to some MSs, e.g. "+
      "such as those produced by the uvfits2ms converter. "
  );
  parser.add_option("-f","--force",action="store_true",
                    help="proceed non-interactively");

  (options,msnames) = parser.parse_args();

  if len(msnames) != 1:
    parser.error("MS name not supplied.");
    
  if not Owlcat.find_exec("addtiledmscol"):
    print "Can't find the 'addtiledmscol' utility, which is required by this script."
    print "You may need to install or build the MeqTrees package.";
    sys.exit(1);

  msname = msnames[0];
  ms = table(msname,readonly=False);

  if not options.force:
    print """
  This script will reinitialize the DATA, FLAG and BITFLAG columns
  of the MS using a tiled storage manager and a fixed shape. It will attempt
  to preserve the data in these columns. The MODEL_DATA, CORRECTED_DATA and
  IMAGING_WEIGHT columns will be reset and their contents discarded. If this
  operation is interrupted, the MS may be left in a corrupted state, so
  please make sure you have a backup of the MS.
  """;
    if raw_input("Proceed (y/n)? ").strip().upper()[0] != "Y":
      print "Aborted by user.";
      sys.exit(1);

  SAVECOLS = [ "DATA","FLAG","BITFLAG" ];
  KILLCOLS = SAVECOLS + [ "MODEL_DATA","CORRECTED_DATA","IMAGING_WEIGHT" ];

  # read in DATA and FLAG cols
  savedcols = {};
  for colname in "DATA","FLAG","BITFLAG":
    if colname in ms.colnames():
      progress("Reading %s column"%colname);
      savedcols[colname] = ms.getcol(colname);
    elif colname is "DATA":
      print "Error: this MS does not appear to have a DATA column";
      sys.exit(1);

  datacol = savedcols['DATA'];
  progress("DATA column shape is %s"%(datacol.shape,));
  flagcol = savedcols.get('FLAG');
  if flagcol is not None and datacol.shape != flagcol.shape:
    print "Error: FLAG shape (%s) does not match DATA shape"%(flagcol.shape,);

  # kill columns
  cols = [ col for col in KILLCOLS if col in ms.colnames() ];
  if cols:
    progress("Removing columns %s"%", ".join(cols));
    ms.removecols(cols);

  # now reinsert DATA and FLAG columns
  ms.close();
  ncorr = datacol.shape[-1];
  nfreq = datacol.shape[-2];
  tilerow = 512;
  for colname,coltype in ("DATA","complex"),("FLAG","bool"),("BITFLAG","int"):
    if colname in savedcols:
      progress("Inserting tiled %s column"%colname);
      if not call_system("addtiledmscol","%s %s %s %d %d %d %d %d 4"%
                          (msname,colname,coltype,ncorr,nfreq,min(ncorr,4),min(nfreq,8),tilerow)):
        sys.exit(1);

  # reopen table
  ms = table(msname,readonly=False);
  for colname in SAVECOLS:
    col = savedcols.get(colname);
    if col is not None:
      progress("Rewriting %s column"%colname);
      ms.putcol(colname,col);

  progress("Closing MS");
  ms.close();
  progress("Adding imaging columns");
  pyrap.tables.addImagingColumns(msname);
  progress("Done");
