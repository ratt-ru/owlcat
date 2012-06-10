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

import os.path
import sys
import traceback
import gzip
import cPickle
import Owlcat

flagger = parser = ms = msname = None;

def error (message):
  print "%s: %s"%(os.path.basename(sys.argv[0]),message);
  sys.exit(1);

def get_ms ():
  global ms;
  global msname;
  if not ms:
    ms = Owlcat.table(msname);
  return ms;

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser,OptionGroup
  parser = OptionParser(usage="""%prog: [actions] [options] MS COLUMN FILENAME[.gz]""",
      description="Exports MS column to an external file, which can be reloaded with import-ms-column."
  );
  parser.add_option("-r","--row-step",type="int",
                    help="how many rows to step over at a time, default is %default");
#  parser.add_option("-z","--gzip",action="store_true",
#                    help="compress data with gzip");
#  parser.add_option("-v","--verbose",metavar="LEVEL",type="int",
#                    help="verbosity level for messages. Higher is more verbose, default is 0.");
  parser.set_defaults(verbose=0,row_step=200000);

  # parse args
  (options,args) = parser.parse_args();
  if len(args) != 3:
    parser.error("Incorrect number of arguments. Use '-h' for help.");
  msname,colname,filename = args;

  try:
    gzf = gzip.GzipFile(filename,"w") if filename.endswith(".gz") else file(filename,"w");
    # write stuff
    ms = get_ms();
    print "Opened MS %s"%msname;
    if colname not in ms.colnames():
      error("Column %s not found"%colname);
    nrows = ms.nrows();
    print "Exporting %s to %s:"%(colname,filename);
    cPickle.dump((options.row_step,nrows),gzf);
    for r0 in range(0,nrows,options.row_step):
      print "%d/%d"%(r0,nrows);
      col = ms.getcol(colname,r0,options.row_step);
      cPickle.dump(col,gzf);
    ms.close();
    gzf.close();
  except:
    traceback.print_exc();
    error("Error exporting %s to %s"%(colname,filename));
  print "Column exported OK.";
