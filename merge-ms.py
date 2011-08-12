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

time_start = time.time();

from Owlcat import table,tablecopy,tableexists,tabledelete
import Owlcat.Console;

progress = Owlcat.Console.Reporter(timestamp=True);

SKIP_COLUMNS = set(("FLAG_CATEGORY","NFRA_AVERAGEDATA","PEELED_DATA","UVMODEL_DATA"));

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser
  parser = OptionParser(usage="""%prog: [options] MS_output MS_in1 MS_in2 ...""",
    description="Concatenates several MSs into one. This is done on a row-by-row basis "
    "with no sanity checking, so it's up to you to insure that the input MSs actually "
    "belong together. All subtables and keywords of the output MS are inherited from "
    "the first input MS, unless the --renumber-spws option is in effect, in which case "
    "the SPECTRAL_WINDOW and DATA_DESCRIPTION tables are also merged (and spwids and "
    "ddids renumbered.)");
  parser.add_option("-f","--force",dest="force",action="store_true",
                    help="proceed without confirmation, and overwrite output MS if it already exists");
  parser.add_option("-s","--renumber-spws",dest="renumber",action="store_true",
                    help="treat each MS as a separate spectral window");

  (options,msnames) = parser.parse_args();

  if len(msnames) < 3:
    parser.error("Insufficient number of arguments. Use '-h' for help.");

  msout = msnames[0];
  msins = msnames[1:];

  if tableexists(msout):
    if not options.force:
      print "Output MS",msout,"already exists. Use the -f switch to overwrite.";
      sys.exit(1);

  if not options.force:
    print "Will create merged MS %s from %d input MSs:"%(msout,len(msins));
    for ms in msins:
      print "  ",ms;
    if options.renumber:
      print "Each input MS will be put into a separate spectral window, spws will be renumbered.";
    else:
      print "Spectral windows will not be renumbered.";
    if raw_input("Proceed (y/n)? ").strip().upper()[0] != "Y":
      print "Aborted by user.";
      sys.exit(1);

  if tableexists(msout) and options.force:
    progress("Deleting existing copy of %s"%msout);
    tabledelete(msout);

  # copy first MS to output as-is
  progress("Copying %s to %s"%(msins[0],msout));
  tablecopy(msins[0],msout,deep=True);

  # open output for writing
  tab0 = table(msout,readonly=False);

  # get the number of DDIDs and SPWIDs
  if options.renumber:
    ddid_tab0 = table(tab0.getkeyword("DATA_DESCRIPTION"),readonly=False);
    spw_tab0 = table(tab0.getkeyword("SPECTRAL_WINDOW"),readonly=False);

  for msname in msins[1:]:
    tab = table(msname);
    nr0 = tab0.nrows();
    progress("Current row count: %d"%nr0);
    progress("Merging in %d rows from %s"%(tab.nrows(),msname));
    tab0.addrows(tab.nrows());
    for col in tab.colnames():
      if col not in SKIP_COLUMNS:
        try:
          data = tab.getcol(col);
        except:
          print "WARNING: MS %s does not appear to contain a valid %s column"%(msname,col);
          print "   This is not necessarily fatal, proceeding with the merge anyway."
          continue;
        # if renumbering DDIDs, increment the DDID column by the # of rows in the DDID table --
        # this ensures uniqueness of DDIDs
        if options.renumber and col == "DATA_DESC_ID":
          data += ddid_tab0.nrows();
        progress.overprint("Writing column %s, shape %s"%(col,data.shape));
        tab0.putcol(col,data,nr0);
    # if renumbering, need to concatenate the DDID and SPW tables
    if options.renumber:
      progress.overprint("Updating DATA_DESCRIPTION subtable");
      # append content of the DATA_DESCRIPTION table, while renumbering spectral window IDs
      ddid_tab = table(tab.getkeyword("DATA_DESCRIPTION"),readonly=False);
      nr0 = ddid_tab0.nrows();
      ddid_tab0.addrows(ddid_tab.nrows());
      for col in ddid_tab.colnames():
        data = ddid_tab.getcol(col);
        if col == "SPECTRAL_WINDOW_ID":
          data += spw_tab0.nrows();
        progress.overprint("Writing column %s, shape %s"%(col,data.shape));
        ddid_tab0.putcol(col,data,nr0);
      # append content of the SPECTRAL_WINDOW
      progress.overprint("Updating SPECTRAL_WINDOW subtable");
      spw_tab = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False);
      nr0 = spw_tab0.nrows();
      spw_tab0.addrows(spw_tab.nrows());
      for col in spw_tab.colnames():
        data = spw_tab.getcol(col);
        spw_tab0.putcol(col,data,nr0);


  progress.overprint("Closing output MS %s\n"%msout);
  nr0 = tab0.nrows();
  tab0.close();
  progress("Wrote %d rows to output MS"%nr0);
  if options.renumber:
    progress("Output MS contains %d spectral windows and %d DDIDs"%
      (spw_tab0.nrows(),ddid_tab0.nrows()));
    ddid_tab0.close();
    spw_tab0.close();
