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

from Owlcat import table,tablecopy,tableexists,tabledelete,addImagingColumns
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
  parser.add_option("-a","--append",action="store_true",
                    help="append input MSs to output MS");
  parser.add_option("-f","--force",dest="force",action="store_true",
                    help="proceed without confirmation, and overwrite output MS if it already exists");
  parser.add_option("-s","--renumber-spws",action="store_true",
                    help="treat each MS as a separate spectral window");
  parser.add_option("--renumber-obs",action="store_true",
                    help="renumber OBSERVATION_ID of each MS");
  parser.add_option("-v","--verbose",action="store_true",
                    help="be more verbose");

  (options,msnames) = parser.parse_args();

  if len(msnames) < (2 if options.append else 3):
    parser.error("Insufficient number of arguments. Use '-h' for help.");

  msout = msnames[0];
  msins = msnames[1:];

  if tableexists(msout) and not options.append and not options.force:
    print "Output MS",msout,"already exists. Use the -f switch to overwrite.";
    sys.exit(1);

  if not options.force:
    if options.append:
      print "Will extend merged MS %s with %d input MSs:"%(msout,len(msins));
    else:
      print "Will create merged MS %s from %d input MSs:"%(msout,len(msins));
    for ms in msins:
      print "  ",ms;
    if options.renumber_spws:
      print "Each input MS will be put into a separate spectral window, spws will be renumbered.";
    else:
      print "Spectral windows will not be renumbered.";
    if options.renumber_obs:
      print "OBSERVATION_ID will be incremented for each MS";
    if raw_input("Proceed (y/n)? ").strip().upper()[0] != "Y":
      print "Aborted by user.";
      sys.exit(1);
      
  obsid = 0;

  if not options.append:
    if tableexists(msout) and options.force:
      progress("Deleting existing copy of %s"%msout);
      tabledelete(msout);
    # copy first MS to output as-is
    progress("Copying %s to %s"%(msins[0],msout));
    tablecopy(msins[0],msout,deep=True);
    # when renumbering spectral windows, need to have an addImagingColumns() call
    # unfortunately, this wipes CORRECTED_DATA! I see no way around this (copynorows=True is no help,
    # as then it fails to copy subtables) apart from the ugly: copy first MS deeply,
    # add imaging columns, then re-copy columns
    if options.renumber_spws:
      addImagingColumns(msout);
    # otherwise, first MS is already copied, so remove it from list
    # (but reset OBSERVATION_ID)
    else:
      if options.renumber_obs:
        tab = table(msout,readonly=False);
        obs = tab.getcol("OBSERVATION_ID");
        obs[:] = obsid;
        obsid += 1;
        tab.putcol("OBSERVATION_ID",obs);
        tab.close();
      msins = msins[1:];

  # open output for writing
  tab0 = table(msout,readonly=False);

  # get the number of DDIDs and SPWIDs
  ddid_offsets = [0]*len(msins);
  if options.renumber_spws:
    ddid_tab0 = table(tab0.getkeyword("DATA_DESCRIPTION"),readonly=False);
    spw_tab0 = table(tab0.getkeyword("SPECTRAL_WINDOW"),readonly=False);
    # now, somehow lwimager will still try to re-init MODEL_DATA and CORRECTED_DATA when dealing with merged
    # MSs. The only way around it that I can see is to fill up the ddid and spw tabs first
    progress("Filling new DATA_DESCRIPTION and SPECTRAL_WINDOW table");
    for num_ms,ms in enumerate(msins[1:]):
      tab = table(ms);
      ddid_tab = table(tab.getkeyword("DATA_DESCRIPTION"),readonly=False);
      nr0 = ddid_tab0.nrows();
      ddid_offsets[num_ms+1] = nr0;
      ddid_tab0.addrows(ddid_tab.nrows());
      for col in ddid_tab.colnames():
        data = ddid_tab.getcol(col);
        if col == "SPECTRAL_WINDOW_ID":
          data += spw_tab0.nrows();
        ddid_tab0.putcol(col,data,nr0);
      # append content of the SPECTRAL_WINDOW
      spw_tab = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False);
      nr0 = spw_tab0.nrows();
      spw_tab0.addrows(spw_tab.nrows());
      for col in spw_tab.colnames():
        data = spw_tab.getcol(col);
        spw_tab0.putcol(col,data,nr0);
    num_spws,num_ddids = spw_tab0.nrows(),ddid_tab0.nrows();
    ddid_tab0.close();
    spw_tab0.close();
    tab0.close();
    addImagingColumns(msout);
    tab0 = table(msout,readonly=False);

  overprint = progress if options.verbose else progress.overprint;

  for num_ms,msname in enumerate(msins):
    tab = table(msname);
    nr0 = tab0.nrows();
    # in spw renumbering mode, because of the ugliness explained above,
    # we need to re-copy the data from the first MS
    if options.renumber_spws and not num_ms:
      nr0 = 0;
    # else, insert more rows to accommodate added MS
    else:
      tab0.addrows(tab.nrows());
    progress("Current row count: %d"%nr0);
    progress("Merging in %d rows from %s"%(tab.nrows(),msname));
    for col in tab.colnames():
      if col not in SKIP_COLUMNS:
        overprint("Reading column %s"%col);
        try:
          data = tab.getcol(col);
        except:
          print "WARNING: MS %s does not appear to contain a valid %s column"%(msname,col);
          print "   This is not necessarily fatal, proceeding with the merge anyway."
          continue;
        if options.renumber_obs and col == "OBSERVATION_ID":
          data[:] = obsid;
          obsid += 1;
        # if renumbering DDIDs, increment the DDID column by the # of rows in the DDID table --
        # this ensures uniqueness of DDIDs
        if options.renumber_spws and col == "DATA_DESC_ID":
          progress("Adding %d to DATA_DESC_ID"%ddid_offsets[num_ms]);
          data += ddid_offsets[num_ms];
        overprint("Writing column %s, shape %s"%(col,data.shape));
        tab0.putcol(col,data,nr0);


  overprint("Closing output MS %s\n"%msout);
  nr0 = tab0.nrows();
  progress("Wrote %d rows to output MS"%nr0);
  if options.renumber_spws:
    progress("Output MS contains %d spectral windows and %d DDIDs"%(num_spws,num_ddids));
  tab0.close();
