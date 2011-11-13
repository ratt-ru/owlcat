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
import numpy
import Owlcat

ROWCHUNK = 50000;

time_start = time.time();

try:
  from pyrap_tables import table,tablecopy,tableexists,tabledelete
except:
  from pyrap.tables import table,tablecopy,tableexists,tabledelete

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
      description="""Finds redundant baselines and sets the WEIGHT (or
IMAGING_WEIGHT) column of the MS, weighting the redundant baselines as 1/n. For a redundant array such as WSRT, this considerably reduces sidelobes in the image.""");
  #parser.add_option("-o","--output",dest="output",type="string",
                    #help="name of output FITS file");
  #parser.add_option("-z","--zoom",dest="output_quad",type="string",
                    #help="name of zoomed output FITS file");
  parser.add_option("-t","--tolerance",dest="tolerance",type="float",
                    help="How close (in meters) two baselines need to be to each other to be considered redundant (default .1)");
  parser.add_option("-w","--weight",action="store_true",
                    help="If this option is given, the IMAGING_WEIGHT column will be "
                    "scaled by 1/n (where n is the number of redundant baselines). The "
                    "default behaviour is to simply set the WEIGHT column to 1/n.")
  parser.add_option("-I","--ifrs",dest="ifrs",type="string",
                    help="subset of interferometers to use.");
  parser.add_option("-s","--select",dest="select",action="store",
                    help="additional TaQL selection string. Note that redundant baselines are counted only within the subset "
                         "given by the --ifrs and --select options.");
  parser.add_option("-l","--list",action="store_true",
                    help="list all baselines and exit.");
  parser.add_option("-r","--reset",action="store_true",
                    help="reset WEIGHT column to unity. All other options are then ignored.");
  parser.set_defaults(tolerance=.1,select="",ifrs="");

  (options,msnames) = parser.parse_args();

  if len(msnames) != 1:
    parser.error("MS name not supplied.");

  msname = msnames[0];
  ms = table(msname,readonly=False);
  ncorr = ms.getcol('DATA',0,1).shape[2];

  if options.reset:
    # apply weights
    nrows = ms.nrows();
    for row0 in range(0,nrows,ROWCHUNK):
      progress("Resetting weights, row %d of %d"%(row0,nrows),newline=False);
      row1 = min(nrows,row0+ROWCHUNK);
      w = numpy.ones((row1-row0,ncorr),float);
      ms.putcol('WEIGHT',w,row0,row1-row0);
    progress("Weights reset to unity.");
  else:
    taqls = [];
    # get IFR set
    import Meow.IfrSet
    ifrset = Meow.IfrSet.from_ms(ms);
    if options.ifrs:
      ifrset = ifrset.subset(options.ifrs);
      taqls.append(ifrset.taql_string());

    if options.select:
      taqls.append(options.select);

    if taqls:
      select = "( " + " ) &&( ".join(taqls) + " )";
      progress("Applying TaQL selection %s"%select,newline=True);
      ms = ms.query(select);
    progress("Looking for redundant baselines",newline=True);
    ant1 = ms.getcol('ANTENNA1');
    ant2 = ms.getcol('ANTENNA2');

    IFRS = sorted(set([ (p,q) for p,q in zip(ant1,ant2) ]));
    print "%d baselines"%len(IFRS);
    groups = [];

    for i,(p,q) in enumerate(IFRS):
      bl = ifrset.baseline_vector(p,q);
      # see if this baseline is within the tolerance of a previous group's baseline
      for ig,(bl0,members) in enumerate(groups):
        if abs(bl-bl0).max() < options.tolerance:
          members.append((p,q));
          break;
      # if none, start a new group
      else:
        members = [(p,q)];
        ig = len(groups);
        groups.append([]);
      # update baseline (as mean baseline of group)
      length = reduce(lambda x,y:x+y,[ifrset.baseline_vector(*ifr) for ifr in members])/len(members);
      groups[ig] = (length,members);

    # convert to list of length,members, and sort by length
    groups = [ (math.sqrt((baseline**2).sum()),members) for baseline,members in groups ];
    groups.sort();

    if options.list:
      baselist = [ "%dm (%s)"%(round(length)," ".join(["%s-%s"%(p,q) for p,q in mem]))
          for length,mem in groups ];
      print "Found %d non-redundant baselines:"%len(baselist),", ".join(baselist);
      sys.exit(0);

    # make a dictionary of per-IFR weights
    have_redundancy = False;
    weight = dict([((p,q),1.) for p,q in IFRS]);
    for baseline,members in groups:
      if len(members) > 1:
        print "Baseline %dm, %d ifrs: %s"%(round(baseline),len(members),
          " ".join(["%d-%d"%(p,q) for p,q in members]));
        have_redundancy = True;
        for p,q in members:
          weight[p,q] = 1.0/len(members);

    if not have_redundancy:
      print "No redundant baselines found, nothing to do."
      sys.exit(0);

    # apply weights
    nrows = ms.nrows();
    for row0 in range(0,nrows,ROWCHUNK):
      progress("Applying weights, row %d of %d"%(row0,nrows),newline=False);
      row1 = min(nrows,row0+ROWCHUNK);
      if not options.weight:
        w = numpy.zeros((row1-row0,ncorr),float);
        for i in range(row0,row1):
          w[i-row0,:] = weight[ant1[i],ant2[i]];
        ms.putcol('WEIGHT',w,row0,row1-row0);
      else:
        w = ms.getcol('IMAGING_WEIGHT',row0,row1-row0);
        for i in range(row0,row1):
          w[i-row0,:] *= weight[ant1[i],ant2[i]];
        ms.putcol('IMAGING_WEIGHT',w,row0,row1-row0);

  progress("Closing MS.");
  ms.close();
