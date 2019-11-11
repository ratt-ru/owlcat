#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# % $Id$
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

from Owlcat import table, tableexists, tabledelete

if __name__ == "__main__":

    # setup some standard command-line option parsing
    #
    from optparse import OptionParser

    parser = OptionParser(usage="""%prog: [options] MS""",
                          description="Splits MS into one or more output MSs by DATA_DESC_ID (which typically corresponds "
                                      "to spectral window.)")
    parser.add_option("-D", "--ddid", type="int",
                      help="specific DATA_DESC_ID to extract. If not given, then all DDIDs are extracted, each into a separate MS.")
    parser.add_option("-o", "--output", type="string",
                      help="Name of output MS. Use '%(ms)s' to insert basename of input MS. Use '%(ddid)d' to insert DDID number. Use '%(msext)s' to insert extension of input MS. Default is %default")
    parser.add_option("-f", "--force", dest="force", action="store_true",
                      help="overwrites output MS if it already exists")

    parser.set_defaults(ddid=-1, output="%(ms)s_spw%(ddid)d%(msext)s")

    (options, msname) = parser.parse_args()

    if len(msname) != 1:
        parser.error("Incorrect number of arguments. Use '-h' for help.")

    # open input MS
    msname = msname[0]
    print(("Reading input MS %s" % msname))
    if not tableexists(msname):
        parser.error("MS %s not found." % msname)
    ms = table(msname)
    # check DDID option
    num_ddids = table(ms.getkeyword("DATA_DESCRIPTION")).nrows()
    print(("MS contains %d DATA_DESC_IDs" % num_ddids))
    if options.ddid < 0:
        ddids = list(range(num_ddids))
    elif options.ddid >= num_ddids:
        parser.error("DATA_DESC_ID %d is out of range" % options.ddid)
    else:
        ddids = [options.ddid]

    # setup outputs
    msname, msext = os.path.splitext(os.path.basename(os.path.normpath(msname)))

    for ddid in ddids:
        msout = options.output % dict(ms=msname, msext=msext, ddid=ddid)
        print(("Extracting DATA_DESC_ID %d into MS %s" % (ddid, msout)))
        # delete if exists
        if tableexists(msout):
            if not options.force:
                parser.error("Output MS %s already exists. Use the -f switch to overwrite." % msout)
            print(("Deleting existing copy of %s" % msout))
            tabledelete(msout)
        # extract
        ms1 = ms.query("DATA_DESC_ID==%d" % ddid)
        ms1.copy(msout, deep=True)
        ms1.close()
