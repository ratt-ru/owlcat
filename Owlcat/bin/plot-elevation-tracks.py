#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

if __name__ == "__main__":

    # setup some standard command-line option parsing
    #
    from optparse import OptionParser, OptionGroup
    from Owlcat.Plotting import LSTElevationPlot
    import casacore.tables

    parser = OptionParser(usage="""%prog: [plots & options] MS""",
                          description="""Plots elevation tracks from an MS.""")

    parser.add_option("-l", "--list", action="store_true",
                      help="lists fields found in MS, then exits")

    plotgroup = OptionGroup(parser, "Plotting options")
    outputgroup = OptionGroup(parser, "Output options")
    LSTElevationPlot.init_options(plotgroup, outputgroup)

    outputgroup.add_option("-o", "--output-name", type=str, help="Output filename", default="lst-elev.png")
    outputgroup.add_option("-d", "--display", action="store_true", help="Display plot on screen")

    parser.add_option_group(plotgroup)
    parser.add_option_group(outputgroup)

    (options, args) = parser.parse_args()

    if not args:
        parser.error("No MS specified")

    if options.display:
        options.output_name = None

    skyplot = LSTElevationPlot(options, output_type='x11' if options.display else 'png')

    for msname in args:
        field_time, field_radec, obs_xyz, scans = LSTElevationPlot.load_ms_fields(msname, verbose=1)
        if options.list:
            continue

        skyplot.make_figure(field_time, field_radec, obs_xyz, scans,
                            suptitle=msname, save=options.output_name, display=options.display)

    if options.display and not options.list:
        from pylab import plt
        plt.show()
