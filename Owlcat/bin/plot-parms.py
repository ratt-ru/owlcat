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
import sys

import numpy
import numpy.ma
import Owlcat.Console

from Owlcat.ParmTables import verbosity
from functools import reduce

progress = Owlcat.Console.Reporter(timestamp=True)


class PlottableEntityProps(object):
    """A PlottableEntityProps objects describes the properties of a plottable entity.
    Instances of this will be e.g. RealValue, ComplexValue, JonesMatrix, etc."""

    def __init__(self, name, plottables):
        """Creates a plottable entity description. 'name' is a short name, e.g. "real values".
        'plottables' is a list of rules for entities into real numbers than can be plotted, i.e.
        a list of ('rulename',func) tuples, where func transforms entities into float arrays."""
        self.name = name
        self.plottables = plottables
        self.plottables_dict = dict(plottables)


RealValue = PlottableEntityProps("real values", [("", lambda x: x)])

ComplexValue = PlottableEntityProps("complex values", [
    ('ampl', numpy.ma.absolute),
    ('phase', lambda x: numpy.ma.masked_array(numpy.angle(x), x.mask, fill_value=x.fill_value)),
    ('real', lambda x: x.real),
    ('imag', lambda x: x.imag),
])

JonesMatrix = PlottableEntityProps("Jones matrices", [
    ('norm', lambda matrix: numpy.sqrt(sum([x * numpy.ma.conjugate(x) for x in matrix])))
])


def _matrix22(a11, *args):
    """helper function for aggregators below: makes 2x2 matrix. When called with 1 args, makes
    scalar-diagonal matrix, with 2 makes diagonal matrix, with 4 makes full matrix."""
    if not args:
        return (a11, 0, 0, a11)
    elif len(args) == 1:
        return (a11, 0, 0, args[0])
    elif len(args) == 3:
        return (a11, args[0], args[1], args[2])
    raise TypeError("incorrect number of arguments in call to _matrix22")


# 'Aggregators' is a list of rules for aggregating individual real parm coefficients into compound
# entities. E.g. ':r' and ':i' parms can be aggregated to form a complex number, then xx,xy,yx,yy
# complex values can be aggregated into a 2x2 matrix, etc.
# Each entry in the list has the form:
#   suffix_list,func,prop
# where 'suffix_list' gives the suffixes of the elements ('r','i'), 'func' is a callable for transforming
# the elements into an aggregate value, and 'prop' is a PlottableEntityProps object describing the aggregate.
#
Aggregators = [
    (('r', 'i'), (lambda r, i: r + 1j * i), ComplexValue),
    (('r',), (lambda r: r), RealValue),
    (('ampl', 'phase'), (lambda ampl, phase: ampl * numpy.exp(1j * phase)), ComplexValue),
    (('ampl',), (lambda r: r), RealValue),
    (('xx', 'xy', 'yx', 'yy'), _matrix22, JonesMatrix),
    (('XX', 'XY', 'YX', 'YY'), _matrix22, JonesMatrix),
    (('rr', 'rl', 'lr', 'll'), _matrix22, JonesMatrix),
    (('RR', 'RL', 'LR', 'LL'), _matrix22, JonesMatrix),
    (('xx', 'yy'), _matrix22, JonesMatrix),
    (('rr', 'll'), _matrix22, JonesMatrix),
    (('XX', 'YY'), _matrix22, JonesMatrix),
    (('RR', 'LL'), _matrix22, JonesMatrix)
]


def _make_nested_func(agg_func, comp_func, n):
    """Helper method to make a nested aggregator function. Returns callable which takes any number of args a[],
    and returns agg_func(comp_func(a[0],...,a[n-1]),comp_func(a[n],...,a[2n-1]),...).
    If comp_func is None, returns simply arg_func."""
    if comp_func:
        return lambda *a: agg_func(*[comp_func(*a[i:i + n]) for i in range(0, len(a), n)])
    else:
        return agg_func


def makeParmEntityList(parmnames):
    """Processes list of parmnames and makes a list of plottable entities.
    Returns list of tuples:
      (entity_dict,ncomp,func,props), where:
    'desc' is a short description of this group of entities
    'entity_dict' is a dict of {name:components},
      'name' is an entity name (original parmname, or shortened aggregated name)
      'components' is a flat list of parmnames making up that aggregate
    'ncomp' is the number of components (each component list in the dict must have that many elements)
    'func' is a callable for transforming components into an aggregate value (None for individual parms)
    'props' is a PlottableEntityProps object.

    Each entry in the list corresponds to a group of entities that matched a particular aggregation rule.
    """
    # start with one dict: that of all parms
    entity_list = []
    new_entity_list = [(dict([(name, [name]) for name in parmnames]), 1, lambda x: x, RealValue)]
    aggregated = set()
    # entities are moved from new_entity_list to entity_list as they're processed.
    while new_entity_list:
        entity_dict, ncomp, func, props = new_entity_list[0]
        # see if we have recipe for making compound objects based on last parmname qualifier
        nqlist = [name.rsplit(':', 1) for name in list(entity_dict.keys())]
        # set of basenames and set of all occurring qualifiers
        basenames = set([nq[0] for nq in nqlist])
        quals = set([nq[1] for nq in nqlist])
        for quallist, agg_func, agg_props in Aggregators:
            # look for an aggregator whose qualifiers are a subset of what we have
            if set(quallist) <= quals:
                new_entity_dict = {}
                # for each basename, attempt to find entities for every qualifier in the list
                for basename in basenames:
                    fqnames = [basename + ':' + qq for qq in quallist]
                    components = [name not in aggregated and entity_dict.get(name, None) for name in fqnames]
                    # If every qualified name gives a valid component list, make an aggregated entity
                    # Components have to be transformed into a flat list.
                    if all(components):
                        new_entity_dict[basename] = reduce(lambda x, y: x + y, components)
                        aggregated.update(fqnames)
                # if any new entities were aggregated, add to new list for firther reprocessing
                if new_entity_dict:
                    new_entity_list.append((new_entity_dict,
                                            ncomp * len(quallist),
                                            _make_nested_func(agg_func, func, ncomp),
                                            agg_props))
        # move current entity to processed list
        entity_list.append(new_entity_list.pop(0))
    # flip list around, so latest match (presumably most complex entity) comes first
    entity_list.reverse()
    return entity_list


def namesToSetNotation(namelist):
    # build set of possible values for each field
    fsets = {}
    for name in namelist:
        for i, f in enumerate(name.split(':')):
            fsets.setdefault(i, set()).add(f)
    # return compound string
    return ":".join([(fsets[i].pop() if len(fsets[i]) == 1 else "{%s}" % ",".join(sorted(fsets[i])))
                     for i in range(len(fsets))])


if __name__ == "__main__":
    # setup some standard command-line option parsing
    #
    from optparse import OptionParser, OptionGroup

    parser = OptionParser(usage="""%prog: [options] <meptable> <plots...>""",
                          description="Makes plots of parmtables. Each plot is specified as '[+]name/what'. "
                                      "E.g. 'G:*:xx/phase G:*:xx/ampl' plots all G XX phases in one plot, and amplitudes "
                                      "in a second plot (assuming your MEP table actually contains G XX solutions.) By contrast, "
                                      "'G:*:xx/phase +G:*:xx/ampl' will put phases and amplitudes into one plot, with amplitudes "
                                      "using the right Y-axis (NB: the \"+\" feature currently doen't work). "
                                      "A second form produces scatterplots: 'name/whatX-whatY'. For example, G:*:xx/real-imag "
                                      "will produce a complex scatterplot. "
                                      "Names may contain shell-style wildcards and {a,b,c} constructs; for more flexibility, "
                                      "use the -x option to interpret them as regexes instead. "
                                      "Use the -l (--list) option to see what can "
                                      "actually be plotted from your table. If you don't specify anything, the default is to plot "
                                      "whatever is the first group printed by -l."
                          )
    parser.add_option("-x", "--regex", dest="regex", action="store_true",
                      help="parse plot specifications as regexes instead of shell-style wildcards."
                      )
    parser.add_option("-l", "--list", dest="list", action="store_true",
                      help="queries parmtable and lists what can be plotted")
    #  parser.add_option("-a","--axis-select",dest="axis_select",metavar="axis=start[:end[:[+]step]]",type="string",action="append",
    #                    help="Selects subset of data along an axis. E.g. 'time=0' selects the first timeslot, "
    #                    "'freq=0:99' selects the first 100 frequencies, 'freq:0:99:2' selects every second frequency "
    #                    "from the first 100, "
    #                    "'time=::100' selects every 100th timeslot, 'time=::+100' averages every block of 100 "
    #                    "timeslots. 'axis=all' or 'axis=::' selects everything. 'axis=avg' or 'axis=::+' selects and "
    #                    "averages everything. The -a option may be given multiple times for different axes, or "
    #                    "you can given one comma-separated list. The default for an axis that is not explicitly "
    #                    "selected is 'all' for time and 'avg' for the rest.")
    #  parser.add_option("-x","--xaxis",dest="xaxis",metavar="axis",type="string",
    #                  help="Selects which data axis to use for the plotted X axis. Default is 'time'.")
    parser.add_option("--average-before", dest="average_before", action="store_true",
                      help="apply averaging to parameters BEFORE conversion to plottable values.")
    parser.add_option("--average-after", dest="average_after", action="store_true",
                      help="apply averaging to plottable values, AFTER conversion from parameters (default).")
    parser.add_option("--offset", type="float",
                      help="use a fixed offset between plot tracks. Default is automatically chosen.")
    parser.add_option("--min", type="float", default=None,
                      help="filter out values <MIN")
    parser.add_option("--max", type="float", default=None,
                      help="filter out values >MAX")

    group = OptionGroup(parser, "Plot label options")
    group.add_option("--label-mean", dest="label_mean", action="store_true",
                     help="include mean value in plot labels")
    group.add_option("--no-label-mean", dest="label_mean", action="store_false",
                     help="do not include mean value in plot labels")
    group.add_option("--label-stddev", dest="label_stddev", action="store_true",
                     help="include standard deviation in plot labels")
    group.add_option("--no-label-stddev", dest="label_stddev", action="store_false",
                     help="do not include standard deviation in plot labels")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output options")
    group.add_option("-o", "--output", dest="output", type="string",
                     help="name of output file (e.g. 'foo.ps' or 'bar.png'). If not specified, plot will be displayed in a window.")
    group.add_option("--dpi", dest="resolution", type="int", metavar="DPI",
                     help="plot resolution for .ps/.png output (default is %default)")
    group.add_option("--plots-per-page", dest="ppp", metavar="N", type="int",
                     help="max plots per page. Default is 120.")
    group.add_option("--title", dest="title", type="string",
                     help="plot title (default is '<parmtab> <desc>')")
    group.add_option("-v", "--verbose", dest="verbose", type="int",
                     help="set verbosity level for debugging messages")
    parser.add_option_group(group)

    parser.set_defaults(output="", xaxis="time", resolution=300, ppp=120,
                        label_mean=True, label_stddev=True, offset=None,
                        )

    (options, args) = parser.parse_args()
    if not args:
        parser.error("MEP table not specified. Use '-h' for help.")
    tabname = args[0]
    plotspecs = args[1:]

    verbosity.set_verbose(options.verbose)

    # read table, figure out what's in it
    progress("Reading MEP table %s, please wait..." % tabname)
    from Owlcat import ParmTables
    from Timba import mequtils

    pt = ParmTables.open(tabname)

    name_list = pt.parmtable().name_list()
    entity_list = makeParmEntityList(name_list)

    if options.list:
        print(("\nTable contains %d funklets. We can plot the following things:\n" % len(name_list)))
        for i, (entity_dict, ncomp, func, props) in enumerate(entity_list):
            names = namesToSetNotation(iter(list(entity_dict.keys())))
            plts = ",".join([p[0] for p in props.plottables])
            plts = plts and ("/{%s}" % plts if len(props.plottables) > 1 else "/%s" % plts)
            print(("  %s: %s%s%s" % (props.name, names, plts, ("" if i else " (plotted by default)"))))
        print("\nValid axes are:\n")
        cells = pt.envelope_cells()
        for iaxis in range(mequtils.max_axis):
            axis = mequtils.get_axis_id(iaxis)
            grid = cells.grid.get(axis)
            if grid is not None:
                print(("  %s: %d points from %g to %g" % (axis, len(grid), grid[0], grid[-1])))
        print("")
        sys.exit(0)

    #
    # parse plot specifications
    #
    import fnmatch
    import re

    # This will contain a definitive list of all plots. Each plot is specified as a
    # [description,entity_list,entity_list2,scatterplot]. entity_list2 is valid for dual plots,
    # and is None for a single plot.
    # Scatterplot is None for normal plots, and a 2-tuple of axis labels, for scatter plots.
    Plots = []

    # go through plotlist and find matching entities. regexes must match field-by-field
    if plotspecs:
        for plotspec in plotspecs:
            # '+' sign specifies second plot on previous page
            secondplot = (plotspec[0] == '+')
            if secondplot:
                if not Plots:
                    parser.error("Can't use '+' in '%s': no preceding plot." % plotspec)
                if Plots[-1][2] is not None:
                    parser.error("Can't use '+' in '%s': preceding plot is already dual." % plotspec)
                if Plots[-1][3]:
                    parser.error("Can't use '+' in '%s': preceding plot is a scatterplot." % plotspec)
                plotspec = plotspec[1:]
            # split into name/what components
            names, what = plotspec.rsplit("/", 1) if "/" in plotspec else (plotspec, None)
            # what may be an "X-Y" combo (for scatterplots)
            # in this case scatterplot is set to "X","Y", else to None.
            scatterplot = what.split("-", 1) if what else []
            if len(scatterplot) != 2:
                scatterplot = None
            else:
                if secondplot:
                    parser.error("Can't use '+' in '%s': scatterplots cannot be dual." % plotspec)
            # split name into fields at each ":", and translate into regexes
            if options.regex:
                regexes = list(map(re.compile, names.split(':')))
            else:
                # glob mode: translate the {a,b,c} case to a regex, and use
                # fnmatch.translate() for the rest      else
                regexes = [
                    re.compile("(%s)" % ('|'.join(ss[1:-1].split(","))) if (ss and ss[0] == '{' and ss[-1] == "}")
                               else fnmatch.translate(ss)) for ss in names.split(":")]
            # build list of what to plot, as (name,components,func,plotfunc,plotfunc2) tuples
            # plotfunc2 is None for nroaml plots, or a plotfunc for scatterplots
            plot_entlist = []
            for entity_dict, ncomp, func, props in entity_list:
                for name, comps in sorted(entity_dict.items()):
                    # name must have same number of fields as the original regex, and all must match
                    nf = name.split(':')
                    if len(nf) == len(regexes) and all([r.match(f) for r, f in zip(regexes, nf)]):
                        # check that 'what' is a plottable
                        if what:
                            if scatterplot:
                                plotfuncs = [props.plottables_dict.get(wh) for wh in scatterplot]
                                plotfunc, plotfunc2 = plotfuncs
                            else:
                                plotfuncs = [props.plottables_dict.get(what)]
                                plotfunc, plotfunc2 = plotfuncs[0], None
                            if not all(plotfuncs):
                                print(("%s/%s cannot be plotted. Use the -l option to see what can." % (name, what)))
                                sys.exit(2)
                        # if 'what' not given, use first item in list
                        else:
                            what, plotfunc = props.plottables[0]
                            plotfunc2 = None
                        # add to list
                        # print "Plotting %s %s '%s' %s %s"%(name,comps,what,func,plotfunc)
                        plot_entlist.append((name, comps, func, plotfunc, plotfunc2))
            # have we found something to plot?
            if plot_entlist:
                if secondplot:
                    Plots[-1][2] = plot_entlist
                else:
                    Plots.append([plotspec, plot_entlist, None, scatterplot])
                progress("Plot %d: %s%s%s" % (len(Plots), namesToSetNotation([e[0] for e in plot_entlist]),
                                              "/" + what if what else "", " (second axis)" if secondplot else ""))
            else:
                print(("Nothing found to match '%s', ignoring." % plotspec))
    # if no plot options supplied, default is to take first group of plottables.
    else:
        entity_dict, ncomp, func, props = entity_list[0]
        what, plotfunc = props.plottables[0]
        progress("Plot 1: %s%s" % (namesToSetNotation(iter(list(entity_dict.keys()))), "/" + what if what else ""))
        Plots.append(
            ["", [(name, comps, func, plotfunc, None) for name, comps in sorted(entity_dict.items())], None, False])

    if not Plots:
        print("Found nothing to plot. Use the -l option to see what may be plotted.")
        sys.exit(0)

    # set non-interactive backend, if -o is in effect
    import matplotlib

    if options.output:
        matplotlib.use('agg')
    else:
        matplotlib.use('qt4agg')

    # ok, make the plots
    from Owlcat.Plotting import PlotCollection, ScatterPlot

    # make a page counter, if multiple pages are expected
    if len(Plots) > 1 or any(
            [len(entlist) > options.ppp for (desc, entlist, entlist2, scatterplot) in Plots if not scatterplot]):
        ipage = 0
    else:
        ipage = None

    for iplot, (desc, entlist, entlist2, scatterplot) in enumerate(Plots):
        # setup title
        title = options.title or "%s: %s" % (tabname, desc)

        # make X-Y scatterplot
        if scatterplot:
            plot = ScatterPlot(options)
            plot.desc = desc
            progress("Starting plot %d of %d, %d data points" % (iplot + 1, len(Plots), len(entlist)))
            # loop over all plottables and compute data for each
            for name, comps, func, plotfunc, plotfunc2 in entlist:
                progress.overprint("Getting data for %s" % name)
                # get per-component arrays from ParmTable
                arrays = [pt.funkset(comp).array() for comp in comps]
                if options.average_before:
                    arrays = list(map(numpy.ma.mean, arrays))
                dd = func(*arrays)
                x, y = plotfunc(dd), plotfunc2(dd)
                if options.average_after:
                    x, y = list(map(numpy.ma.mean, (x, y)))
                # add point
                plot.add_point(name, x, y)

            # make filename, if saving file
            savefile = options.output
            if ipage is not None:
                if savefile:
                    basename, ext = os.path.splitext(savefile)
                    savefile = "%s.%d%s" % (basename, ipage, ext)
                    progress("Plot will be written to file %s\n" % savefile)
                ipage += 1
            progress("Generating plot")
            plot.make_figure(save=savefile, suptitle=title, dpi=options.resolution, xlabel=scatterplot[0],
                             ylabel=scatterplot[1])

        # normal "collections" plot
        else:
            # figure out label format
            labels = ["%(name)s"]
            if options.label_mean:
                labels.append("mean=%(mean).3g")
            if options.label_stddev:
                labels.append("std=%(stddev).3g")
            label_format = " ".join(labels)
            # make collection object
            coll = PlotCollection(options)
            if options.offset:
                coll.offset = options.offset
            coll.desc = desc
            progress("Starting plot %d of %d, %d data tracks" % (iplot + 1, len(Plots), len(entlist)))
            # loop over all plottables and compute data for each
            for name, comps, func, plotfunc, plotfunc2 in entlist:
                progress.overprint("Getting data for %s" % name)
                # get per-component arrays from ParmTable
                arrays = [pt.funkset(comp).array() for comp in comps]
                if options.average_before:
                    arrays = [arr.mean(1) for arr in arrays]
                # apply functions to get plottable values
                dd = plotfunc(func(*arrays))
                # apply after-averaging, if needed
                if not options.average_before:
                    dd = dd.mean(1)
                mask = None
                if options.min is not None:
                    mask = dd < options.min
                if options.max is not None:
                    mask = mask & (dd > options.max) if mask is not None else (dd > options.max)
                if mask is not None:
                    dd = numpy.ma.masked_array(dd, mask)
                # put together label
                labelattrs = dict(name=name, mean=dd.mean(), stddev=dd.std())
                # add track
                coll.add_track(name, dd, label=label_format % labelattrs)

            for np0 in range(0, len(entlist), options.ppp):
                np1 = min(np0 + options.ppp, len(entlist))
                names = [ent[0] for ent in entlist[np0:np1]]
                # make filename, if saving file
                savefile = options.output
                if ipage is not None:
                    if savefile:
                        basename, ext = os.path.splitext(savefile)
                        savefile = "%s.%d%s" % (basename, ipage, ext)
                        progress("Plot will be written to file %s\n" % savefile)
                    ipage += 1
                # dual plot?
                dual = (len(names) > options.ppp / 2)
                progress("Generating plot")
                coll.make_figure(names, save=savefile, suptitle=title, dpi=options.resolution, dual=dual)

    if not options.output:
        from pylab import plt

        plt.show()
