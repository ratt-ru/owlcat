# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy
import numpy.ma
import math
import sys
import traceback
import os.path
import matplotlib
import matplotlib.dates
import casacore.measures

DEG = math.pi / 180
ARCMIN = math.pi / (180 * 60)


def init_pyplot(output_type="X11"):
    """Initializes the pyplot singleton. Inits the matplotlib back-end
    corresponding to the specified plot type"""
    if 'pyplot' not in globals():
        if output_type.upper() != "X11":
            matplotlib.use('agg')
        else:
            matplotlib.use('qt4agg')
        global pyplot
        import matplotlib.pyplot as pyplot


def _int_or_0(x):
    try:
        return int(x)
    except:
        return 0


_version = list(map(_int_or_0, matplotlib.__version__.split('.')))
_version_major = (_version and _version[0]) or 0
_version_minor = (len(_version) > 1 and _version[1]) or 0


class PlotCollection(object):
    """PlotCollection plots a collection of data tracks in one plot"""

    def __init__(self, options):
        self.data = {}
        self.mean = {}
        self.stddev = {}
        self.label = {}
        self.track_counts = {}
        # inter-plot offset -- may be set from outside
        self.offset = 0

    def num_tracks(self):
        return len(self.data)

    def add_track(self, key, data, mean=None, stddev=None, label=None, count=1):
        self.data[key] = data
        self.mean[key] = mean if mean is not None else data.mean()
        self.stddev[key] = stddev if stddev is not None else data.std()
        self.label[key] = label if label is not None else str(key)
        self.track_counts[key] = count

    def get_track_data(self, key):
        return self.data.get(key, None)

    def get_track_count(self, key):
        return self.track_counts.get(key, 0)

    # function to make single-page plot
    def make_figure(self, keylist=None, suptitle=None, save=None,
                    dual=True, offset_std=10, xgrid=False, ygrid=False,
                    figsize=(210, 290), dpi=100, papertype='a4', landscape=False):
        import matplotlib.pyplot as pyplot
        # setup sizes
        keylist = keylist or sorted(self.data.keys())
        # compute inter-plot offset, based on the median stddev
        offset = self.offset
        if not offset:
            stddevs = sorted(self.stddev.values())
            if stddevs:
                offset = self.offset = offset_std * stddevs[len(stddevs) // 2]
        # if still not set, force to 1
        if not offset:
            offset = self.offset = 1
        # partition keys into two sets, if asked to
        if dual:
            # split ifrs into two sets (if number is odd, make sure first subset has the extra 1)
            n = len(keylist)
            n2 = len(keylist) // 2
            if n % 2:
                n2 += 1
            keysets = [keylist[0:n2], keylist[n2:]]
        else:
            keysets = [keylist]
        # create Figure object of given size and resolution
        figsize_in = (figsize[0] / 25.4, figsize[1] / 25.4)
        fig = pyplot.figure(figsize=figsize_in, dpi=100)
        # margin sizes. We want to keep them fixed in absolute terms,
        # regardless of plot size. The relative numbers here are good for a 210x290 plot, so
        # we rescale them accordingly.
        mleft = 0.05 * 210. / figsize[0]
        mbottom = 0.03 * 290. / figsize[1]
        mright = 0.01 * 210. / figsize[0]
        mtop = 0.03 * 290. / figsize[1]
        ytitle = 1 - 0.01 * 290. / figsize[1]
        width = (1 - mleft - mright)
        height = (1 - mtop - mbottom)
        # now plot all tracks
        for iplot, keys in enumerate(keysets):
            if dual:
                plt = fig.add_subplot(1, 2, iplot + 1)
            else:
                plt = fig.add_axes([mleft, mbottom, width, height])
            nx = max([len(self.data.get(key, [])) for key in keys])
            # first and last key
            key0 = keys[0]
            key1 = keys[-1]
            mindata, maxdata = 1e+99, -1e+99
            # make plots for all ifrs
            firstkey = None
            for key in keys:
                data = self.data.get(key)
                if data is None or data.mask.all():
                    continue
                else:
                    # print data.min(),data.max()
                    # first plot goes at its natural coordinates.
                    if firstkey is None:
                        firstkey = key
                        y0 = 0
                        y0text = self.mean[key] - offset / 2
                        mindata = data.min() - offset / 2
                        #            mindata = y0text = self.mean[key] - offset/2
                        dd = data
                    # subsequent plots offset accordingly
                    else:
                        y0 += offset
                        dd = data + y0
                    #          maxdata = self.mean[key] + y0 + offset/2
                    mindata = min(mindata, data.min() + y0 - offset / 10)
                    maxdata = max(maxdata, self.mean[key] + y0 + offset)
                    # plot data
                    try:
                        if len(dd) == 1:
                            dd = numpy.array([dd[0], dd[0]])
                        plt.plot(dd)
                        if ygrid:
                            plt.axhline(y0, color='0.8', zorder=-10)
                        # print dd
                    except:
                        traceback.print_exc()
                        print("Error plotting data for", key)
                        continue
                    if self.label[key]:
                        # figure out where to place text label
                        ytext = dd[:nx // 10].mean()
                        if numpy.isnan(ytext):
                            ytext = dd.mean()
                        # do not put too close to previous label
                        y0text = max(y0text + offset / 4, ytext)
                        # annotate() seems to have a bug in matplotlib 0.98, so use text() instead
                        if _version_major > 0 or _version_minor > 98:
                            plt.annotate(self.label[key], xy=(nx / 100, ytext), xytext=(nx / 100, y0text),
                                         size=5, horizontalalignment='left', verticalalignment='center',
                                         bbox=dict(facecolor='white', edgecolor='none', alpha=0.6),
                                         arrowprops=dict(fc='none', width=0, headwidth=0, ec="0.8", alpha=0.6))
                        else:
                            plt.text(nx / 100, ytext, self.label[key],
                                     size=5, horizontalalignment='left', verticalalignment='center',
                                     bbox=dict(facecolor='white', edgecolor='none', alpha=0.6))

            # set plot limits. First panel is determined by data. Scale of second panel is fixed to first.
            plt.set_xbound(0, nx)
            if iplot == 0:
                ylim = mindata, maxdata
                plt.set_ybound(*ylim)
                for lab in plt.get_yticklabels():
                    lab.set_fontsize(5)
            else:
                for lab in plt.get_yticklabels():
                    lab.set_fontsize(0)
                # plt.set_yticklabels([])
                plt.set_ybound(*ylim)
            for lab in plt.get_xticklabels():
                lab.set_fontsize(5)
            # set x grid
            if _version_major > 0 or _version_minor > 98:
                plt.minorticks_on()
            if xgrid and xgrid != [0]:
                if isinstance(xgrid, (list, tuple)):
                    if len(xgrid) > 0:
                        plt.set_xticks(list(range(0, nx, xgrid[0])))
                    if len(xgrid) > 1:
                        plt.set_xticks(list(range(0, nx, xgrid[1])), True)
                for x in plt.get_xticks():
                    plt.axvline(x, color='0.7', ls='-', lw=1, zorder=-10)
                for x in plt.get_xticks(minor=True):
                    plt.axvline(x, color='0.9', ls='-', lw=1, zorder=-10)
                plt.set_xlim(0, nx)
        fig.subplots_adjust(left=mleft, right=1 - mright, top=1 - mtop, bottom=mbottom, wspace=0.01)
        # plot title if asked to
        if suptitle:
            fig.suptitle(suptitle, y=ytitle, size=8)
        if save:
            fig.savefig(save, papertype=papertype, dpi=dpi,
                        orientation='portrait' if not landscape else 'landscape')
            print("===> Wrote", save)
        return fig


class PlotCollectionSep(PlotCollection):
    """PlotCollectionSep is a PlotCollection that separates data tracks into individual plots,
    and optionally groups them by key"""

    def __init__(self, options):
        PlotCollection.__init__(self, options)
        self.groups = {}

    def num_tracks(self):
        return len(self.data)

    def add_track(self, key, data, mean=None, stddev=None, label=None, count=1, group=None):
        self.groups.setdefault(group, set()).add(key)
        PlotCollection.add_track(self, key, data, mean=mean, stddev=stddev, label=label, count=count)

    # function to make single-page plot
    def make_figure(self, grouplist=None, suptitle=None, save=None,
                    dual=True, xgrid=False, ygrid=False,
                    figsize=(210, 290), dpi=100, papertype='a4', landscape=False):
        import matplotlib.pyplot as pyplot
        # setup sizes
        grouplist = grouplist or sorted(self.groups.keys())
        # partition keys into two sets, if asked to
        if dual:
            # split ifrs into two sets (if number is odd, make sure first subset has the extra 1)
            n = len(grouplist)
            nrow = len(grouplist) / 2
            if n % 2:
                nrow += 1
            # swap things around, so plots start at the bottom
            plotnums = numpy.arange(1, nrow * 2 + 1).reshape(nrow, 2)[-1::-1].transpose().ravel()
        else:
            nrow = len(grouplist)
            plotnums = numpy.arange(nrow, 0, -1)
        # create Figure object of given size and resolution
        figsize_in = (figsize[0] / 25.4, figsize[1] / 25.4)
        fig = pyplot.figure(figsize=figsize_in, dpi=100)
        # margin sizes. We want to keep them fixed in absolute terms,
        # regardless of plot size. The relative numbers here are good for a 210x290 plot, so
        # we rescale them accordingly.
        mleft = 0.05 * 210. / figsize[0]
        mbottom = 0.03 * 290. / figsize[1]
        mright = 0.01 * 210. / figsize[0]
        mtop = 0.03 * 290. / figsize[1]
        ytitle = 1 - 0.01 * 290. / figsize[1]
        width = (1 - mleft - mright)
        height = (1 - mtop - mbottom)
        colors = ("blue", "green", "red", "cyan", "purple", "magenta", "black", "orange", "grey")
        # now plot all tracks
        for iplot, groupname in enumerate(grouplist):
            plt = fig.add_subplot(nrow, 2 if dual else 1, plotnums[iplot])
            keys = list(self.groups[groupname])
            nx = max([len(self.data.get(key, [])) for key in keys])
            # first and last key
            key0 = keys[0]
            key1 = keys[-1]
            mindata, maxdata = 1e+99, -1e+99
            # make plots for all ifrs
            for ikey, key in enumerate(keys):
                dd = self.data.get(key)
                if dd is None or dd.mask.all():
                    continue
                else:
                    try:
                        if len(dd) == 1:
                            dd = numpy.array([dd[0], dd[0]])
                        plt.plot(dd, ',', color=colors[ikey % len(colors)], label=self.label[key])
                    except:
                        traceback.print_exc()
                        print("Error plotting data for", key)
                        continue
            # add text labels
            # legend is not quite functional in matplotlib 0.98, so use text() instead
            y0, y1 = plt.get_ylim()
            if _version_major > 0 or _version_minor > 98:
                legend = plt.legend(ncol=len(keys), numpoints=1, prop=dict(size=5), handletextpad=0, columnspacing=0)
                #      legend.set_frame_on(False)
                for i, txt in enumerate(legend.get_texts()):
                    txt.set_color(colors[i % len(colors)])
                legend.get_frame().set_edgecolor('None')
                legend.get_frame().set_facecolor('white')
                legend.get_frame().set_alpha(.6)
            else:
                x0 = 0
                for ikey, key in enumerate(keys):
                    txt = plt.text(x0, y1, self.label[key], color=colors[ikey % len(colors)],
                                   size=5, horizontalalignment='left', verticalalignment='top',
                                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.6))
                    x0 += nx / 20
            plt.set_ylabel(groupname, size=8)
            # set x grid
            if _version_major > 0 or _version_minor > 98:
                plt.minorticks_on()
            for lab in plt.get_xmajorticklabels():
                lab.set_fontsize(5)
            for lab in plt.get_xminorticklabels():
                lab.set_fontsize(0)
            for lab in plt.get_ymajorticklabels():
                lab.set_fontsize(5)
            for lab in plt.get_yminorticklabels():
                lab.set_fontsize(0)
            if xgrid and xgrid != [0]:
                if isinstance(xgrid, (list, tuple)):
                    if len(xgrid) > 0:
                        plt.set_xticks(list(range(0, nx, xgrid[0])))
                    if len(xgrid) > 1:
                        plt.set_xticks(list(range(0, nx, xgrid[1])), True)
                for x in plt.get_xticks():
                    plt.axvline(x, color='0.7', ls='-', lw=1, zorder=-10)
                for x in plt.get_xticks(minor=True):
                    plt.axvline(x, color='0.9', ls='-', lw=1, zorder=-10)
                plt.set_xlim(0, nx)
            if ygrid:
                for y in plt.get_yticks():
                    plt.axhline(y, color='0.7', ls='-', lw=1, zorder=-10)
                for y in plt.get_yticks(minor=True):
                    plt.axhline(y, color='0.9', ls='-', lw=1, zorder=-10)
                plt.set_ylim(y0, y1)

        fig.subplots_adjust(left=mleft, right=1 - mright, top=1 - mtop, bottom=mbottom, wspace=0.12)
        # plot title if asked to
        if suptitle:
            fig.suptitle(suptitle, y=ytitle, size=8)
        if save:
            fig.savefig(save, papertype=papertype, dpi=dpi,
                        orientation='portrait' if not landscape else 'landscape')
            print("===> Wrote", save)
        return fig


class ComplexCirclePlot(PlotCollection):
    """ComplexCirclePlot plots a complex circle plot"""
    # plot colors: xx xy yx yy
    colors = ("blue", "red", "purple", "green")

    # function to make single-page plot
    def make_figure(self, keylist=None, suptitle=None, save=None, xgrid=False, ygrid=False,
                    dual=True, offset_std=10,
                    figsize=(210, 290), dpi=100, papertype='a4', landscape=False):
        import matplotlib.pyplot as pyplot
        figsize_in = (figsize[0] / 25.4, figsize[1] / 25.4)
        fig = pyplot.figure(figsize=figsize_in, dpi=100)
        if keylist:
            datalist = [(key, self.data.get(key)) for key in keylist]
        else:
            datalist = iter(list(data.items()))
        # margin sizes. We want to keep them fixed in absolute terms,
        # regardless of plot size. The relative numbers here are good for a 210x210 plot, so
        # we rescale them accordingly.
        mleft = 0.06 * 210. / figsize[0]
        mbottom = 0.06 * 210. / figsize[1]
        mright = 0.01 * 210. / figsize[0]
        mtop = 0.04 * 210. / figsize[1]
        ytitle = 1 - 0.01 * 210. / figsize[1]
        width = (1 - mleft - mright)
        height = (1 - mtop - mbottom)
        plt = fig.add_axes([mleft, mbottom, width, height])
        # find length of common prefix of labels
        labels = [lab for lab in list(self.label.values()) if lab]
        if labels:
            prefix = 1
            while len(labels[0]) > prefix and all([lab.startswith(labels[0][:prefix]) for lab in labels]):
                prefix += 1
        # plot data over text labels
        ncorr = 0
        rad = 0
        for key, data in datalist:
            label = self.label.get(key)
            ncorr = max(ncorr, data.shape[-1])
            for icorr in range(data.shape[-1]):
                color = self.colors[icorr]
                d = data[..., icorr]
                rad = max(rad, abs(d).max())
                if d.size > 1 or not label:
                    plt.plot(d.real, d.imag, ".-", mec=color, mfc=color, zorder=0)
                if label:
                    d0 = d[-1] if d.ndim else d
                    plt.text(d0.real, d0.imag, label[prefix - 1:], horizontalalignment='center',
                             verticalalignment='center', zorder=1, fontsize=8, color=color)
        # plot arrow to means
        for icorr in range(ncorr):
            color = self.colors[icorr]
            datameans = [d[..., icorr].mean() for d in list(self.data.values()) if d.shape[-1] > icorr]
            meanval = sum(datameans) / len(datameans) if datameans else 0
            plt.arrow(0, 0, meanval.real, meanval.imag, color=color)
            # plot circle
            phi = numpy.arange(0, 1.001, 0.001) * 2 * math.pi
            plt.plot(numpy.cos(phi) * abs(meanval), numpy.sin(phi) * abs(meanval), "-", color=color)
        for lab in plt.get_xticklabels():
            lab.set_fontsize(10)
        for lab in plt.get_yticklabels():
            lab.set_fontsize(10)
            lab.set_rotation('vertical')
        if xgrid or ygrid:
            plt.grid(color="0.4", linestyle='-')
        # plot title if asked to
        plt.set_xlim(-rad, rad)
        plt.set_ylim(-rad, rad)
        if suptitle:
            fig.suptitle(suptitle, y=ytitle, size=10)
        if save:
            fig.savefig(save, papertype=papertype, dpi=dpi,
                        orientation='portrait' if not landscape else 'landscape')
            print("===> Wrote", save)
        return fig


class ScatterPlot(object):
    """ScatterPlot plots a complex scatterplot"""

    def __init__(self, options):
        self.data = {}
        self.label = {}
        self.track_counts = {}
        # inter-plot offset -- may be set from outside
        self.offset = 0

    def add_point(self, key, x, y, label=None):
        if isinstance(x, float):
            x = [x]
        if isinstance(y, float):
            y = [y]
        self.data[key] = x, y
        self.label[key] = label if label is not None else str(key)

    def get_point(self, key):
        return self.data.get(key, None)

    # function to make single-page plot
    def make_figure(self, suptitle=None, save=None, xlabel=None, ylabel=None,
                    offset_std=10, xgrid=False, ygrid=False,
                    figsize=(210, 210), dpi=100, papertype='a4', landscape=False):
        import matplotlib.pyplot as pyplot
        # create Figure object of given size and resolution
        figsize_in = (figsize[0] / 25.4, figsize[1] / 25.4)
        fig = pyplot.figure(figsize=figsize_in, dpi=100)
        # margin sizes. We want to keep them fixed in absolute terms,
        # regardless of plot size. The relative numbers here are good for a 210x210 plot, so
        # we rescale them accordingly.
        mleft = 0.06 * 210. / figsize[0]
        mbottom = 0.06 * 210. / figsize[1]
        mright = 0.01 * 210. / figsize[0]
        mtop = 0.04 * 210. / figsize[1]
        ytitle = 1 - 0.01 * 210. / figsize[1]
        width = (1 - mleft - mright)
        height = (1 - mtop - mbottom)
        plt = fig.add_axes([mleft, mbottom, width, height])
        # remove common prefix from labels
        if self.label:
            lab0 = list(self.label.values())[0]
            prefix = 1
            while len(lab0) >= prefix and all([lab.startswith(lab0[:prefix]) for lab in list(self.label.values())]):
                prefix += 1
            labels = dict([(key, label[prefix - 1:]) for key, label in list(self.label.items())])
        else:
            labels = None
        # plot data over text labels
        for key, (x, y) in list(self.data.items()):
            # plot data
            try:
                plt.plot(x, y, "-x", zorder=0)
            except:
                traceback.print_exc()
                print("Error plotting data for", key)
                continue
            # plot text labels
            if labels.get(key):
                plt.text(x[0], y[0], " " + labels[key], horizontalalignment='left', verticalalignment='top', zorder=1,
                         fontsize=8)
        #    plt.set_xbound(min([min(x) for x,y in self.data.itervalues()])[0],
        #                   max([max(x) for x,y in self.data.itervalues()])[0])
        #    plt.set_ybound(min([min(y) for x,y in self.data.itervalues()])[0],
        #                   max([max(y) for x,y in self.data.itervalues()])[0])
        for lab in plt.get_xticklabels():
            lab.set_fontsize(10)
        for lab in plt.get_yticklabels():
            lab.set_fontsize(10)
            lab.set_rotation('vertical')
        # plot title if asked to
        if xlabel:
            plt.set_xlabel(xlabel, size=10)
        if ylabel:
            plt.set_ylabel(ylabel, size=10)
        if suptitle:
            fig.suptitle(suptitle, y=ytitle, size=10)
        if save:
            fig.savefig(save, papertype=papertype, dpi=dpi,
                        orientation='portrait' if not landscape else 'landscape')
            print("===> Wrote", save)
        return fig


PLOT_SINGLE = 'single'
PLOT_MULTI = 'multi'
PLOT_ERRORBARS = 'errorbars'
PLOT_BARPLOT = 'barplot'


class AbstractBasePlot(object):
    """Abstract base class for the various plotting objects"""

    def __init__(self, options=None, figsize=(290, 210), output_type=None):
        if options is None:
            from optparse import OptionParser, OptionGroup
            parser = OptionParser()
            plotgroup = OptionGroup(parser, "Plotting options")
            outputgroup = OptionGroup(parser, "Output options")
            self.init_options(plotgroup, outputgroup)
            parser.add_option_group(plotgroup)
            parser.add_option_group(outputgroup)
            (options, args) = parser.parse_args([])

        init_pyplot(output_type or getattr(options, 'output_type', 'x11'))
        self.options = options
        if options.width and options.height:
            self._default_figsize = self._user_figsize = (options.width, options.height)
        else:
            self._user_figsize = None
            self._default_figsize = figsize

    @classmethod
    def init_options(self, plotgroup, outputgroup):
        """Initializes options relevant to this plotting class"""
        self._plotgroup = plotgroup
        self._outputgroup = outputgroup

    @classmethod
    def add_plot_option(self, opt, *args, **kw):
        """Helper method: adds a plot option, if it's not already defined"""
        if not self._plotgroup.has_option(opt):
            self._plotgroup.add_option(opt, *args, **kw)

    @classmethod
    def add_output_option(self, opt, *args, **kw):
        """Helper method: adds an output option, if it's not already defined"""
        if not self._outputgroup.has_option(opt):
            self._outputgroup.add_option(opt, *args, **kw)


class MultigridPlot(AbstractBasePlot):
    """MultigridPlot makes NxM plots on a single page"""

    def __init__(self, options, figsize=(290, 210)):
        AbstractBasePlot.__init__(self, options, figsize)
        self._borders = [0.05, 0.99, 0.05, 0.95]
        if options.borders and isinstance(options.borders, str):
            try:
                self._borders = [float(x) for x in options.borders.split(",")]
            except:
                pass
        if options.y_lock:
            try:
                y0, y1 = list(map(float, options.y_lock.split(",")))
            except:
                print("Error parsing --y-lock option, --y-lock MIN,MAX expected.")
                sys.exit(1)
            self.ylock = y0, y1
        else:
            self.ylock = None

    @classmethod
    def init_options(self, plotgroup, outputgroup):
        AbstractBasePlot.init_options(plotgroup, outputgroup)
        self.add_plot_option("--title-fontsize", metavar="POINTS", type="int", default=12,
                             help="Set plot title font size, 0 for no title. Default is %default.")
        self.add_plot_option("--subtitle-fontsize", metavar="POINTS", type="int", default=10,
                             help="Set plot subtitle font size, 0 for none. Default is %default.")
        self.add_plot_option("--label-fontsize", metavar="POINTS", type="int", default=8,
                             help="Set plot label font size, 0 for no labels. Default is %default.")
        self.add_plot_option("--text-fontsize", metavar="POINTS", type="int", default=6,
                             help="Set plot text font size, 0 for no labels. Default is %default.")
        self.add_plot_option("--axis-fontsize", metavar="POINTS", type="int", default=5,
                             help="Set axis label font size, 0 for no axis labels. Default is %default.")
        self.add_plot_option("--text-spacing", metavar="VALUE", type="float", default=0.15,
                             help="Set spacing between lines of text in text comments. Default is %default.")
        self.add_plot_option("--y-ticks", metavar="N", type="int", default=4,
                             help="Maximum number of major ticks along the Y axis. Default is %default.")
        self.add_plot_option("--y-minor-ticks", metavar="INTERVAL", type="float", default=0,
                             help="Add minor ticks along the Y axis at the given intervals.")
        self.add_plot_option("--y-lock", metavar="MIN,MAX", type="str", default="",
                             help="Locks the Y scale of the plot to the given min/max values.")
        self.add_plot_option("--borders", metavar="L,R,B,T", type="str",
                             help="Set left/right/bottom/top plot borders, in normalized page coordinates. Default is %default.")
        self.add_plot_option("--subplot-wspace", metavar="W", type="float", default=0,
                             help="Set spacing between subplots, as fraction of subplot width. Default is automatic.")
        self.add_plot_option("-S", "--subtitle", type="str", default="",
                             help="Subtitle for plot, added (in parentheses) after plot title")

        self.add_output_option("-o", "--output-type", metavar="TYPE", type="string", default="png",
                               help="File format, see matplotlib documentation "
                                    "for supported formats. At least 'png', 'pdf', 'ps', 'eps' and 'svg' are supported, or use 'x11' to display "
                                    "plots interactively. Default is '%default.'")
        self.add_output_option("-r", "--refresh", action="store_true",
                               help="Refresh plots even if they already exist (default is to keep existing plots.)")
        self.add_output_option("--output-prefix", metavar="PREFIX", type="string", default="",
                               help="Prefix output filenames with PREFIX_")
        self.add_output_option("--papertype", dest="papertype", type="string", default="a4",
                               help="set paper type (for .ps output only.) Default is '%default', but can also use e.g. 'letter', 'a3', etc.")
        self.add_output_option("-W", "--width", type="int", default=0,
                               help="set explicit plot width, in mm. (Useful for .eps output)")
        self.add_output_option("-H", "--height", type="int", default=0,
                               help="set explicit plot height, in mm. (Useful for .eps output)")
        self.add_output_option("--dpi", type="int", default=300,
                               help="figure resolution. Default is %default.")
        self.add_output_option("--scale", type="float", default=1,
                               help="scale plot sizes up by the given factor.")
        self.add_output_option("--portrait", action="store_true",
                               help="Force portrait orientation. Default is to select orientation based on plot size")
        self.add_output_option("--landscape", action="store_true",
                               help="Force landscape orientation.")

    @staticmethod
    def get_plot_data(data, xaxis=None):
        """Helper function, interprets the plot data returned by a datafunc.
        'data' is input data as returned by the datafunc argument to make_figure() below.
        'xaxis' is default X axis, Note to use ordinal numbering.
        Return value is tuple of X,Y,Yerr. Yerr is None if no error bars are provided.
        """
        # a 2- or 3-tuple is interpreted as x,y[,yerr]
        if isinstance(data, tuple):
            if len(data) == 2:
                x, y = data
                yerr = None
            elif len(data) == 3:
                x, y, yerr = data
            else:
                raise TypeError("incorrect datum returned: 2- or 3-tuple expected, got %d-tuple""" % len(data))
        # else interpret data as Y
        else:
            y = data
            x = yerr = None
        # set X axis
        if x is None:
            x = list(range(len(y))) if xaxis is None else xaxis
        # check lengths
        if len(x) != len(y):
            raise TypeError("misshaped datum returned: %d X elements, %d Y elements""" % (len(x), len(y)))
        if yerr is not None and len(yerr) != len(y):
            raise TypeError("misshaped datum returned: %d Y elements, %d Yerr elements""" % (len(y), len(yerr)))
        return x, y, yerr

    def make_figure(self, rows, cols,  # (irow,row) and (icol,col) list
                    datafunc,  # datafunc(irow,icol) returns plot data for plot i,j
                    mode=PLOT_SINGLE, xaxis=None,
                    hline=None,  # plot horizontal line at Y position, None for none
                    ylock="row",
                    # lock Y scale. "row" locks across rows, "col" locks across columns, True locks across whole plot, (ymin,ymax) sets an explicit scale
                    mean_format="%.2f",
                    suptitle=None,  # title of plot
                    save=None,  # filename to save to
                    format=None,  # format: use options.output_type by default
                    figsize=None  # figure width,height in mm
                    ):
        from matplotlib import ticker
        if save and (format or self.options.output_type.upper()) != 'X11':
            save = "%s.%s" % (save, format or self.options.output_type)
            if self.options.output_prefix:
                save = "%s_%s" % (self.options.output_prefix, save)
            # exit if figure already exists, and we're not refreshing
            if os.path.exists(
                    save) and not self.options.refresh:  # and os.path.getmtime(save) >= os.path.getmtime(__file__):
                print(save, "exists, not redoing")
                return save
        else:
            save = None

        figsize = numpy.array(self._user_figsize or figsize or self._default_figsize) / 25.4
        figsize *= self.options.scale
        fig = pyplot.figure(figsize=figsize, dpi=10)
        iplot = 0
        rows = list(rows)
        cols = list(cols)
        if ylock and not self.ylock:
            dummy = numpy.array([0.])
            # form up ymin, ymax: NROWxNCOL arrays of min/max values per each plot
            if mode is PLOT_SINGLE or mode is PLOT_BARPLOT:
                datafilter = lambda x: x if isinstance(x, numpy.ndarray) else dummy
                ymin = numpy.array([[datafilter(datafunc(row[0], col[0])).min() for col in cols] for row in rows])
                ymax = numpy.array([[datafilter(datafunc(row[0], col[0])).max() for col in cols] for row in rows])
            elif mode is PLOT_MULTI:
                ymin = numpy.array(
                    [[min([self.get_plot_data(dd, xaxis)[1].min() for dd in datafunc(row[0], col[0])]) for col in cols]
                     for row in rows])
                ymax = numpy.array(
                    [[max([self.get_plot_data(dd, xaxis)[1].max() for dd in datafunc(row[0], col[0])]) for col in cols]
                     for row in rows])
            elif mode is PLOT_ERRORBARS:
                dummy = numpy.array([0.])
                datafilter1 = lambda x: x[0] - x[1] if isinstance(x[0], numpy.ndarray) else dummy
                datafilter2 = lambda x: x[0] + x[1] if isinstance(x[0], numpy.ndarray) else dummy
                ymin = numpy.array([[datafilter1(datafunc(row[0], col[0])).min() for col in cols] for row in rows])
                ymax = numpy.array([[datafilter2(datafunc(row[0], col[0])).max() for col in cols] for row in rows])
            # now collapse them according to ylock mode
            if isinstance(ylock, (list, tuple)) and len(ylock) == 2:
                ymin[...], ymax[...] = ylock
            elif ylock == "row":
                ymin[...] = ymin.min(1)[:, numpy.newaxis]
                ymax[...] = ymax.max(1)[:, numpy.newaxis]
            elif ylock == "col":
                ymin[...] = ymin.min(0)[numpy.newaxis, :]
                ymax[...] = ymax.max(0)[numpy.newaxis, :]
            else:
                ymin[...] = ymin.min()
                ymax[...] = ymax.max()
        # make legend across the top
        for icol, col in cols:
            iplot += 1
            plt = fig.add_subplot(len(rows) + 1, len(cols) + 1, iplot)
            plt.axis("off")
            plt.text(.5, 0, col, fontsize=self.options.subtitle_fontsize * self.options.scale,
                     horizontalalignment='center', verticalalignment='bottom')
        iplot += 1

        # helper function to put lines of text instead of real plot
        def plot_text(plt, lines):
            y0 = 0.95
            for i, text in enumerate(lines):
                plt.text(0, y0 - self.options.text_spacing * i, text,
                         fontsize=self.options.text_fontsize * self.options.scale, horizontalalignment='left',
                         verticalalignment='top')
            plt.axis("off")
            plt.set_ylim(0, 1)

        # now make plots
        for irow, row in rows:
            for icol, col in cols:
                iplot += 1
                plt = fig.add_subplot(len(rows) + 1, len(cols) + 1, iplot)
                # plt.axis("off")
                plt.yaxis.set_major_locator(ticker.MaxNLocator(self.options.y_ticks))
                if self.options.y_minor_ticks:
                    plt.yaxis.set_minor_locator(ticker.MultipleLocator(self.options.y_minor_ticks))
                if ylock and ylock != "col" and icol:
                    plt.set_yticklabels([])
                plt.set_xticks([])
                data = datafunc(irow, icol)
                realplot = True
                # if data is a sequence of strings, plot them
                if isinstance(data, (list, tuple)) and isinstance(data[0], str):
                    plot_text(plt, data)
                    realplot = False
                elif mode is PLOT_SINGLE:
                    plt.plot(list(range(len(data))) if xaxis is None else xaxis, data)
                    if xaxis is None:
                        plt.set_xlim(-1, len(data))
                elif mode is PLOT_BARPLOT:
                    plt.bar(list(range(len(data))) if xaxis is None else xaxis, data, linewidth=0)
                    if xaxis is None:
                        plt.set_xlim(-1, len(data))
                elif mode is PLOT_MULTI:
                    for dd in data:
                        x, y, yerr = self.get_plot_data(dd, xaxis)
                        if yerr is None:
                            plt.plot(x, y)
                        else:
                            plt.errorbar(x, y, yerr, fmt=None, capsize=1)
                        if x == list(range(len(y))):
                            plt.set_xlim(-1, len(y))
                elif mode is PLOT_ERRORBARS:
                    y, yerr = data
                    if len(y) == 1:
                        plot_text(plt, (mean_format % y[0], "+/-" + (mean_format % yerr[0])))
                        realplot = False
                    else:
                        plt.errorbar(list(range(len(y))) if xaxis is None else xaxis, y, yerr, fmt=None, capsize=1)
                        if xaxis is None:
                            plt.set_xlim(-1, len(y))
                if realplot:
                    if self.ylock:
                        plt.set_ylim(*self.ylock)
                        if icol:
                            plt.set_yticklabels([])
                    elif ylock:
                        plt.set_ylim(ymin[irow, icol], ymax[irow, icol])
                        if ylock is not "col" and icol:
                            plt.set_yticklabels([])
                    if hline is not None:
                        plt.axhline(y=hline, color='black')
                for lab in plt.get_yticklabels():
                    lab.set_fontsize(self.options.axis_fontsize * self.options.scale)
            # add row label
            iplot += 1
            plt = fig.add_subplot(len(rows) + 1, len(cols) + 1, iplot)
            plt.text(0, .5, row, fontsize=self.options.label_fontsize * self.options.scale, horizontalalignment='left',
                     verticalalignment='center')
            plt.axis("off")
        # adjust layout
        fig.subplots_adjust(left=self._borders[0], right=self._borders[1],
                            top=self._borders[3], bottom=self._borders[2],
                            wspace=self.options.subplot_wspace or None)
        # plot title if asked to
        if suptitle and self.options.title_fontsize:
            if self.options.subtitle:
                suptitle = "%s (%s)" % (suptitle, self.options.subtitle)
            fig.suptitle(suptitle, fontsize=self.options.title_fontsize * self.options.scale)
        if save:
            if self.options.portrait:
                orientation = 'portrait'
            elif self.options.landscape:
                orientation = 'landscape'
            else:
                orientation = 'portrait' if figsize[0] < figsize[1] else 'landscape'
            fig.savefig(save, papertype=self.options.papertype, dpi=self.options.dpi,
                        orientation=orientation)
            print("Wrote", save, "in", orientation, "orientation")
            fig = None
            pyplot.close("all")
        return save


class SkyPlot(AbstractBasePlot):
    """SkyPlot plots things in an l,m plot"""

    def __init__(self, options, figsize=(290, 210)):
        AbstractBasePlot.__init__(self, options, figsize)
        self._borders = [.05, .9, .05, .9]
        if options.circle_borders and isinstance(options.circle_borders, str):
            try:
                self._borders = [float(x) for x in options.circle_borders.split(",")]
            except:
                pass

    @classmethod
    def init_options(self, plotgroup, outputgroup):
        AbstractBasePlot.init_options(plotgroup, outputgroup)
        self.add_plot_option("--radius", type="int", default=0,
                             help="radius of circle plots, in arcmin. Default is auto-sizing.")
        self.add_plot_option("-G", "--grid-circle", metavar="ARCMIN", type="int", action="append",
                             help="sets radius of grid circle(s) in circle plots. May be given multiple times. Default is 30 and 60.")
        self.add_plot_option("--circle-title-fontsize", metavar="POINTS", type="int", default=12,
                             help="Set plot title font size in circle plots, 0 for no title. Default is %default.")
        self.add_plot_option("--circle-label-fontsize", metavar="POINTS", type="int", default=8,
                             help="Set plot label font size in circle plots, 0 for no labels. Default is %default.")
        self.add_plot_option("--circle-axis-fontsize", metavar="POINTS", type="int", default=5,
                             help="Set axis label font size in circle_plots, 0 for no axis labels. Default is %default.")
        self.add_plot_option("--names-only", action="store_true",
                             help="Use only source names for labels in circle plots. Default is name, mean and std values.")
        self.add_plot_option("--circle-label-offset", metavar="FRAC", type="float", default=0.02,
                             help="Vertical offset of circle labels. Use 0 to center labels on source position. Default is %default.")
        self.add_plot_option("--circle-borders", metavar="L,R,B,T", type="str",
                             help="Set left/right/bottom/top circle-plot borders, in normalized page coordinates. Default is %default.")
        self.add_plot_option("-S", "--subtitle", type="str", default="",
                             help="Subtitle for plot, added (in parentheses) after plot title")

        self.add_output_option("-o", "--output-type", metavar="TYPE", type="string", default="png",
                               help="File format, see matplotlib documentation "
                                    "for supported formats. At least 'png', 'pdf', 'ps', 'eps' and 'svg' are supported, or use 'x11' to display "
                                    "plots interactively. Default is '%default.'")
        self.add_output_option("-r", "--refresh", action="store_true",
                               help="Refresh plots even if they already exist (default is to keep existing plots.)")
        self.add_output_option("--output-prefix", metavar="PREFIX", type="string", default="",
                               help="Prefix output filenames with PREFIX_")
        self.add_output_option("--papertype", dest="papertype", type="string", default="a4",
                               help="set paper type (for .ps output only.) Default is '%default', but can also use e.g. 'letter', 'a3', etc.")
        self.add_output_option("-W", "--width", type="int",
                               help="set explicit plot width, in mm. (Useful for .eps output)")
        self.add_output_option("-H", "--height", type="int",
                               help="set explicit plot height, in mm. (Useful for .eps output)")
        self.add_output_option("--dpi", type="int", default=300,
                               help="figure resolution. Default is %default.")
        self.add_output_option("--portrait", action="store_true",
                               help="Force portrait orientation. Default is to select orientation based on plot size")
        self.add_output_option("--landscape", action="store_true",
                               help="Force landscape orientation.")

    """Plotting function"""

    def make_figure(self, ll, mm,
                    markers,  # marker is a list of strings or dicts or tuples.
                    # For each marker, if str: axes.plot(x,y,str) is called
                    # if dict: axes.plot(x,y,**dict)
                    # else tuple is func,arg,kwargs
                    # and axes.func(*args,**kwargs) is called
                    labels=None,
                    radius=None,  # plot radius. If None, then it's set automatically
                    zero_lines=True,  # plot l=0 and m=0 lines
                    suptitle=None,  # title of plot
                    save=None,  # filename to save to
                    format=None,  # format: use self.options.output_type by default
                    figsize=(210, 210),  # figure width,height in mm
                    ):
        if save and (format or self.options.output_type.upper()) != 'X11':
            save = "%s.%s" % (save, format or self.options.output_type)
            if self.options.output_prefix:
                save = "%s_%s" % (self.options.output_prefix, save)
            # exit if figure already exists, and we're not refreshing
            if os.path.exists(
                    save) and not self.options.refresh:  # and os.path.getmtime(save) >= os.path.getmtime(__file__):
                print(save, "exists, not redoing")
                return save
        else:
            save = None

        figsize = numpy.array(figsize or self._default_figsize) / 25.4
        fig = pyplot.figure(figsize=figsize, dpi=600)
        plt = fig.add_axes([self._borders[0], self._borders[2], self._borders[1] - self._borders[0],
                            self._borders[3] - self._borders[2]])
        if zero_lines:
            plt.axhline(y=0, color='black', linestyle=':')
            plt.axvline(x=0, color='black', linestyle=':')

        # if ll is None, use length of markers for everything
        if ll is None:
            ll = mm = numpy.zeros(len(markers), float)

        if not isinstance(markers, (list, tuple)):
            markers = [markers] * len(ll)

        if not labels:
            labels = [None] * len(ll)

        # set limits
        if self.options.radius:
            lim = self.options.radius
        elif radius:
            lim = radius
        else:
            lim = max(max(abs(ll)), max(abs(mm))) * 1.05 / ARCMIN

        # plot markers
        for l, m, mark, label in zip(ll / ARCMIN, mm / ARCMIN, markers, labels):
            if isinstance(mark, str):
                plt.plot(l, m, mark)
            elif isinstance(mark, dict):
                plt.plot(l, m, **mark)
            elif isinstance(mark, tuple):
                func, args, kwargs = mark
                getattr(plt, func)(*args, **kwargs)
            if label and self.options.circle_label_fontsize:
                plt.text(l, m - lim * self.options.circle_label_offset,
                         label, fontsize=self.options.circle_label_fontsize,
                         horizontalalignment='center',
                         verticalalignment='top' if self.options.circle_label_offset else 'center')

        x = numpy.cos(numpy.arange(0, 361) * math.pi / 180)
        y = numpy.sin(numpy.arange(0, 361) * math.pi / 180)

        gc = self.options.grid_circle or (30, 60)
        for R in gc:
            plt.plot(x * R, y * R, linestyle=':', color='black')

        plt.set_xlim(lim, -lim)
        plt.set_ylim(-lim, lim)
        for lab in list(plt.get_xticklabels()) + list(plt.get_yticklabels()):
            lab.set_fontsize(self.options.circle_axis_fontsize)

        if suptitle and self.options.title_fontsize:
            if self.options.subtitle:
                suptitle = "%s (%s)" % (suptitle, self.options.subtitle)
            fig.suptitle(suptitle, fontsize=self.options.circle_title_fontsize)
        if save:
            if self.options.portrait:
                orientation = 'portrait'
            elif self.options.landscape:
                orientation = 'landscape'
            else:
                orientation = 'portrait' if figsize[0] < figsize[1] else 'landscape'
            fig.savefig(save, papertype=self.options.papertype, dpi=self.options.dpi,
                        orientation=orientation)
            print("Wrote", save, "in", orientation, "orientation")
            fig = None
            pyplot.close("all")
        return save

class LSTElevationPlot(AbstractBasePlot):
    """Plots tracks of sky objects in an LST vs elevation plot"""

    def __init__(self, options=None, figsize=(200, 100), output_type='x11'):
        AbstractBasePlot.__init__(self, options, figsize, output_type)
        self._borders = [.1, .9, .1, .85]
        self._default_figsize = (200,200)

    @classmethod
    def init_options(self, plotgroup, outputgroup):
        AbstractBasePlot.init_options(plotgroup, outputgroup)
        self.add_plot_option("--title-fontsize", metavar="POINTS", type="int", default=8,
                             help="Set plot title font size, 0 for no title. Default is %default.")
        # self.add_plot_option("--label-fontsize", metavar="POINTS", type="int", default=8,
        #                      help="Set plot label font size in circle plots, 0 for no labels. Default is %default.")
        self.add_plot_option("--tick-fontsize", metavar="POINTS", type="int", default=5,
                             help="Set axis tick label font size, 0 for no tick labels. Default is %default.")
        self.add_plot_option("--axis-fontsize", metavar="POINTS", type="int", default=5,
                             help="Set axis label font size, 0 for no axis labels. Default is %default.")
        self.add_plot_option("--scan-fontsize", metavar="POINTS", type="int", default=8,
                             help="Set scan label font size, 0 for no scan labels. Default is %default.")
        self.add_plot_option("--legend-fontsize", metavar="POINTS", type="int", default=5,
                             help="Set legend label font size, 0 for no legend labels. Default is %default.")
        self.add_plot_option("--marker-size", metavar="POINTS", type="int", default=2,
                             help="Set plot marker size. Default is %default.")
        self.add_plot_option("-S", "--subtitle", type="str", default="",
                             help="Subtitle for plot, added (in parentheses) after plot title")
        self.add_output_option("--papertype", dest="papertype", type="string", default="a4",
                               help="set paper type (for .ps output only.) Default is '%default', but can also use e.g. 'letter', 'a3', etc.")
        self.add_output_option("-W", "--width", type="int",
                               help="set explicit plot width, in mm. (Useful for .eps output)")
        self.add_output_option("-H", "--height", type="int",
                               help="set explicit plot height, in mm. (Useful for .eps output)")
        self.add_output_option("--dpi", type="int", default=300,
                               help="figure resolution. Default is %default.")
        self.add_output_option("--portrait", action="store_true",
                               help="Force portrait orientation. Default is to select orientation based on plot size")
        self.add_output_option("--landscape", action="store_true",
                               help="Force landscape orientation.")

    """Plotting function"""

    def make_figure(self, field_time, field_radec, obs_xyz, scan_starts,
                    suptitle=None,  # title of plot
                    save=None,  # filename to save to
                    display=None, # True if display is on
                    figsize=(200, 200),  # figure width,height in mm
                    ):
        fields = list(field_radec.keys())

        timestamps = set()
        for field in fields:
            timestamps.update(field_time[field])
        tm_uniq = sorted(timestamps)  # unique timestamps

        import numpy as np
        from casacore.measures import dq
        dm = casacore.measures.measures()


        # create an epoch measure per each unique timestamp
        tm_epoch = { t: dm.epoch('utc', dq.quantity(t, 's')) for t in tm_uniq }

        tm_d = { t: dm.getvalue(ep)[0].get_value('d') for t,ep in tm_epoch.items() }

        # create a direction measure per each field
        dir_field = { field: dm.direction('j2000', dq.quantity(ra, 'rad'), dq.quantity(dec, 'rad'))
                      for field,(ra,dec) in field_radec.items() }

        fields += ["Sun"]
        dir_field['Sun'] = dm.direction("SUN")
        field_time['Sun'] = tm_uniq

        # get starting date
        import datetime
        start_day = datetime.datetime.fromtimestamp(dq.quantity(tm_uniq[0], 's').to_unix_time()).strftime("%d %b %Y")

        # create a position measure for the antenna and set it as the frame
        pos_ant0 = dm.position('itrf', *[dq.quantity(x, 'm') for x in obs_xyz])
        dm.doframe(pos_ant0)

        # convert to LST
        tm_lst0 = dm.getvalue(dm.measure(tm_epoch[tm_uniq[0]], 'LAST'))[0].get_value('d')
        offset = tm_lst0 - tm_d[tm_uniq[0]]  # offset from LST to UCT
        def lst2utc(lst):
            return lst - offset
        def utc2lst(utc):
            return utc + offset

        # get AZ/EL per each field and its timeslots
        field_lst_azel = {}
        scan_labels_el = {}
        scan_labels_az = {}

        for field in fields:
            num_ts = len(field_time[field])
            field_azel = np.zeros((num_ts, 3), float)
            for its, t in enumerate(field_time[field]):
                dm.doframe(pos_ant0)
                dm.doframe(tm_epoch[t])
                azel_meas = dm.measure(dir_field[field], 'AZELGEO')
                azel_quant = dm.get_value(azel_meas)
                field_azel[its, 0] = x = tm_d[t]
                field_azel[its, 1] = az = azel_quant[0].get_value('deg')
                field_azel[its, 2] = el = azel_quant[1].get_value('deg')
                scan = scan_starts.get(t, None)
                if scan is not None and field != "Sun":
                    scan_labels_el[x, el] = str(scan)
                    scan_labels_az[x, az] = str(scan)
            mask = np.zeros_like(field_azel, bool)
            mask[:] = (field_azel[:, 2]<0)[:, np.newaxis]
            field_lst_azel[field] = numpy.ma.masked_array(field_azel, mask[np.newaxis, np.newaxis, :])
        figsize = numpy.array(figsize or self._default_figsize)/25.4
        # fig, (plt, plt_az) = pyplot.figure(figsize=figsize, dpi=self.options.dpi)
        # plt = fig.add_axes([self._borders[0], self._borders[2], self._borders[1] - self._borders[0],
        #                     self._borders[3] - self._borders[2]])
        fig, (plt, plt_az) = pyplot.subplots(2, sharex=True, figsize=figsize, dpi=self.options.dpi)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        for field in fields:
            plt_az.plot(field_lst_azel[field][:,0], field_lst_azel[field][:,1], '.', ms=self.options.marker_size)
        plt_az.yaxis.set_tick_params(labelsize=self.options.tick_fontsize)

        for field in fields:
            plt.plot(field_lst_azel[field][:,0], field_lst_azel[field][:,2], '.', ms=self.options.marker_size, label=field)

        if self.options.scan_fontsize:
            for (x, y), label in scan_labels_el.items():
                plt.text(x, y, label, horizontalalignment='left', verticalalignment='bottom',
                         fontdict=dict(fontsize=self.options.scan_fontsize))
            for (x, y), label in scan_labels_az.items():
                plt_az.text(x, y, label, horizontalalignment='left', verticalalignment='bottom',
                         fontdict=dict(fontsize=self.options.scan_fontsize))
        plt_az.xaxis.set_major_locator(matplotlib.dates.HourLocator())
        plt_az.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        plt_az.xaxis.set_tick_params(labelsize=self.options.tick_fontsize)
        plt.yaxis.set_tick_params(labelsize=self.options.tick_fontsize)
        secax = plt.secondary_xaxis('top', functions=(utc2lst, lst2utc))
        secax.set_xlabel("LST", fontdict=dict(fontsize=self.options.axis_fontsize))
        secax.xaxis.set_major_locator(matplotlib.dates.HourLocator())
        secax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        secax.xaxis.set_tick_params(labelsize=self.options.tick_fontsize)
        if self.options.axis_fontsize:
            pyplot.xlabel('UTC ({})'.format(start_day), fontdict=dict(fontsize=self.options.axis_fontsize))
            plt.set_ylabel('Elevation, deg', fontdict=dict(fontsize=self.options.axis_fontsize))
            if plt_az:
                plt_az.set_ylabel('Azimuth, deg', fontdict=dict(fontsize=self.options.axis_fontsize))
        if self.options.legend_fontsize:
            plt.legend(loc='best', #bbox_to_anchor=(1, 0.5),
                              fontsize=self.options.legend_fontsize)

        # for lab in list(plt.get_xticklabels()) + list(plt.get_yticklabels()) + list(secax.get_xticklabels()):
        #     lab.set_fontsize(self.options.tick_fontsize)
        #
        if suptitle and self.options.title_fontsize:
            if self.options.subtitle:
                suptitle = "%s (%s)" % (suptitle, self.options.subtitle)
            fig.suptitle(suptitle, fontsize=self.options.title_fontsize)
        if save:
            if self.options.portrait:
                orientation = 'portrait'
            elif self.options.landscape:
                orientation = 'landscape'
            else:
                orientation = 'portrait' if figsize[0] < figsize[1] else 'landscape'
            fig.savefig(save, papertype=self.options.papertype, dpi=self.options.dpi,
                        orientation=orientation)
            print("Wrote", save, "in", orientation, "orientation")
            fig = None
            pyplot.close("all")

        return fig

    @staticmethod
    def load_ms_fields(msname, verbose=1):
        from casacore.tables import table
        from collections import OrderedDict
        import numpy as np

        ms = table(msname, ack=False)
        timestamps = ms.getcol("TIME")
        field_id = ms.getcol("FIELD_ID")
        scan = ms.getcol("SCAN_NUMBER")

        verbose and print("Fields in {}:".format(msname))

        tm_uniq = np.array(sorted(set(timestamps)))          # unique timestamps
        timeslot_index =  { t:i for i, t in enumerate(tm_uniq)}   # timestamp -> timeslot number
        ts = [timeslot_index[t] for t in timestamps]   # timeslot[row]
        num_ts = len(tm_uniq)
        import bisect
        first_row_per_ts = [bisect.bisect_left(ts, its) for its in range(num_ts)]  # first row of any given timestamp
        time_ts = np.array([timestamps[row0] for row0 in first_row_per_ts])        # time per timestamp
        fid_ts = np.array([field_id[row0] for row0 in first_row_per_ts])           # field_id per timestamp
        scan_ts = np.array([scan[row0] for row0 in first_row_per_ts])              # scan per timestamp

        # figure out starting timeslots of scans
        scan_start = np.zeros_like(scan_ts, bool)
        scan_start[0] = True
        scan_start[1:] = (scan_ts[1:] != scan_ts[:-1])
        scan_start_ts = {time_ts[ts]: scan for ts, scan in enumerate(scan_ts) if scan_start[ts]}

        ftab = table(msname+"::FIELD", ack=False)
        field_names = ftab.getcol("NAME")
        field_dirs = ftab.getcol("PHASE_DIR")

        field_time = OrderedDict()
        field_radec = OrderedDict()

        for fid, field in enumerate(field_names):
            field_name = "{}: {}".format(fid, field)
            field_radec[field_name] = fdir = field_dirs[fid,0]
            field_time[field_name] = ftime = time_ts[fid_ts == fid]
            verbose and print("    [{:2}] {:16s} {:8.2f}d {:+8.2f}d     {} timeslots".format(
                                fid, field, fdir[0]*180/np.pi, fdir[1]*180/np.pi, len(ftime)
                              ))

        anttab = table(msname+"::ANTENNA", ack=False)
        ant_xyz = anttab.getcol("POSITION", 0 , 1)[0]
        verbose and print("  Antenna 0 ITRF position is {} {} {}".format(*ant_xyz))

        return field_time, field_radec, ant_xyz, scan_start_ts
