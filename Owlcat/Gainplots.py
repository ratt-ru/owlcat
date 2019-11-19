#!/usr/bin/env python

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

import matplotlib
matplotlib.use('Agg')
import pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from six.moves import cPickle
import numpy
import math
import cmath
from past.builtins import cmp
from functools import cmp_to_key
from future.utils import raise_with_traceback


def _normifrgain(rr):
    """Converts gains to offsets with std"""
    if type(rr) in (float, complex):
        return abs(rr), 0
    else:
        offset = abs(rr[rr != 1])
        return float(offset.mean()), float(offset.std())


def _complexifrgain(rr):
    """Converts gains to complex offsets with std"""
    if type(rr) in (float, complex):
        return rr, 0
    else:
        vals = rr[rr != 1]
        offset = float(abs(vals).mean())
        mean = vals.mean().ravel()[0]
        mean = cmath.rect(offset, cmath.phase(mean))
        return mean, float((vals - mean).std())


def _is_unity(rr, ll):
    return type(rr) in (float, complex) and rr == 1 and type(ll) in (float, complex) and ll == 1


def _cmp_antenna(sa, sb):
    """Helper function to sort antenna names. Try numeric compare first, fall back to text compare if failed"""
    try:
        return cmp(int(sa), int(sb))
    except:
        return cmp(sa, sb)


def make_gain_plots(filename, prefix=None, gain_label="G",
                    ylim=None, ylim_offdiag=None,
                    ant=None, outname=None,
                    feed_type="rl",
                    gain_plot_style="dot"):
    """Makes a set of gain plots from the specified saved file ('filename').
    'ylim' can be used to override GAIN_PLOT_AMPL_YLIM setting.
    'ant' can be set to a whitespace-separated list of antennas. Wildcard patterns are allowed."""

    print(("loading gain solutions from %s" % filename))
    G0 = cPickle.load(open(filename))

    G = G0['gains'][gain_label]['solutions']
    solkeys = list(G.keys())
    # solutions are stored either as [antenna,{0,1}] for diagonal Jones, or [antenna][{0,1,2,3}] for 2x2 Jones
    diagonal = all([type(k) is tuple and len(k) == 2 and k[1] in (0, 1) for k in solkeys])
    if diagonal:
        antennas = antennas0 = sorted(set([k[0] for k in solkeys]), key=cmp_to_key(_cmp_antenna))
    else:
        antennas = antennas0 = sorted(solkeys, key=cmp_to_key(_cmp_antenna))

    NANT = len(antennas)
    print(("loaded solutions for %d antennas, diagonal=%s" % (NANT, diagonal)))

    ylim1 = ylim_offdiag

    if ant:
        antpatts = str(ant).split()
        antennas = [x for x in antennas if any([fnmatch.fnmatch(str(x), patt) for patt in antpatts])]
        print(("applied subset '%s' to antennas %s, selecting %d" % (nant, " ".join(map(str, antennas0)), NANT)))

    if not antennas:
        print("no antennas to plot")
        return

    # feed labels
    feeds = ("RR", "RL", "LR", "LL") if feed_type.upper() == "RL" else ("XX", "XY", "YX", "YY")

    ncols = 4 if diagonal else 8
    nrows = len(antennas)

    # find amplitude axis ranges for parallel and cross-hand plots
    ppminmax = 1e+99, -1e+99
    xhminmax = 1e+99, -1e+99
    for ant in antennas:
        for j in ((0, 1) if diagonal else list(range(4))):
            if diagonal:
                gg = G[ant, j]
                xhand = False
            else:
                gg = G[ant][j]
                xhand = j in (1, 2)
            valid = (gg != 0) & (gg != 1);  # mask trivial or unfilled solutions
            # TODO: Find out what's causing these entries to be scalers, and fix the issue properly
            if isinstance(gg, numpy.ndarray) and valid.any():
                a = abs(gg)
                amin, amax = a[valid].min(), a[valid].max()
                if not xhand:
                    ppminmax = (min(ppminmax[0], amin), max(ppminmax[1], amax))
                else:
                    xhminmax = (min(xhminmax[0], amin), max(xhminmax[1], amax))

    for xaxis, yaxis, label in (0, 1, "timeslot"), (1, 0, "chan"):
        pylab.figure(figsize=(5 * ncols, 3 * nrows))
        for row, ant in enumerate(antennas):
            for icol, j in enumerate([0, 1] if diagonal else [0, 3, 1, 2]):
                jonesterm = [0, 3][j] if diagonal else j
                feed = feeds[j]
                # TODO: Find out what's causing these entries to be scalers, and fix the issue properly
                if not isinstance(G[ant, j] if diagonal else G[ant][j], numpy.ndarray):
                    continue
                gg = numpy.transpose(G[ant, j] if diagonal else G[ant][j], (xaxis, yaxis))
                # get plot axis and averaging axis
                nx, ny = gg.shape
                x = numpy.zeros((nx, ny))
                x[...] = numpy.arange(nx)[:, numpy.newaxis]
                ax = pylab.subplot(nrows, ncols, row * ncols + icol * 2 + 1)
                pylab.subplots_adjust(bottom=.01, left=.01, top=.99, right=.99)
                valid = (gg != 0) & (gg != 1);  # mask trivial or unfilled solutions
                amp = abs(gg)
                ampwh = numpy.ma.masked_array(amp, mask=~valid)
                amid = ampwh.mean(1, dtype='float64')
                xvalid = ~amid.mask
                if gain_plot_style == 'fill':
                    pylab.fill_between(x[:, 0], ampwh.min(1), ampwh.max(1), where=xvalid, color='grey')
                    pylab.plot(x[xvalid, 0], amid[xvalid], '.', mec='blue', mfc='blue')
                else:
                    pylab.plot(x[valid], amp[valid], '.', ms=0.5, mec='grey', mfc='grey')
                    pylab.plot(x[xvalid, 0], amid[xvalid], '-', ms=0.5, mec='blue', mfc='blue', color='blue')
                pylab.title("%s:%s:ampl - %s" % (ant, feed, label))
                tickstep = 10 ** int(math.log10(nx) - 1)
                labstep = 10 ** int(math.log10(nx))
                if nx / labstep < 5:
                    labstep /= 2.
                ax.xaxis.set_major_locator(MultipleLocator(labstep))
                ax.xaxis.set_minor_locator(MultipleLocator(tickstep))
                ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
                ax.xaxis.set_tick_params(which='major', length=6)
                ax.xaxis.set_tick_params(which='minor', length=3)

                pylab.xlim(-1, nx)
                pylab.ylim(*((ylim or ppminmax) if jonesterm in (0, 3) else (ylim1 or xhminmax)))

                ax = pylab.subplot(nrows, ncols, row * ncols + icol * 2 + 2)
                ph0 = numpy.angle(gg) * 180 / math.pi
                pylab.plot(x[valid], ph0[valid], '.', ms=0.5, mec='grey', mfc='grey')
                ph = ph0[:, ny / 2]
                xvalid = valid[:, ny / 2]

                if xvalid.any():
                    pylab.plot(x[xvalid, 0], ph[xvalid], '-', ms=0.5, mec='blue', mfc='blue', color='blue')
                pylab.title("%s:%s:phase (deg) - %s" % (ant, feed, label))
                pylab.xlim(-1, nx)
                ax.xaxis.set_major_locator(MultipleLocator(labstep))
                ax.xaxis.set_minor_locator(MultipleLocator(tickstep))
                ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
                ax.xaxis.set_tick_params(which='major', length=6)
                ax.xaxis.set_tick_params(which='minor', length=3)

        prefix = prefix or filename + "_gains"
        outname = "%s-%s.png" % (prefix, label)
        pylab.savefig(outname, dpi=75)
        print(("generated plot %s" % outname))


def make_diffgain_plots(filename, prefix=None, ylim=None, subset=None, ant=None, feed_type="rl"):
    """Makes a set of diffgain plots from the specified saved file ('filename'). Plots are placed into directory
    'dir' and filenames are prefixed by 'prefix'.
    'ylim' can be used to override DIFFGAIN_PLOT_AMPL_YLIM setting.
    'subset' can be set to a whitespace-separated list of parameter names
    'ant' can be set to a whitespace-separated list of antennas. Wildcard patterns are allowed in both cases"""

    print(("loading diffgain solutions from %s" % filename))
    DG0 = cPickle.load(open(filename))

    DG = DG0['gains']
    srcnames = sorted(DG.keys())
    from past.builtins import cmp
    from functools import cmp_to_key
    antennas = sorted(list(DG[srcnames[0]]['solutions'].keys()), key=cmp_to_key(_cmp_antenna))

    NANT = len(antennas)
    prefix = prefix or filename + "_diffgains"

    # apply subsets
    if subset:
        src = subset.split()
        print("applying subset '%s' to sources" % (src, srcnames))
        srcnames = [x for x in srcnames if any([fnmatch.fnmatch(x, patt) for patt in src])]
    if ant:
        ant = ant.split()
        print("applying subset '%s' to antennas" % (ant, antennas))
        antennas = [x for x in antennas if any([fnmatch.fnmatch(x, patt) for patt in ant])]

    if not srcnames or not antennas:
        print("no parameters or antennas to plot")
        return

    ncols = len(srcnames)

    print()
    "making diffgain plots for", srcnames
    print(("and %d antennas" % NANT))
    # feed labels
    feeds = ("RR", "RL", "LR", "LL") if feed_type.upper() == "RL" else ("XX", "XY", "YX", "YY")

    for ant in antennas:
        filename = "%s-%s.png" % (prefix, ant)
        print(("making plot %s" % filename))
        pylab.figure(figsize=(5 * ncols, 3 * 8))
        for i, src in enumerate(sorted(srcnames)):
            sols = DG[src]['solutions'][ant]
            for j, xx in enumerate(feeds):
                pylab.subplot(8, ncols, j * 2 * ncols + i + 1)
                pylab.title("%s:%s:%s:ampl" % (src, ant, xx))
                if not hasattr(sols[j], 'shape'):
                    continue
                ntime, nfreq = sols[j].shape
                pylab.xticks([])
                pylab.xlim(0, ntime - 1)
                if ylim:
                    pylab.ylim(*ylim)
                if nfreq > 1:
                    x = list(range(ntime))
                    pylab.fill_between(x, abs(sols[j][:, 0]), abs(sols[j][:, -1]), color='grey')
                    pylab.plot(abs(sols[j][:, nfreq / 2]))
                else:
                    pylab.plot(abs(sols[j][:, 0]))
                pylab.subplot(8, ncols, (j * 2 + 1) * ncols + i + 1)
                ph = numpy.angle(sols[j]) * 180 / math.pi
                if nfreq > 1:
                    pylab.plot(ph[:, 0], '.', ms=0.5, mec='0.2')
                    pylab.plot(ph[:, -1], '.', ms=0.5, mec='0.2')
                    pylab.plot(ph[:, nfreq / 2], '.b', ms=0.5)
                else:
                    pylab.plot(ph[:, 0], '.', ms=0.5)
                pylab.title("%s:%s:%s:phase (deg)" % (src, ant, xx))
                pylab.xticks([])
                pylab.xlim(0, ntime - 1)
        pylab.savefig(filename, dpi=150)
        pylab.close()

    ncols = len(srcnames)
    nrows = len(antennas)
    filename = "%s-ampl-summary.png" % prefix

    print(("making plot %s" % filename))
    pylab.figure(figsize=(5 * ncols, 3 * nrows))
    for row, ant in enumerate(antennas):
        for isrc, src in enumerate(sorted(srcnames)):
            sols = DG[src]['solutions'][ant][0]
            pylab.subplot(nrows, ncols, row * ncols + isrc + 1)
            pylab.title("%s:%s:%s:ampl" % (src, ant, "RR"))
            if not hasattr(sols, 'shape'):
                continue
            ntime, nfreq = sols.shape
            if nfreq > 1:
                x = list(range(ntime))
                pylab.fill_between(x, abs(sols[:, 0]), abs(sols[:, -1]), color='grey')
                pylab.plot(abs(sols[:, nfreq / 2]))
            else:
                pylab.plot(abs(sols[:, 0]))
            pylab.xticks([])
            pylab.xlim(0, ntime - 1)
            if ylim:
                pylab.ylim(*ylim)

    pylab.savefig(filename, dpi=75)
    pylab.close()


def make_ifrgain_plots(filename,
                       prefix=None, gain_label="IG",
                       feed=None,
                       feed_type="rl",
                       msname=None):
    """Makes a set of ifrgain plots from the specified saved file."""
    import pylab

    prefix = prefix or filename + "_ifrgain"

    # load baseline info, if MS is available
    baseline = None
    if msname:
        try:
            from pyrap.tables import table
            anttab = table(msname + "/ANTENNA")
            antpos = anttab.getcol("POSITION")
            antname = anttab.getcol("NAME")
            antpos = [(name, antpos[i]) for i, name in enumerate(antname)]
            # make dict of p,q to baseline length
            baseline = dict(
                [("%s-%s" % (p, q), math.sqrt(((ppos - qpos) ** 2).sum())) for p, ppos in antpos for q, qpos in antpos])
        except:
            error = "WARNING:: can't access %s antenna table, so IFR gain versus baseline length plots will not be made" % msname
            raise_with_traceback(ValueError(error))
    else:
        print(
            "WARNING:: MS not specified, thus no antenna info. IFR gain versus baseline length plots will not be made")

    # feed labels
    feeds = ("RR", "LL", "RL", "LR") if feed_type.upper() == "RL" else ("XX", "YY", "XY", "YX")

    ig = cPickle.load(open(filename))

    def plot_xy(content, title):
        """Plots x vs y"""
        pylab.errorbar(
            [x for l, (x, xe), (y, ye) in content], [y for l, (x, xe), (y, ye) in content],
            [ye for l, (x, xe), (y, ye) in content], [xe for l, (x, xe), (y, ye) in content],
            fmt=None, ecolor="lightgrey"
        )
        # plot labels
        for label, (x, xe), (y, ye) in content:
            pylab.text(x, y, label, horizontalalignment='center', verticalalignment='center', size=8)
        pylab.title(title)

    def plot_baseline(content, baseline, title, feeds):
        """Plots offset versus baseline"""
        bl = []
        xx = []
        xxe = []
        lab = []
        col = []
        for l, (x, xe), (y, ye) in content:
            b = baseline.get(l, None)
            if b is None:
                print(("WARNING:: baseline %s not found in MS ANTENNA table" % l))
            else:
                lab += ["%s:%s" % (l, feeds[0]), "%s:%s" % (l, feeds[1])]
                col += ["blue", "red"]
                xx += [x, y]
                xxe += [xe, ye]
                bl += [b, b]
        pylab.axhline(1, color='lightgrey')
        pylab.errorbar(bl, xx, yerr=xxe, fmt=None, ecolor="lightgrey")
        # plot labels
        for l, x, y, c in zip(lab, bl, xx, col):
            pylab.text(x, y, l, horizontalalignment='center', verticalalignment='center', size=8, color=c)
        pylab.xlabel("Baseline, m.")
        pylab.title(title)

    def plot_hist(content, title):
        """Plots histogram"""
        values = [x for l, (x, xe), (y, ye) in content] + [y for l, (x, xe), (y, ye) in content]
        hist = pylab.hist(values)
        pylab.xlim(min(values), max(values))
        pylab.title(title)

    def plot_complex(content, title):
        """Plots x vs y"""
        # plot errors bars, if available
        pylab.axhline(0, color='lightgrey')
        pylab.axvline(1, color='lightgrey')
        pylab.errorbar(
            [x.real for l1, l2, (x, xe), (y, ye) in content], [x.imag for l1, l2, (x, xe), (y, ye) in content],
            [xe for l1, l2, (x, xe), (y, ye) in content], [xe for l1, l2, (x, xe), (y, ye) in content],
            fmt=None, ecolor="lightgrey"
        )
        pylab.errorbar(
            [y.real for l1, l2, (x, xe), (y, ye) in content], [y.imag for l1, l2, (x, xe), (y, ye) in content],
            [ye for l1, l2, (x, xe), (y, ye) in content], [ye for l1, l2, (x, xe), (y, ye) in content],
            fmt=None, ecolor="lightgrey"
        )
        # max plot amplitude -- max point plus 1/4th of the error bar
        maxa = max([max(abs(x), abs(y)) for l1, l2, (x, xe), (y, ye) in content])
        # plotlim = max([ abs(numpy.array([
        #                  getattr(v,attr)+sign*e/4 for v,e in (x,xe),(y,ye) for attr in 'real','imag' for sign in 1,-1
        #                ])).max()
        #   for l1,l2,(x,xe),(y,ye) in content ])
        minre, maxre, minim, maxim = 2, -2, 2, -2
        for l1, l2, (x, xe), (y, ye) in content:
            offs = numpy.array([getattr(v, attr) + sign * e / 4 for v, e in [(x, xe), (y, ye)]
                                for attr in ['real', 'imag'] for sign in [1, -1]])
            minre, maxre = min(x.real - xe / 4, y.real - ye / 4, minre), max(x.real + xe / 4, y.real + ye / 4, maxre)
            minim, maxim = min(x.imag - xe / 4, y.imag - ye / 4, minim), max(x.imag + xe / 4, y.imag + ye / 4, maxim)
        # plot labels
        for l1, l2, (x, xe), (y, ye) in content:
            pylab.text(x.real, x.imag, l1, horizontalalignment='center', verticalalignment='center', color='blue',
                       size=8)
            pylab.text(y.real, y.imag, l2, horizontalalignment='center', verticalalignment='center', color='red',
                       size=8)
        # pylab.xlim(-plotlim,plotlim)
        # pylab.ylim(-plotlim,plotlim)
        pylab.xlim(minre, maxre)
        pylab.ylim(minim, maxim)
        pylab.title(title + " (max %.5g)" % maxa)

    def plot_ants(content, title):
        """Plots x vs y"""
        # plot labels
        for i, (p, gainlist) in enumerate(content):
            for label, color, (value, std) in gainlist:
                if value:
                    pylab.plot(i, value, 'w.')
                    pylab.text(i, value, label, horizontalalignment='center', verticalalignment='center', size=8,
                               color=color)
        pylab.xlim(-1, len(content))
        pylab.xticks(list(range(len(content))), [p for p, gainlist in content])
        pylab.title(title)

    antennas = set([p for p, q in list(ig.keys())])
    # plot RR vs LL offsets
    crl_diag = [("%s-%s" % (p, q), _normifrgain(rr), _normifrgain(ll))
                for (p, q), (rr, rl, lr, ll) in list(ig.items()) if not _is_unity(rr, ll)]
    have_diag = bool(crl_diag)
    # plot RL vs LR offsets, if present
    crl_offdiag = [("%s-%s" % (p, q), _normifrgain(rl), _normifrgain(lr))
                   for (p, q), (rr, rl, lr, ll) in list(ig.items()) if not _is_unity(rl, lr)]
    have_offdiag = bool(crl_offdiag)

    # plot size and layout
    FS = (16, 12) if not baseline else (16, 18)
    NR = 3 if baseline else 2
    NC = 2

    if have_diag:
        pylab.figure(figsize=FS)
        pylab.subplot(NR, NC, 3)
        plot_xy(crl_diag, "IFR gain amplitude (%s vs. %s)" % feeds[:2])
        pylab.subplot(NR, NC, 4)
        plot_hist(crl_diag, "IFR gain histogram for %s and %s" % feeds[:2])
        crl = [("%s-%s:%s" % (p, q, feeds[0].upper()), "%s-%s:%s" % (p, q, feeds[1].upper()), _complexifrgain(rr),
                _complexifrgain(ll))
               for (p, q), (rr, rl, lr, ll) in list(ig.items()) if not _is_unity(rr, ll)]
        pylab.subplot(NR, NC, 1)
        plot_complex(crl, "IFR complex %s %s gains" % feeds[:2])
        igpa = {}
        igpa0 = {}
        igpa0_means = []
        for (p, q), (rr, rl, lr, ll) in list(ig.items()):
            rr0 = abs(numpy.array(rr) - 1)
            ll0 = abs(numpy.array(ll) - 1)
            rr0 = numpy.ma.masked_array(rr0, rr0 == 0)
            ll0 = numpy.ma.masked_array(ll0, ll0 == 0)
            if rr0.count():
                igpa0_means += [rr0.mean()]
            if ll0.count():
                igpa0_means += [ll0.mean()]
            igpa0.setdefault(p, {})[q] = rr0, ll0
            igpa0.setdefault(q, {})[p] = rr0, ll0
            rr, ll = _normifrgain(rr), _normifrgain(ll)
            igpa.setdefault(p, []).append(("%s:%s" % (q, feeds[0]), 'blue', rr))
            igpa.setdefault(q, []).append(("%s:%s" % (p, feeds[0]), 'blue', rr))
            igpa.setdefault(p, []).append(("%s:%s" % (q, feeds[1]), 'red', ll))
            igpa.setdefault(q, []).append(("%s:%s" % (p, feeds[1]), 'red', ll))
        from past.builtins import cmp
        from functools import cmp_to_key
        content = [(p, igpa[p]) for p in sorted(list(igpa.keys()), key=cmp_to_key(_cmp_antenna))]
        pylab.subplot(NR, NC, 2)
        plot_ants(content, "IFR %s %s gain amplitudes per antenna" % feeds[:2])
        if baseline:
            pylab.subplot(NR, NC, 5)
            plot_baseline(crl_diag, baseline, "IFR gain amplitude vs. baseline length", feeds[:2])
        outname = "%s-%s.png" % (prefix, "".join(feeds[:2]))
        pylab.savefig(outname, dpi=75)
        print(("generated plot %s" % outname))
        # make per-antenna figure
        antennas = sorted(list(igpa0.keys()), key=cmp_to_key(_cmp_antenna))
        NC = 4
        NR = int(math.ceil(len(antennas) / float(NC)))
        offset = numpy.median(igpa0_means)
        pylab.figure(figsize=(8 * NC, 6 * NR))
        for iant, pant in enumerate(antennas):
            pylab.subplot(NR, NC, iant + 1)
            ig = igpa0[pant]
            ants1 = sorted(list(ig.keys()), key=cmp_to_key(_cmp_antenna))
            rr, ll = ig[ants1[0]]
            for i, qant in enumerate(ants1):
                rr, ll = ig[qant]
                if rr.count() > 1:
                    a1, a2 = numpy.ma.flatnotmasked_edges(rr)
                    baseline, = pylab.plot(rr + i * offset, '-')
                    pylab.text(a1, rr[a1] + i * offset, "%s:%s" % (qant, feeds[0]), horizontalalignment='left',
                               verticalalignment='center', size=8,
                               color=baseline.get_color())
                if ll.count() > 1:
                    a1, a2 = numpy.ma.flatnotmasked_edges(ll)
                    baseline, = pylab.plot(ll + i * offset, '-')
                    pylab.text(a2, ll[a2] + i * offset, "%s:%s" % (qant, feeds[1]), horizontalalignment='right',
                               verticalalignment='center', size=8,
                               color=baseline.get_color())
            pylab.title("antenna %s" % pant)

        outname = "%s-ant.png" % (prefix)
        pylab.savefig(outname, dpi=75)
        print(("generated plot %s" % outname))

    if have_offdiag:
        pylab.figure(figsize=FS)
        pylab.subplot(NR, NC, 3)
        plot_xy(crl_diag, "IFR offset amplitude (%s vs. %s)" % feeds[2:])
        pylab.subplot(NR, NC, 4)
        plot_hist(crl_diag, "IFR offset histogram for %s and %s" % feeds[2:])
        crl = [("%s-%s:%s" % (p, q, feeds[0].upper()), "%s-%s:%s" % (p, q, feeds[1].upper()), _complexifrgain(rl),
                _complexifrgain(lr))
               for (p, q), (rr, rl, lr, ll) in list(ig.items()) if not _is_unity(rl, lr)]
        pylab.subplot(NR, NC, 1)
        plot_complex(crl_diag, "IFR complex %s %s offsets" % feeds[2:])
        igpa = {}
        for (p, q), (rr, rl, lr, ll) in list(ig.items()):
            rr, ll = _normifrgain(rl), _normifrgain(lr)
            igpa.setdefault(p, []).append(("%s:%s" % (q, feeds[2]), 'blue', rr))
            igpa.setdefault(q, []).append(("%s:%s" % (p, feeds[2]), 'blue', rr))
            igpa.setdefault(p, []).append(("%s:%s" % (q, feeds[3]), 'red', ll))
            igpa.setdefault(q, []).append(("%s:%s" % (p, feeds[3]), 'red', ll))
        from past.builtins import cmp
        from functools import cmp_to_key
        content = [(p, igpa[p]) for p in sorted(list(igpa.keys()), key=cmp_to_key(_cmp_antenna))]
        pylab.subplot(NR, NC, 2)
        plot_ants(content, "IFR %s %s offset amplitudes per antenna" % feeds[2:])
        if baseline:
            pylab.subplot(NR, NC, 5)
            plot_baseline(crl_offdiag, baseline, "IFR offset amplitude by baseline length", feeds[:2])

        outname = "%s.png" % (prefix)
        pylab.savefig(outname, dpi=75)
        print(("generated plot %s" % outname))
