# -*- coding: utf-8 -*-
import numpy as np
import re
import tempfile
import os
from past.builtins import cmp
from functools import cmp_to_key

import Owlcat
from Kittens.utils import verbosity
try:
    import Purr.Pipe

    has_purr = True
except:
    has_purr = False

_gli = None # Meow.MSUtils.find_exec('glish')


# Various argument-formatting methods to use with the Flagger.AutoFlagger class
# below. These really should be static methods of the class, but that doesn't work
# with Python (specifically, I cannot include them into static member dicts)
def _format_nbins(bins, argname):
    if isinstance(bins, (list, tuple)) and len(bins) == 2:
        return str(list(bins))
    else:
        return "%d" % bins
    raise TypeError("invalid value for '%s' keyword (%s)" % (argname, bins))


def _format_bool(arg, argname):
    if arg:
        return "T"
    else:
        return "F"


def _format_index(arg, argname):
    if isinstance(arg, int):
        if arg < 0:
            return str(arg)
        else:
            return str(arg + 1)
    raise TypeError("invalid value for '%s' keyword (%s)" % (argname, arg))


def _format_list(arg, argname):
    if isinstance(arg, str):
        return "'%s'" % arg
    elif isinstance(arg, (list, tuple)):
        return "[%s]" % (','.join([_format_list(x, argname) for x in arg]))
    elif isinstance(arg, bool):
        return _format_bool(arg, argname)
    elif isinstance(arg, (int, float)):
        return str(arg)
    raise TypeError("invalid value for '%s' keyword (%s)" % (argname, arg))


def _format_ilist(arg, argname):
    if isinstance(arg, str):
        return "'%s'" % arg
    elif isinstance(arg, (list, tuple)):
        return "[%s]" % (','.join([_format_ilist(x, argname) for x in arg]))
    elif isinstance(arg, bool):
        return _format_bool(arg, argname)
    elif isinstance(arg, int):
        return _format_index(arg, argname)
    elif isinstance(arg, float):
        return str(arg)
    raise TypeError("invalid value for '%s' keyword (%s)" % (argname, arg))


def _format_plotchan(arg, argname):
    if isinstance(arg, bool):
        return _format_bool(arg, argname)
    else:
        return _format_index(arg, argname)


def _format_2N(arg, argname):
    if isinstance(arg, (list, tuple)):
        return "%s" % arg
    elif numpy.isarray(arg):
        if arg.dtype is np.int32:
            a = arg + 1
        else:
            a = arg
        if a.ndim == 1:
            return str(list(a))
        elif a.ndim == 2 and a.shape[1] == 2:
            return "[%s]" % (','.join([str(list(a[i, :])) for i in range(a.shape[0])]))
    raise TypeError("invalid value for '%s' keyword (%s)" % (argname, arg))


def _format_clip(arg, argname):
    if isinstance(arg, (list, tuple)):
        if not isinstance(arg[0], (list, tuple)):
            arg = [arg]
        recfields = []
        for i, triplet in enumerate(arg):
            if not isinstance(triplet, (list, tuple)) or len(triplet) != 3:
                raise TypeError("invalid value for '%s' keyword (%s)" % (argname, arg))
            expr, minval, maxval = triplet
            subfields = ["expr='%s'" % expr]
            if minval is not None:
                subfields.append("min=%g" % minval)
            if maxval is not None:
                subfields.append("max=%g" % maxval)
            recfields.append("a%d=[%s]" % (i, ','.join(subfields)))
        return "[%s]" % ','.join(recfields)
    raise TypeError("invalid value for '%s' keyword (%s)" % (argname, arg))


# mapping from casacore table valueType to bit length
BITFLAG_NBITS = dict(uchar=8, short=16, int=32)

# mapping from casacore table valueType to numpy dtype
BITFLAG_DTYPES = dict(uchar=np.uint8, short=np.int16, int=np.int32)



class Flagsets(object):
    """ Manage an measurement set's flagsets. Pasted from Cattery.Meow.MSUtils. """

    def __init__(self, ms):
        """
        Initialises a Flagsets object for the measurement set.

        Args:
            ms (:obj:`~casacore.tables.table.table`):
                A table object belonging to the measurement set.
        """

        self.ms = ms
        if not 'BITFLAG' in ms.colnames():
            self.order = None
            self.bits = {}
            self.NBITS = 8
            self.dtype = np.uint8
        else:
            value_type = ms.getcoldesc('BITFLAG')['valueType']
            self.NBITS = BITFLAG_NBITS[value_type]
            self.dtype = BITFLAG_DTYPES[value_type]
            kws = ms.colkeywordnames('BITFLAG')
            self.bits = {}
            # scan FLAGSET_xxx keywords and populate name->bitmask mappings
            for kw in kws:
                match = re.match('^FLAGSET_(.*)$', kw)
                if match:
                    name = match.group(1)
                    bit = ms.getcolkeyword('BITFLAG', kw)
                    if isinstance(bit, int):
                        self.bits[name] = bit
                    else:
                        print("Warning: unexpected type (%s) for %s keyword of BITFLAG column," \
                              " ignoring" % (type(bit), kw))
            # have we found any FLAGSET_ specs?
            if self.bits:
                order = 'FLAGSETS' in kws and ms.getcolkeyword('BITFLAG', 'FLAGSETS')
                if isinstance(order, str):
                    order = order.split(',')
                else:
                    print("Warning: unexpected type (%s) for FLAGSETS keyword of BITFLAG column," \
                          " ignoring" % type(order))
                    order = []
                # form up "natural" order by comparing bitmasks
                bitwise_order = list(self.bits.keys())
                bitwise_order.sort(key=cmp_to_key(lambda a, b: cmp(self.bits[a], self.bits[b])))
                # if an order is specified, make sure it is actually valid,
                # and add any elements from bitwise_order that are not present
                self.order = [fs for fs in order if fs in self.bits] + \
                             [fs for fs in bitwise_order if fs not in order]
                # if order was fixed compared to what was in MS, write back to MS
                if ms.iswritable() and self.order != order:
                    ms._putkeyword('BITFLAG', 'FLAGSETS', -1, False, ','.join(self.order))
                    ms.flush()
            # else if no flagsets found, try the old-style NAMES keyword
            elif 'NAMES' in kws:
                names = ms.getcolkeyword('BITFLAG', 'NAMES')
                if isinstance(names, (list, tuple)):
                    self.order = list(map(str, names))
                    bit = 1
                    for name in self.order:
                        self.bits[name] = bit
                        bit <<= 1
                    if ms.iswritable():
                        ms._putkeyword('BITFLAG', 'FLAGSETS', -1, False, ','.join(self.order))
                        for name, bit in self.bits.items():
                            ms._putkeyword('BITFLAG', 'FLAGSET_%s' % name, -1, False, bit)
                        ms.flush()
            else:
                self.order = []
        self.LEGACY = 1 << self.NBITS
        self.BITMASK_ALL = self.LEGACY - 1

    def names(self):
        """
        Convenience function for determining active flagsets.

        Returns:
            list or None:
                A list of flagset names, in the order in which they were created or None if BITFLAG
                column is missing (so flagsets are unavailable.)
        """

        return self.order

    def flagmask(self, name, create=False):
        """
        Flagmask corresponding to named flagset.

        Args:
            name (str):
                Name of flagset.
            create (bool, optional):
                If True and flagset is missing, creates named flagset, else raises exception.

        Raises:
            TypeError:
                If the MS does not contain a BITFLAG column.
            ValueError:
                If the named flagset is not found and create is False.
            ValueError:
                If there are too many flagsets to create a new one.
        """

        # Cludge for python2/3 interoperability.
        name = str(name)

        # lookup flagbit, return if found
        if self.order is None:
            raise TypeError("MS does not contain a BITFLAG column. Please run the addbitflagcol" \
                            " utility on this MS.")
        bit = self.bits.get(name, None)
        if bit is not None:
            return bit
        # raise exception if not allowed to create a new one
        if not create:
            raise ValueError("Flagset '%s' not found" % name)
        # find empty bit
        for bitnum in range(self.NBITS):
            bit = 1 << bitnum
            if bit not in list(self.bits.values()):
                self.order.append(name)
                self.bits[name] = bit
                self.ms._putkeyword('BITFLAG', 'FLAGSETS', -1, False, ','.join(self.order))
                self.ms._putkeyword('BITFLAG', 'FLAGSET_%s' % name, -1, False, bit)
                self.ms.flush()
                return bit
        # no free bit found, bummer
        raise ValueError("Too many flagsets in MS, cannot create another one")

    def remove_flagset(self, *fsnames):
        """
        Removes the named flagset(s).

        Args:
            fsnames (tuple):
                Names of flagsets to be removed.

        Returns:
            int:
                Flagmask corresponding to the removed flagsets.
        """

        # lookup all flagsets, raise error if any not found
        if self.bits is None:
            raise TypeError("MS does not contain a BITFLAG column, cannot use flagsets")
        removing = []
        for fs in fsnames:
            bit = self.bits.get(fs, None)
            if bit is None:
                raise ValueError("Flagset '%s' not found" % fs)
            removing.append((fs, bit))
        if not removing:
            return
        # remove items, form up mask of bitflags to be cleared
        mask = 0
        for name, bit in removing:
            mask |= bit
            del self.bits[name]
            del self.order[self.order.index(name)]
            self.ms.removecolkeyword('BITFLAG', 'FLAGSET_%s' % name)
        # write new list of bitflags
        self.ms._putkeyword('BITFLAG', 'FLAGSETS', -1, False, ','.join(self.order))
        self.ms.flush()

        return mask


class Flagger(verbosity):
    def __init__(self, msname, verbose=0, timestamps=False, chunksize=200000):
        verbosity.__init__(self, name="Flagger")
        self.set_verbose(verbose)
        if timestamps:
            self.enable_timestamps(modulo=10000)
        self.msname = msname
        self.ms = None
        self.readwrite = False
        self.chunksize = chunksize
        self._initializing_bitflags = False  # will be true if initializing
        self._reopen()

    def close(self):
        if self.ms:
            self.dprint(2, "closing MS", self.msname)
            self.ms.close()
        self.ms = self.flagsets = None

    def _reopen(self, readwrite=False):
        if self.ms is None:
            self.ms = ms = Owlcat.table(self.msname, readonly=not readwrite, ack=False)
            self.readwrite = readwrite
            self.dprintf(1, "opened MS %s, %d rows\n", ms.name(), ms.nrows())
            has_bf = 'BITFLAG' in ms.colnames()
            has_bfr = 'BITFLAG_ROW' in ms.colnames()
            self.has_bitflags = has_bf and has_bfr
            if self.has_bitflags:
                self.has_bitflags = self.ms.getcoldesc("BITFLAG")['valueType']
                self.dprint(1, "this MS has proper bitflag columns of type '%s'" % self.has_bitflags)
            else:
                if has_bf or has_bfr:
                    self.dprint(0, "WARNING: bitflag columns only partially initialized. Try --reinit-bitflags?")
            self.flagsets = Flagsets(ms)
            self.dprintf(1, "flagsets are %s\n", self.flagsets.names())
            self.purrpipe = Purr.Pipe.Pipe(self.msname) if has_purr else None
            self.LEGACY = self.flagsets.LEGACY
            self.BITMASK_ALL = self.flagsets.BITMASK_ALL

        elif self.readwrite != readwrite:
            self.ms = Owlcat.table(self.msname, readonly=not readwrite, ack=False)
            self.readwrite = readwrite

        return self.ms

    def _add_column(self, col_name, like_col="DATA", like_type=None, stman=None):
        """
        Inserts a new column into the measurement set.

        Args:
            col_name (str):
                Name of target column.
            like_col (str, optional):
                Column will be patterned on the named column.
            like_type (str or None, optional):
                If set, column type will be changed.

        Returns:
            bool:
                True if a new column was inserted, else False.
        """

        if col_name in self.ms.colnames():
            return False
        self._reopen(True)
        # new column needs to be inserted -- get column description from column 'like_col'
        desc = self.ms.getcoldesc(like_col)
        desc['name'] = col_name
        desc['comment'] = desc['comment'].replace(" ", "_")  # casacore not fond of spaces...
        dminfo = self.ms.getdminfo(like_col)
        dminfo["NAME"] =  "{}-{}".format(dminfo["NAME"], col_name)
        if stman is not None:
            dminfo["TYPE"] = stman
        # if a different type is specified, insert that
        if like_type:
            desc['valueType'] = like_type
        self.dprint(0, "inserting new column %s of type '%s' based on %s" % (col_name, desc['valueType'], like_col))
        if stman is not None:
            self.dprint(0, "  using %s" % stman)
        self.ms.addcols(desc, dminfo)
        return True

    def remove_bitflags(self):
        removecols = [col for col in ("BITFLAG", "BITFLAG_ROW") if col in self.ms.colnames()]
        if removecols:
            self._reopen(True)
            self.dprint(0, "removing columns", ", ".join(removecols))
            self.ms.removecols(removecols)
            self.close()
            self._reopen()

    def add_bitflags(self, wait=True, purr=True, bits=32, stman=None):
        if not self.has_bitflags:
            self._reopen(True)
            if bits == 32:
                dtype = 'int'
                self.has_bitflags = np.int32
            elif bits == 16:
                dtype = 'short'
                self.has_bitflags = np.int16
            elif bits == 8:
                dtype = 'uchar'
                self.has_bitflags = np.uint16
            else:
                raise ValueError("invalid bits setting. 8, 16 or 32 expected.")
            self._add_column("BITFLAG", "DATA", dtype, stman=stman)
            self._add_column("BITFLAG_ROW", "FLAG_ROW", dtype, stman=stman)
            self.close()
            self._reopen()
            self._initializing_bitflags = True

    def remove_flagset(self, *fsnames, **kw):
        self._reopen(True)
        mask = self.flagsets.remove_flagset(*fsnames)
        # report to purr
        if kw.get('purr', True) and self.purrpipe:
            self.purrpipe.title("Flagging").comment("Removing flagset(s) %s." % (','.join(fsnames)))
        self.unflag(mask, purr=False)


    def _get_submss(self, ms, ddids=None):
        """Helper method. Splits MS into subsets by DATA_DESC_ID.
        Returns list of (ddid,nrows,subms) tuples, where subms is a subset of the MS with the given DDID,
        and nrows is the number of rows in preceding submss (which is needed for progress stats)
        """
        # build up list of (ddid,nrows,subms) pairs.
        # subms is a subset of the MS with the given DDID
        # nrows is the number of rows in preceding submss (need for progress stats)
        sub_mss = []
        nrows = 0
        if ddids is None:
            ddids = list(range(TABLE(ms.getkeyword('DATA_DESCRIPTION'), ack=False).nrows()))
        for ddid in ddids:
            subms = ms.query("DATA_DESC_ID==%d" % ddid)
            sub_mss.append((ddid, nrows, subms))
            nrows += subms.nrows()
        return sub_mss


    def lookup_flagmask(self, flagset, create=False):
        """helper function: converts a flagset name into an integer flagmask"""
        if flagset is None:
            return None
        # convert flagset to list
        elif isinstance(flagset, int):
            return int(flagset)
        elif isinstance(flagset, str):
            flagset = flagset.split(",")
        elif not isinstance(flagset, (list, tuple)):
            raise TypeError("invalid flagset of type %s" % type(flagset))
        # loop over list and accumulate flagmask
        flagmask = 0
        for fset in flagset:
            if isinstance(fset, int):
                flagmask |= fset
            elif isinstance(fset, str):
                if fset == "ALL":
                    flagmask |= self.BITMASK_ALL | self.LEGACY
                elif fset.upper() == "ALL":
                    flagmask |= self.BITMASK_ALL
                else:
                    # +L suffix includes legacy flags
                    if fset[-2:].upper() == "+L":
                        flagmask |= self.LEGACY
                        fset = fset[:-2]
                    # lookup name or number
                    if re.match('^\d+$', fset):
                        flagmask |= int(fset)
                    elif re.match('^0[xX][\dA-Fa-f]+$', fset):
                        flagmask |= int(fset, 16)
                    elif fset:
                        flagmask |= self.flagsets.flagmask(fset, create=create)
            else:
                raise TypeError("invalid flagset of type %s in list of flagsets" % type(fset))
        return flagmask

    def flagmaskstr(self, flagmask):
        """helper function: converts an integer flagmask into a printable str"""
        if not flagmask:
            return "0"
        legacy = flagmask & self.LEGACY
        bits = flagmask & self.BITMASK_ALL
        return "0x%04X%s" % (bits, "+L" if legacy else "")

    def xflag(self,
              # These options determine what to do with the subset determined below. If none are set,
              # the xflag() simply returns stats without doing anything
              flag=None,  # set the flagmask (or flagset name, optionally "+L")
              unflag=None,  # clear the flagmask (or flagset name, optionally "+L")
              copy=None,    # copy to flagmask (or flagset name, optionally "+L")
              create=False,  # if True and 'flag' is a string, creates new flagset as needed
              fill_legacy=None,  # if not None, legacy flags will be filled (within all of subset A) using this mask

              # The following options restrict the subset to be flagged. Effect is cumulative.
              # Row subset. Data selection by whole rows
              ddid=None, fieldid=None,  # DATA_DESC_ID or FIELD_ID, single int or list
              antennas=None,  # list of antennas
              baselines=None,  # list of baselines (as ip,iq pairs)
              time=None,  # time range as (min,max) pair (can use None for either min or max)
              reltime=None,  # relative time (rel. to start of MS) range as (min,max) pair
              taql=None,  # additional TaQL string
              # Subset A. Freq/corr slices within the row subset.
              channels=None,  # channel subset (index, or slice, or list of index/slices)
              corrs=None,  # correlation subset (index, or slice, or list of index/slices)
              # Subset B. Subset within subset A based on flags
              flagmask=None,  # any bitflag set in the given flagmask (or flagset name)
              flagmask_all=None,  # all bitflags set in the given flagmask (or flagset name)
              flagmask_none=None,  # no bitflags set in the given flagmask (or flagset name)
              # Subset C. Subset within subset B based on data clipping
              data_nan=False,  # restrict flagged subset to NAN/INF data
              data_above=None,  # restrict flagged subset to abs(data)>X
              data_below=None,  # and abs(data)<X
              data_fm_above=None,  # same as data_above/_below, but flags based on the mean
              data_fm_below=None,  # amplitude across all frequencies
              data_column='CORRECTED_DATA',  # data column for clip_above and clip_below
              data_flagmask=-1,  # flagmask to apply to data column when computing mean
              # if True, dict of {name: mask} bitflags for which statistics will be collected
              get_stats_only=None,
              # other options
              flag_allcorr=True,  # flag all correlations if at least one is flagged
              progress_callback=None,  # callback, called with (n,nmax) to report progress
              purr=False  # if True, writes comments to purrpipe
              ):
        """Alternative flag interface, works on the in/out principle."""
        if not self.purrpipe:
            purr = False
        ms = self._reopen(flag or unflag or copy or fill_legacy is not None)
        # lookup flagset names
        for var in 'flag', 'unflag', 'copy', 'flagmask', 'flagmask_all', 'flagmask_none', 'data_flagmask', 'fill_legacy':
            flagval = locals()[var]
            locals()[var] = fm = self.lookup_flagmask(flagval, create=(var == 'flag'))
            if flagval is not None:
                self.dprintf(2, "%s=%s corresponds to bitmask %s\n", var, flagval, self.flagmaskstr(fm))
        # for these two masks, it's more convenient that they're set to 0 if missing
        flag = flag or 0
        unflag = unflag or 0
        copy = copy or 0
        if not self.has_bitflags and (flag | unflag | copy) & self.BITMASK_ALL:
            raise RuntimeError("no BITFLAG column in this MS, can't change bitflags")
        # stats accumulator
        if get_stats_only:
            stats = { name: 0 for name in get_stats_only.keys() }
        # total number of rows
        totrows = ms.nrows()
        # rows and visibilities in the row subset
        sel_nrow = sel_nvis = 0
        # visibilities selected in subsets A, B and C
        nvis_A = nvis_B = nvis_C = 0
        # get DDIDs and FIELD_IDs
        if ddid is None:
            ddids = list(range(Owlcat.table(ms.getkeyword('DATA_DESCRIPTION'), ack=False, readonly=True).nrows()))
        elif isinstance(ddid, int):
            ddids = [ddid]
        elif isinstance(ddid, (tuple, list)):
            ddids = ddid
        else:
            raise TypeError("invalid ddid argument of type %s" % type(ddid))
        # form up list of TaQL expressions for subset selectors
        queries = []
        if taql:
            queries.append(taql)
        if fieldid is not None:
            if isinstance(fieldid, int):
                fieldid = [fieldid]
            elif not isinstance(fieldid, (tuple, list)):
                raise TypeError("invalid fieldid argument of type %s" % type(fieldid))
            queries.append(" || ".join(["FIELD_ID==%d" % f for f in fieldid]))
        if antennas is not None:
            antlist = str(list(antennas))
            queries.append("ANTENNA1 in %s || ANTENNA2 in %s" % (antlist, antlist))
        if time is not None:
            t0, t1 = time
            if t0 is not None:
                queries.append("TIME>=%g" % t0)
            if t1 is not None:
                queries.append("TIME<=%g" % t1)
        if reltime is not None:
            t0, t1 = reltime
            time0 = self.ms.getcol('TIME', 0, 1)[0]
            if t0 is not None:
                queries.append("TIME>=%f" % (time0 + t0))
            if t1 is not None:
                queries.append("TIME<=%f" % (time0 + t1))
        # form up TaQL string, and extract subset of table
        if queries:
            query = "( " + " ) && ( ".join(queries) + " )"
            purr and self.purrpipe.comment("; effective MS selection is \"%s\"" % query, endline=False)
            self.dprintf(2, "selection string is %s\n", query)
            ms = ms.query(query)
            self.dprintf(2, "query reduces MS to %d rows\n", ms.nrows())
        else:
            self.dprintf(2, "no selection applied\n")
        # check list of baselines
        if baselines:
            baselines = [(int(p), int(q)) for p, q in baselines]
            purr and self.purrpipe.comment("; baseline subset is %s" %
                                           " ".join(["%d-%d" % (p, q) for p, q in baselines]),
                                           endline=False)

        # helper func to parse the channels/corrs/timeslots arguments
        def make_slice_list(selection, parm):
            if not selection:
                return [np.s_[:]]
            if isinstance(selection, (int, slice)):
                return make_slice_list([selection], parm)
            if not isinstance(selection, (list, tuple)):
                raise TypeError("invalid %s selection: %s" % (parm, selection))
            sellist = []
            for sel in selection:
                if isinstance(sel, int):
                    sellist.append(slice(sel, sel + 1))
                elif isinstance(sel, slice):
                    sellist.append(sel)
                else:
                    raise TypeError("invalid %s selection: %s" % (parm, selection))
            return sellist

        # parse the arguments
        channels = make_slice_list(channels, 'channels')
        corrs = make_slice_list(corrs, 'corrs')
        purr and self.purrpipe.comment("; channels are %s" % channels, endline=False)
        purr and self.purrpipe.comment("; correlations are %s" % corrs, endline=False)
        # put comment into purrpipe
        purr and self.purrpipe.comment(".")
        #
        flagsubsets = flagmask is not None or flagmask_all is not None or flagmask_none is not None
        dataclip = data_above is not None or data_below is not None or data_nan or \
                   data_fm_above is not None or data_fm_below is not None
        # make list of sub-MSs by DDID
        sub_mss = self._get_submss(ms, ddids)
        nrow_tot = ms.nrows()
        # go through rows of the MS in chunks
        for ddid, irow_prev, ms in sub_mss:
            self.dprintf(2, "processing MS subset for ddid %d\n", ddid)
            if progress_callback:
                progress_callback(irow_prev, nrow_tot)
            for row0 in range(0, ms.nrows(), self.chunksize):
                if progress_callback:
                    progress_callback(irow_prev + row0, nrow_tot)
                nrows = min(self.chunksize, ms.nrows() - row0)
                self.dprintf(2, "processing rows %d:%d (%d rows total)\n", row0, row0 + nrows - 1, nrows)
                # apply baseline selection to the mask
                if baselines:
                    # rowmask will be True for all selected rows
                    rowmask = np.zeros(nrows, bool)
                    a1 = ms.getcol('ANTENNA1', row0, nrows)
                    a2 = ms.getcol('ANTENNA2', row0, nrows)
                    # update mask
                    for p, q in baselines:
                        rowmask |= (a1 == p) & (a2 == q)
                    self.dprintf(2, "baseline selection leaves %d rows\n", nr)
                # else select all rows
                else:
                    # rowmask will be True for all selected rows
                    rowmask = np.ones(nrows, bool)
                # read legacy flags to get a datashape
                legacy_flag_column = ms.getcol('FLAG', row0, nrows)
                legacy_flag_column |= ms.getcol('FLAG_ROW', row0, nrows)[:, np.newaxis, np.newaxis]
                datashape = legacy_flag_column.shape
                nv_per_row = datashape[1] * datashape[2]
                # rowflags and visflags will be constructed on-demand below. Make helper functions for this
                self._visflags = None

                def bitflags():
                    """Reads BITFLAG column on demand. If no such column, forms up zero array"""
                    if self._visflags is None:
                        if not self.has_bitflags or self._initializing_bitflags:
                            self._visflags = np.zeros_like(legacy_flag_column, self.flagsets.dtype)
                            self.dprint(2, "formed empty BITFLAGs (type {})".format(self.flagsets.dtype))
                        elif all([ms.iscelldefined('BITFLAG', r) for r in range(row0, row0+nrows)]):
                            self._visflags = ms.getcol('BITFLAG', row0, nrows)
                            self.dprint(2, "read BITFLAG column")
                            self._visflags |= ms.getcol('BITFLAG_ROW', row0, nrows)[:, np.newaxis, np.newaxis]
                            self.bf_type = self._visflags.dtype.type
                        else:
                            self._visflags = np.zeros_like(legacy_flag_column, self.flagsets.dtype)
                            self.dprint(2, "formed empty BITFLAGs (type {}) because part (or all) of the column is undefined".format(self.flagsets.dtype))
                        self._bf_type = self._visflags.dtype.type
                    return self._visflags

                # apply stats
                nr = rowmask.sum()
                sel_nrow += nr
                nv = nr * nv_per_row
                sel_nvis += nv
                self.dprintf(2, "Row subset (data selection) leaves %d rows and %d visibilities\n", nr, nv)
                # get subset C
                # vismask will be True for all selected visibilities
                vismask = np.zeros(datashape, bool)
                for channel_slice in channels:
                    for corr_slice in corrs:
                        vismask[rowmask, channel_slice, corr_slice] = True
                self.dprintf(3, "  finished applying freq/corr slicing\n")
                nv = vismask.sum()
                nvis_A += nv                
                self.dprintf(2, "subset A (freq/corr slicing) leaves %d visibilities\n", nv)
                # stats-only mode
                if get_stats_only:
                    for name, mask in get_stats_only.items():
                        if mask is None:
                            nv = np.count_nonzero(legacy_flag_column)
                        else:
                            vf1 = bitflags()&mask
                            nv = np.count_nonzero(vf1)
                        stats[name] += nv
                        self.dprintf(3, "flagset %s flags %d visibilities\n", name, nv)
                    continue  # to top of loop to get next chunk

                # select based on flag state
                if flagsubsets:
                    # select using the "any" scheme
                    if flagmask is not None:
                        if flagmask == self.LEGACY:
                            vismask &= legacy_flag_column
                        else:
                            bf1 = (bitflags()&flagmask) != 0
                            if flagmask&self.LEGACY:
                                bf1 |= legacy_flag_column
                            vismask &= bf1
                    # select using the "all" scheme
                    if flagmask_all is not None:
                        if flagmask_all == self.LEGACY:
                            vismask &= legacy_flag_column
                        else:
                            fm = flagmask_all & ~self.LEGACY
                            bf1 = (bitflags()&fm) == fm
                            if flagmask_all&self.LEGACY:
                                bf1 &= legacy_flag_column
                            vismask &= bf1
                    # select using the "None" scheme
                    if flagmask_none is not None:
                        if flagmask_none == self.LEGACY:
                            vismask &= ~legacy_flag_column
                        else:
                            fm = flagmask_none & ~self.LEGACY
                            bf1 = (bitflags()&flagmask_none) == 0
                            if flagmask_none&self.LEGACY:
                                bf1 &= ~legacy_flag_column
                            vismask &= bf1
                nv = vismask.sum()
                nvis_B += nv
                self.dprintf(2, "subset B (flag-based selection) leaves %d visibilities\n", nv)
                # now apply clipping
                if dataclip:
                    datacol = ms.getcol(data_column, row0, nrows)
                    # make it a masked array: mask out stuff not in vismask
                    datamask = ~vismask
                    # and mask stuff in data_flagmask
                    if data_flagmask is not None:
                        datamask |= ((bitflags() & data_flagmask) != 0)
                    datacol = np.ma.masked_array(datacol, datamask)
                    self.dprintf(4, "datamask contains %d masked visibilities\n", datamask.sum())
                    self.dprintf(3, "At start of clipping we have %d visibilities\n", vismask.sum())
                    self.dprintf(3, "of which %d are finite\n", np.isfinite(datacol).sum())
                    # clip on NANs
                    if data_nan:
                        vismask &= ~np.isfinite(datacol)
                        self.dprintf(3, "NAN filtering leaves %d visibilities\n", vismask.sum())
                    # clip on amplitudes
                    abscol = abs(datacol)
                    if data_above is not None:
                        vismask &= abscol > data_above
                        self.dprintf(3, "data_above filtering leaves %d visibilities\n", vismask.sum())
                    if data_below is not None:
                        vismask &= abscol < data_below
                        self.dprintf(3, "data_below filtering leaves %d visibilities\n", vismask.sum())
                    # clip on freq-mean amplitudes
                    if data_fm_above is not None or data_fm_below is not None:
                        datacol = abscol.mean(1)
                        #            print datacol.max(),data_fm_above
                        if data_fm_above is not None:
                            vismask &= (datacol > data_fm_above)[:, np.newaxis, ...]
                        if data_fm_below is not None:
                            vismask &= (datacol < data_fm_below)[:, np.newaxis, ...]
                        self.dprintf(3, "data_fm_above/below filtering leaves %d visibilities\n", vismask.sum())
                # finally, subset E is ready
                nv = vismask.sum()
                self.dprintf(2, "subset C (data clipping) leaves %d visibilities\n", nv)
                # extending flagging to all correlations
                if flag_allcorr:
                    vismask |= np.logical_or.reduce(vismask, 2)[:, :, np.newaxis]
                    nv = vismask.sum()
                    self.dprintf(2, "which extends to %d visibilities with flag_allcorr in effect\n", nv)
                nvis_C += nv

                # ok, at this point vismask is a boolean array of things to flag or unflag.
                # BITFLAGS would only have been read in if necessary

                # now, do the actual flagging
                if flag or unflag or copy or fill_legacy is not None:
                    changing_legacy_flags = fill_legacy is not None or (flag|unflag|copy)&self.LEGACY
                    # flag/unflag visibilities
                    if flag:
                        self.dprint(4, "doing flag")
                        if flag&self.LEGACY:
                            legacy_flag_column[vismask] = True
                        if flag&~self.LEGACY:
                            bitflags()[vismask] |= self._bf_type(flag&~self.LEGACY)
                    if copy:
                        self.dprint(4, "doing flag-copy")
                        if copy&self.LEGACY:
                            legacy_flag_column = vismask
                        if copy&~self.LEGACY:
                            bf = bitflags()
                            bf[vismask] |= self._bf_type(copy&~self.LEGACY)
                    if unflag:
                        self.dprint(4, "doing unflag")
                        if unflag&self.LEGACY:
                            legacy_flag_column[unflag] = False
                        if unflag&~self.LEGACY:
                            bitflags()[vismask] &= ~self._bf_type(unflag&~self.LEGACY)
                    # fill legacy flags
                    if fill_legacy is not None:
                        self.dprint(4, "filling legacy flags")
                        bf = bitflags()
                        legacy_flag_column[rowmask,...] = (bf[rowmask,...] & self._bf_type(fill_legacy)) != 0
                    # adjust the rowflags
                    self.dprint(4, "adjusting rowflags")

                    # were any bitflag manipulations done -- write bitflags
                    if self.has_bitflags and (flag | unflag | copy) & ~self.LEGACY:
                        bf_row = np.bitwise_and.reduce(bitflags(), axis=(1,2))
                        ms.putcol('BITFLAG', bitflags(), row0, nrows)
                        ms.putcol('BITFLAG_ROW', bf_row, row0, nrows)
                        self.dprint(4, "wrote BITFLAG")
                    # write legacy flags
                    if changing_legacy_flags:
                        self.dprintf(4, "filling legacy flags for rows %d:%d\n" % (row0, row0 + nrows))
                        ms.putcol('FLAG', legacy_flag_column, row0, nrows)
                        ms.putcol('FLAG_ROW', legacy_flag_column.all(axis=(1, 2)), row0, nrows)
                        self.dprint(4, "wrote FLAG")
                self.dprint(4, "done with this chunk")
        if progress_callback:
            progress_callback(99, 100)
        self._rowflags = self._visflags = None
        # print collected stats
        self.dprint(1, "xflag stats:")
        self.dprintf(1, "total MS size:           %8d rows\n", totrows)
        self.dprintf(1, "data selection leaves    %8d rows, %8d visibilities\n", sel_nrow, sel_nvis)
        self.dprintf(1, "chan/corr slicing leaves %8s       %8d visibilities\n", '', nvis_A)
        self.dprintf(1, "flag selection leaves    %8s       %8d visibilities\n", '', nvis_B)
        self.dprintf(1, "data clipping leaves     %8s       %8d visibilities\n", '', nvis_C)

        if get_stats_only:
            for name, nv in stats.items():
                get_stats_only[name] = nv

        return totrows, sel_nrow, sel_nvis, nvis_A, nvis_B, nvis_C

