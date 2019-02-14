# -*- coding: utf-8 -*-

from six.moves import cPickle
import gzip


class TableDump(object):
    # initializes a table dump file
    # format:
    #   header dict with version=2,
    #   followed by colname,keywords,nchunks(int),chunk1,...,chunkN for each column
    def __init__(self, filename, version=1, compress=True, write=False):
        if compress:
            self.stream = gzip.GzipFile(filename, "w" if write else "r")
        else:
            self.stream = open(filename, "wb" if write else "rb")
        if write:
            cPickle.dump(dict(format='Owlcat MS column dump', version=version), self.stream)
        else:
            self.header = cPickle.load(self.stream)

    def __del__(self):
        self.close()

    def close(self):
        self.stream.close()

    def dump_column(self, tab, colname, verbose=False):
        """Exports the specified column of the table"""
        nrows = tab.nrows()
        kws = tab.getcolkeywords(colname)
        # dump in <.5Gb chunks, so figure out size of one row
        row0 = tab.getcol(colname, 0, 1)
        maxrows = (2 ** 27) // (row0.size * row0.itemsize)
        rng = list(range(0, nrows, maxrows))
        # dump name,keywords,num_chunks
        for obj in colname, kws, len(rng):
            cPickle.dump(obj, self.stream)
        # dump chunks one by one
        for i0 in rng:
            cPickle.dump(tab.getcol(colname, i0, maxrows), self.stream)
        if verbose:
            print("  %s: %d chunk(s)%s%s" % (colname, len(rng),
                                       ", array column shape %s" % "x".join(
                                           map(str, row0.shape[1:])) if row0.ndim > 2 else "scalar column",
                                       ", %d keywords" % len(kws) if kws else ""))

    def load(self, tab, verbose):
        """Loads all columns from the dump file"""
        colnames = set(tab.colnames())
        while True:
            try:
                colname = cPickle.load(self.stream)
            except EOFError:
                break
            if not colname:
                break
            if colname in colnames:
                kws, nchunks = cPickle.load(self.stream), cPickle.load(self.stream)
                if kws is not None:
                    tab.putcolkeywords(colname, kws)
                row0 = 0
                for i in range(nchunks):
                    col = cPickle.load(self.stream)
                    #          print col.shape,row0
                    tab.putcol(colname, col, row0, col.shape[0])
                    row0 += col.shape[0]
                if verbose:
                    print("  %s: %d chunks%s" % (colname, nchunks, ", %d keywords" % len(kws) if kws else ""))
