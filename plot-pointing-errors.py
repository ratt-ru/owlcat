#!/usr/bin/python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser,OptionGroup
  parser = OptionParser(usage="""%prog: [plots & options] parmtables""",
      description="""Makes various plots of dE solutions.""");

  parser.add_option("-l","--list",action="store_true",
                    help="lists stuff found in MEP tables, then exits");
  parser.add_option("-c","--cache",metavar="FILENAME",type="string",
                    help="cache parms to file, which can be fed to plot-de-solutions script");

  group = OptionGroup(parser,"Output options");
  group.add_option("-o","--output-type",metavar="TYPE",type="string",
                  help="File format, see matplotlib documentation "
                  "for supported formats. At least 'png', 'pdf', 'ps', 'eps' and 'svg' are supported, or use 'x11' to display "
                  "plots interactively. Default is '%default.'");
  group.add_option("-r","--refresh",action="store_true",
                  help="Refresh plots even if they already exist (default is to keep existing plots.)");
  group.add_option("--output-prefix",metavar="PREFIX",type="string",
                    help="Prefix output filenames with PREFIX_");
  group.add_option("--papertype",dest="papertype",type="string",
                    help="set paper type (for .ps output only.) Default is '%default', but can also use e.g. 'letter', 'a3', etc.");
  parser.add_option_group(group);

  parser.set_defaults(
    output_prefix="",output_type="png",papertype='a4');

  (options,args) = parser.parse_args();

  if not args:
    parser.error("No parmtables specified.");

  import os.path
  import os
  import sys

  import numpy
  from ParmTables import ParmTab
  import matplotlib
  import math

  if options.output_type.upper() != "X11":
    matplotlib.use('agg');
  else:
    matplotlib.use('qt4agg');
  import matplotlib.pyplot as pyplot

  SPWS = range(len(args));

  # set of all sources, antennas and correlations
  ANTS = set();

  ## dict of WSRT antenna positions, relative to RT0.
  ##
  ANTX_dict = {
    '0':   0.0,        '1':  143.98881006, '2':  287.98154006, '3':   431.97187006,
    '4': 575.96470004, '5':  719.95646011, '6':  863.94757006, '7':  1007.93746007,
    '8':1151.92894011, '9': 1295.9213701 , 'A': 1331.92456019, 'B':  1403.91958018,
    'C':2627.84690046, 'D': 2699.84118052
  };

  # complex array of dEs per each src,antenna,corr tuple
  des = {};

  # scan funklet names to build up sets of keys
  pt = ParmTab(args[0]);
  for name in pt.funklet_names():
    if name.startswith("E::dl"):
      fields = name.split(':');
      ANTS.add(fields[-1]);
  NTIMES = len(pt.funkset(pt.funklet_names()[0]).get_slice());

  ANTS = sorted(ANTS);
  if options.list:
    print "MEP table %s contains pointing offsets for"%args[0];
    print "%d antennas: %s"%(len(ANTS)," ".join(ANTS));
    sys.exit(0);

  # check that antenna positions are known
  unknown_antennas = [ p for p in ANTS if p not in ANTX_dict ];
  if unknown_antennas:
    print "Don't have positions for antenna(s) %s"%",".join(unknown_antennas);
    sys.exit(0);

  # make vector of antenna positions
  ANTX = numpy.array([ ANTX_dict[p] for p in ANTS ]);
  # recenter at middle of the array
  ANTX -= ANTX[-1]/2;

  # check options that specify antennas by name
  for ant in 'circle_ampl_ant','circle_phase_ant':
    antname = getattr(options,ant,None);
    if antname:
      try:
        globals()[ant] = ANTS.index(antname);
      except IndexError:
        parser.error("Antenna name '%s' not found."%antname);
    else:
      globals()[ant] = None;

  def c00 (funklet):
    if numpy.isscalar(funklet.coeff):
      return funklet.coeff;
    else:
      return funklet.coeff.ravel()[0];

  # dl,dm is a 2 x NSPW x  NANT x NTIME array of poitning offsets
  dlm = numpy.zeros((2,len(SPWS),len(ANTS),NTIMES),dtype=float);

  # bsz is a 2 x 2 x NSPW x  NANT x NTIME array of beam sizes
  # first index is x/y, second is l/m
  bsz = numpy.zeros((2,2,len(SPWS),len(ANTS),NTIMES),dtype=float);
  beam_sizes = 0;

  for spw,tabname in enumerate(args):
    print "Reading",tabname;
    pt = ParmTab(tabname);
    for i,ant in enumerate(ANTS):
      # fill dlm
      fsl = pt.funkset('E::dl:%s'%ant).get_slice();
      fsm = pt.funkset('E::dm:%s'%ant).get_slice();
      if len(fsl) != len(fsm) or len(fsl) != NTIMES:
        print "Error: table contains %d funklets for dl and %d for dm; %d expected"%(len(fsl),len(fsm),NTIMES);
        sys.exit(1);
      dlm[0,spw,i,:] = map(c00,fsl);
      dlm[1,spw,i,:] = map(c00,fsm);
      # fill beam sizes
      if 'E::beamshape:%s'%ant in pt.funklet_names():
        beam_sizes = 1;
        fs = pt.funkset('E::beamshape:%s'%ant).get_slice();
        bsz[0,0,spw,i,:] = map(c00,fs);
      elif 'E::beamshape:xy:lm:%s'%ant in pt.funklet_names():
        beam_sizes = 4;
        for ixy,xy in enumerate("xy"):
          for ilm,lm in enumerate("lm"):
            fs = pt.funkset('E::beamshape:%s:%s:%s'%(ant,xy,lm)).get_slice();
            bsz[ixy,ilm,spw,i,:] = map(c00,fs);

  # write cache
  if options.cache:
    import cPickle
    cachefile = options.cache+'.cache';
    cPickle.dump((dlm,bsz,beam_sizes),file(cachefile,'w'));
    print "Cached all structures to file",cachefile;

  # convert dlm to millidegrees
  dlm *= 180*1000/math.pi;

  # take mean and std along freq axis
  # these are now 2 x NANT x NTIME arrays
  dlm_mean = dlm.mean(1);
  dlm_std  = dlm.std(1);
  bsz_mean = bsz.mean(2);
  bsz_std  = bsz.std(2);
  # and along time axis
  # these are now 2 x NSPW x NANT
  dlm_fmean = dlm.mean(3);
  dlm_fstd  = dlm.std(3);
  bsz_fmean = bsz.mean(4);
  bsz_fstd  = bsz.std(4);

  print "Read %d parmtables"%len(args);

  PLOT_SINGLE = 'single';
  PLOT_MULTI  = 'multi';
  PLOT_ERRORBARS = 'errorbars';

  def get_plot_data (data,xaxis=None):
    """Helper function, interprets the plot data returned by a datafunc.
    'data' is input data as returned by the datafunc argument to make_figure() below.
    'xaxis' is default X axis, Note to use ordinal numbering.
    Return value is tuple of X,Y,Yerr. Yerr is None if no error bars are provided.
    """;
    # a 2- or 3-tuple is interpreted as x,y[,yerr]
    if isinstance(data,tuple):
      if len(data) == 2:
        x,y = data;
        yerr = None;
      elif len(data) == 3:
        x,y,yerr = data;
      else:
        raise TypeError,"incorrect datum returned: 2- or 3-tuple expected, got %d-tuple"""%len(data);
    # else interpret data as Y
    else:
      y = data;
      x = yerr = None;
    # set X axis
    if x is None:
      x = range(len(y)) if xaxis is None else xaxis;
    # check lengths
    if len(x) != len(y):
      raise TypeError,"misshaped datum returned: %d X elements, %d Y elements"""%(len(x),len(y));
    if yerr is not None and len(yerr) != len(y):
      raise TypeError,"misshaped datum returned: %d Y elements, %d Yerr elements"""%(len(y),len(yerr));
    return x,y,yerr;

  def make_figure (rows,cols,  # (irow,row) and (icol,col) list
                  datafunc,   # datafunc(irow,icol) returns plot data for plot i,j
                  mode=PLOT_SINGLE,xaxis=None,
                  hline=None, # plot horizontal line at Y position, None for none
                  ylock="row", # lock Y scale. "row" locks across rows, "col" locks across columns, True locks across whole plot
                  mean_format="%.2f",
                  suptitle=None,      # title of plot
                  save=None,          # filename to save to
                  format=None,        # format: use options.output_type by default
                  figsize=(290,210)      # figure width,height in mm
                  ):
    if save and (format or options.output_type.upper()) != 'X11':
      save = "%s.%s"%(save,format or options.output_type);
      if options.output_prefix:
        save = "%s_%s.%s"%(options.output_prefix,save);
      # exit if figure already exists, and we're not refreshing
      if os.path.exists(save) and not options.refresh: # and os.path.getmtime(save) >= os.path.getmtime(__file__):
        print save,"exists, not redoing";
        return save;
    else:
      save = None;

    figsize = (figsize[0]/25.4,figsize[1]/25.4);
    fig = pyplot.figure(figsize=figsize,dpi=600);
    iplot = 0;
    rows = list(rows);
    cols = list(cols);
    if ylock:
      # form up ymin, ymax: NROWxNCOL arrays of min/max values per each plot
      if mode is PLOT_SINGLE:
        ymin = numpy.array([[datafunc(row[0],col[0]).min() for col in cols ] for row in rows ]);
        ymax = numpy.array([[datafunc(row[0],col[0]).max() for col in cols ] for row in rows]);
      elif mode is PLOT_MULTI:
        ymin = numpy.array([[ min([get_plot_data(dd,xaxis)[1].min() for dd in datafunc(row[0],col[0])]) for col in cols ] for row in rows ]);
        ymax = numpy.array([[ max([get_plot_data(dd,xaxis)[1].max() for dd in datafunc(row[0],col[0])]) for col in cols ] for row in rows ]);
      elif mode is PLOT_ERRORBARS:
        ymin = numpy.array([[ (datafunc(row[0],col[0])[0]-datafunc(row[0],col[0])[1]).min() for col in cols ] for row in rows]);
        ymax = numpy.array([[ (datafunc(row[0],col[0])[0]+datafunc(row[0],col[0])[1]).max() for col in cols ] for row in rows]);
      # now collapse them according to ylock mode
      if ylock == "row":
        ymin[...] = ymin.min(1)[:,numpy.newaxis];
        ymax[...] = ymax.max(1)[:,numpy.newaxis];
      elif ylock == "col":
        ymin[...] = ymin.min(0)[numpy.newaxis,:];
        ymax[...] = ymax.max(0)[numpy.newaxis,:];
      else:
        ymin[...] = ymin.min();
        ymax[...] = ymax.max();
    # make legend across the top
    for icol,col in cols:
      iplot += 1;
      plt = fig.add_subplot(len(rows)+1,len(cols)+1,iplot);
      plt.axis("off");
      plt.text(.5,0,col,fontsize='x-small',horizontalalignment='center',verticalalignment='bottom');
    iplot += 1;
    # now make plots
    for irow,row in rows:
      for icol,col in cols:
        iplot += 1;
        plt = fig.add_subplot(len(rows)+1,len(cols)+1,iplot);
        # plt.axis("off");
        plt.set_xticks([]);
        data = datafunc(irow,icol);
        realplot = True;
        if mode is PLOT_SINGLE:
          plt.plot(range(len(data)) if xaxis is None else xaxis,data);
          if xaxis is None:
            plt.set_xlim(-1,len(data));
        elif mode is PLOT_MULTI:
          for dd in data:
            x,y,yerr = get_plot_data(dd,xaxis);
            if yerr is None:
              plt.plot(x,y);
            else:
              plt.errorbar(x,y,yerr,fmt=None,capsize=1);
        elif mode is PLOT_ERRORBARS:
          y,yerr = data;
          if len(y) == 1:
            plt.text(0,0.6,mean_format%y[0],
                fontsize='x-small',horizontalalignment='left',verticalalignment='bottom');
            plt.text(0,0.5,"+/-"+(mean_format%yerr[0]),
                fontsize='x-small',horizontalalignment='left',verticalalignment='top');
            plt.axis("off");
            plt.set_ylim(0,1);
            realplot = False;
          else:
            plt.errorbar(range(len(y)) if xaxis is None else xaxis,y,yerr,fmt=None,capsize=1);
            if xaxis is None:
              plt.set_xlim(-1,len(y));
        if realplot:
          if ylock:
            plt.set_ylim(ymin[irow,icol],ymax[irow,icol]);
            if ylock is not "col" and icol:
              plt.set_yticklabels([]);
          if hline is not None:
            plt.axhline(y=hline,color='black');
        for lab in plt.get_yticklabels():
          lab.set_fontsize(5);
      # add row label
      iplot += 1;
      plt = fig.add_subplot(len(rows)+1,len(cols)+1,iplot);
      plt.text(0,.5,row,fontsize='x-small',horizontalalignment='left',verticalalignment='center');
      plt.axis("off");
    # adjust layout
    fig.subplots_adjust(left=0.05,right=.99,top=0.95,bottom=0.05);
    # plot if asked to
    if suptitle:
      fig.suptitle(suptitle);
    if save:
      fig.savefig(save,papertype=options.papertype,
                  orientation='portrait' if figsize[0]<figsize[1] else 'landscape');
      print "Wrote",save;
      fig = None;
      pyplot.close("all");
    return save;


  funcs = [
    lambda iant:(dlm_mean[0,iant,:],dlm_std[0,iant,:]),
    lambda iant:( numpy.array([dlm_mean[0,iant,:].mean()]),
                  numpy.array([dlm_mean[0,iant,:].std()])),
    lambda iant:(dlm_mean[1,iant,:],dlm_std[1,iant,:]),
    lambda iant:( numpy.array([dlm_mean[1,iant,:].mean()]),
                  numpy.array([dlm_mean[1,iant,:].std()])),
    lambda iant:(dlm_fmean[0,:,iant],dlm_fstd[0,:,iant]),
    lambda iant:(dlm_fmean[1,:,iant],dlm_fstd[1,:,iant])
  ];

  make_figure(enumerate(("dl","","dm","","dl, freq","dm, freq")),enumerate(ANTS),
        lambda i,iant:funcs[i](iant),
      hline=0,ylock=True,figsize=(290,150),mode=PLOT_ERRORBARS,
      suptitle="Pointing offset mean & stddev across all bands (top two plots) and times (bottom two plots), millideg.",
      save="Epnt_mean");

  for iant,ant in enumerate(ANTS):
    print "mean offset %s: %6.2f %6.2f"%(ant,dlm_mean[0,iant,:].mean(),dlm_mean[1,iant,:].mean());


  if beam_sizes == 4:
    funcs = [];
    for i0 in range(2):
      for j0 in range(2):
        funcs += [
          lambda iant,i=i0,j=j0:(bsz_mean[i,j,iant,:],bsz_std[i,j,iant,:]),
          lambda iant,i=i0,j=j0:( numpy.array([bsz_mean[i,j,iant,:].mean()]),
                        numpy.array([bsz_mean[i,j,iant,:].std()]))
        ];

    make_figure(enumerate(("Lx","","Mx","","Ly","","My","")),enumerate(ANTS),
          lambda i,iant:funcs[i](iant),
        hline=1,ylock=True,figsize=(290,210),mode=PLOT_ERRORBARS,
        mean_format="%.4f",
        suptitle="Beam extent in L/M, for X and Y dipoles, mean over time",
        save="Eshape_mean");

    funcs = [];
    for i0 in range(2):
      for j0 in range(2):
        funcs += [
          lambda iant,i=i0,j=j0:(bsz_fmean[i,j,:,iant],bsz_fstd[i,j,:,iant]),
          lambda iant,i=i0,j=j0:( numpy.array([bsz_fmean[i,j,:,iant].mean()]),
                        numpy.array([bsz_fmean[i,j,:,iant].std()]))
        ];
    make_figure(enumerate(("Lx","","Mx","","Ly","","My","")),enumerate(ANTS),
          lambda i,iant:funcs[i](iant),
        mean_format="%.4f",
        hline=1,ylock=True,figsize=(290,210),mode=PLOT_ERRORBARS,
        suptitle="Beam extent in L/M, for X/Y dipoles, mean over frequency",
        save="Eshape_mean_fq");
  elif beam_sizes == 1:
    funcs = [
      lambda iant:(bsz_mean[0,0,iant,:],bsz_std[0,0,iant,:]),
      lambda iant:( numpy.array([bsz_mean[0,0,iant,:].mean()]),
                    numpy.array([bsz_mean[0,0,iant,:].std()])),
      lambda iant:(bsz_fmean[0,0,:,iant],bsz_fstd[0,0,:,iant]),
    ];
    make_figure(enumerate(("size","","size fq")),enumerate(ANTS),
          lambda i,iant:funcs[i](iant),
        mean_format="%.4f",
        hline=1,ylock=True,figsize=(290,75),mode=PLOT_ERRORBARS,
        suptitle="Beam extent",
        save="Eshape");


  if options.output_type.upper() == "X11":
    from pylab import plt
    plt.show();
