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
  parser.add_option("--ampl",action="store_true",
                    help="makes amplitude plots");
  parser.add_option("--phase",action="store_true",
                    help="makes phase plots");
  parser.add_option("--pol",action="store_true",
                    help="makes instrumental polarization plots");
  parser.add_option("--per-corr",action="store_true",
                    help="includes per-correlation plots in the above");
  parser.add_option("--per-band",action="store_true",
                    help="includes combined per-band plots in the above");
  parser.add_option("--circle-pe",action="store_true",
                    help="makes circle plots of average pointing error per antenna");
  parser.add_option("--circle-ampl",action="store_true",
                    help="makes circle plots of average dE-amplitudes per antenna");
  parser.add_option("--circle-phase",action="store_true",
                    help="makes circle plots of average dE-phases per antenna");
  parser.add_option("--circle-ampl-ant",metavar="ANTENNA",type="string",
                    help="makes per-timeslot plots of average dE-amplitudes for the given antenna");
  parser.add_option("--circle-phase-ant",metavar="ANTENNA",type="string",
                    help="makes per-timeslot plots of average dE-phases for the given antenna");
  parser.add_option("--ampl-slope",action="store_true",
                    help="makes plots of a slope fit for dE-amplitudes");
  parser.add_option("--phase-slope",action="store_true",
                    help="makes plots of a slope fit for dE-phases");
  parser.add_option("--lsm",type="string",
                    help="LSM file (for source positions in circle plots.) Default is first *.lsm.html file found.");
  parser.add_option("-s","--sources",type="string",
                    help="source subset (default is ':' meaning all, use --list to list)");

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

  parser.set_defaults(circle_ampl_ant="",circle_phase_ant="",sources="",
    output_prefix="",output_type="png",papertype='a4');

  (options,args) = parser.parse_args();

  if not options.circle_pe and not options.circle_phase and \
      not options.circle_ampl_ant and not options.circle_phase_ant and \
      not options.ampl_slope and not options.phase_slope and \
      not options.ampl and not options.phase and not options.pol and \
      not options.list:
    parser.error("No plots specified.");

  if not args:
    parser.error("No parmtables specified.");

  import os.path
  import os
  import sys

  # load LSM file if needed
  if options.circle_ampl or options.circle_phase or options.circle_ampl_ant or options.circle_phase_ant:
    if not options.lsm:
      lsms = [ filename for filename in os.listdir(".") if filename.endswith(".lsm.html") ];
      if not lsms:
        parser.error("No LSMs found. Use --lsm to specify one explicitly.");
      lsm = lsms[0];
    else:
      lsm = options.lsm;
    # find tigger
    try:
      import Tigger;
    except ImportError:
      # make plot of average dE to model
      sys.path.append(os.getenv('HOME'));
      import Tigger
    print "Using LSM file %s"%lsm;
    model = Tigger.load(lsm);

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

  # vector of frequencies
  freq0 = numpy.array([0.]*len(SPWS));

  # scan funklet names to build up sets of keys
  pt = ParmTab(args[0]);
  for name in pt.funklet_names():
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


  def read_complex_parms (pt,prefix="dE"):
    """Reads complex c00 funklets from the given parmtable (given as filename).
    For each key in the 'keys' list, forms up keystring (by joining the elements
    of key weith ":"), then reads <prefix>:<keystring>:r and reads <prefix>:<keystring>:i,
    and forms up and returns a complex array.
    This assumes that the parm only changes in time.
    """
    arr = numpy.zeros((len(SRCS),len(ANTS),len(CORRS),NTIMES),dtype=complex);
    for i,src in enumerate(SRCS):
      for j,ant in enumerate(ANTS):
        for k,corr in enumerate(CORRS):
          fs_real = pt.funkset(':'.join([prefix,src,ant,corr,"r"])).get_slice();
          fs_imag = pt.funkset(':'.join([prefix,src,ant,corr,"i"])).get_slice();
          if len(fs_real) != len(fs_imag) or len(fs_real) != NTIMES:
            print "Error: table contains %d real and %d imaginary funklets; %d expected"%(len(fs_real),len(fs_imag),NTIMES);
            sys.exit(1);
          arr[i,j,k,:] = numpy.array([complex(c00(fr),c00(fi)) for fr,fi in zip(fs_real,fs_imag)]);
    return arr;


  # dl,dm is a 2 x NSPW x  NANT x NTIME array of poitning offsets
  dlm = numpy.zeros((2,len(SPWS),len(ANTS),NTIMES),dtype=float);

  for spw,tabname in enumerate(args):
    print "Reading",tabname;
    pt = ParmTab(tabname);
    freq0[spw] = sum(pt.envelope_domain().freq)/2;
    for i,ant in enumerate(ANTS):
      fsl = pt.funkset('E::dlm::dl:%s'%ant).get_slice();
      fsm = pt.funkset('E::dlm::dm:%s'%ant).get_slice();
      if len(fsl) != len(fsm) or len(fsl) != NTIMES:
        print "Error: table contains %d funklets for dl and %d for dm; %d expected"%(len(fsl),len(fsm),NTIMES);
        sys.exit(1);
      dlm[0,spw,i,:] = map(c00,fsl);
      dlm[1,spw,i,:] = map(c00,fsm);

  # convert to millidegrees
  dlm *= 180*1000/math.pi;

  # take mean and std along freq axis
  # these are now 2 x NANT x NTIME arrays
  dlm_mean = dlm.mean(1);
  dlm_std  = dlm.std(1);

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
        if mode is PLOT_SINGLE:
          plt.plot(range(len(data)) if xaxis is None else xaxis,data);
        elif mode is PLOT_MULTI:
          for dd in data:
            x,y,yerr = get_plot_data(dd,xaxis);
            if yerr is None:
              plt.plot(x,y);
            else:
              plt.errorbar(x,y,yerr,fmt=None,capsize=1);
        elif mode is PLOT_ERRORBARS:
          y,yerr = data;
          print y,yerr;
          plt.errorbar(range(len(y)) if xaxis is None else xaxis,y,yerr,fmt=None,capsize=1);
        if ylock:
          plt.set_ylim(ymin[irow,icol],ymax[irow,icol]);
          if ylock is not "col" and icol:
            plt.set_yticks([]);
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

  def make_skymap (ll,mm,
    markers,  # marker is a list of strings or dicts. For each marker, axes.plot() is invoked
              # as plot(x,y,str) or plot(x,y,**dict).
    labels=None,
    suptitle=None,      # title of plot
    save=None,          # filename to save to
    format=None,        # format: use options.output_type by default
    figsize=(290,210)       # figure width,height in mm
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
    plt = fig.add_axes([.05,.05,.9,.9]);
    plt.axhline(y=0,color='black',linestyle=':');
    plt.axvline(x=0,color='black',linestyle=':');

    if not isinstance(markers,(list,tuple)):
      markers = [markers]*len(ll);

    if not labels:
      labels = [None]*len(ll);

    ARCMIN = math.pi/(180*60);

    # set limits
    lim = max(max(abs(ll)),max(abs(mm)))*1.05/ARCMIN;

    # plot markers
    for l,m,mark,label in zip(ll/ARCMIN,mm/ARCMIN,markers,labels):
      if isinstance(mark,str):
        plt.plot(l,m,mark);
      elif isinstance(mark,dict):
        plt.plot(l,m,**mark);
      if label:
        plt.text(l,m-lim/50,label,fontsize=8,horizontalalignment='center',verticalalignment='top');

    # make circles at 30' and 1deg
    x = numpy.cos(numpy.arange(0,360)*math.pi/180);
    y = numpy.sin(numpy.arange(0,360)*math.pi/180);

    for R in numpy.arange(2)*30:
      plt.plot(x*R,y*R,linestyle=':',color='black');

    plt.set_xlim(-lim,lim);
    plt.set_ylim(-lim,lim);

    if suptitle:
      fig.suptitle(suptitle);
    if save:
      fig.savefig(save,papertype=options.papertype,
                  orientation='portrait' if figsize[0]<figsize[1] else 'landscape');
      print "Wrote",save;
      fig = None;
      pyplot.close("all");
    return save;

  sz_pnt = ( 290,60 );

  make_figure(enumerate(("dl","dm")),enumerate(ANTS),
    lambda i,iant:(dlm_mean[i,iant,:],dlm_std[i,iant,:]),
      hline=0,ylock=True,figsize=sz_pnt,mode=PLOT_ERRORBARS,
      suptitle="Pointing offset mean & stddev across all bands, millideg.",
      save="Epnt_mean");

  if options.output_type.upper() == "X11":
    from pylab import plt
    plt.show();
