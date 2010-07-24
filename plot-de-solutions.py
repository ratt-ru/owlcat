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

  if not options.circle_ampl and not options.circle_phase and \
      not options.circle_ampl_ant and not options.circle_phase_ant and \
      not options.ampl_slope and not options.phase_slope and \
      not options.ampl and not options.phase and not options.pol and \
      not options.list_sources:
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
  SRCS = set();
  ANTS = set();
  CORRS = set();

  # vector of antenna positions, relative to first antenna
  ANTX = numpy.array([    0.        ,   143.98881006,   287.98154006,   431.97187006,
  ## ant 5 is gone (DIGESTIF)
  #         575.96470004,   719.95646011,   863.94757006,  1007.93746007,
          575.96470004,   863.94757006,  1007.93746007,
          1151.92894011,  1295.9213701 ,  1331.92456019,  1403.91958018,
          2627.84690046,  2699.84118052])
  # recenter at middle of the array
  ANTX -= ANTX[-1]/2;

  # complex array of dEs per each src,antenna,corr tuple
  des = {};

  # vector of frequencies
  freq0 = numpy.array([0.]*len(SPWS));

  # scan funklet names to build up sets of keys
  pt = ParmTab(args[0]);
  for name in pt.funklet_names():
    label,src,ant,corr,ri = name.split(':');
    SRCS.add(src);
    ANTS.add(ant);
    CORRS.add(corr);
  NTIMES = len(pt.funkset(pt.funklet_names()[0]).get_slice());

  SRCS = sorted(SRCS);
  ANTS = sorted(ANTS);
  CORRS = sorted(CORRS);
  if options.list:
    print "MEP table %s contains dE's for:"%args[0];
    print "%d correlations: %s"%(len(CORRS)," ".join(CORRS));
    print "%d antennas: %s"%(len(ANTS)," ".join(ANTS));
    print "%d sources: %s"%(len(SRCS)," ".join(["%d:%s"%(i,src) for i,src in enumerate(SRCS)]));
    sys.exit(0);

  # source subset selection
  if options.sources:
    srcs = set();
    for ss in options.sources.split(","):
      if ss in SRCS:
        srcs.add(ss);
      else:
        try:
          ss = eval("SRCS[%s]"%ss);
        except:
          print "Error parsing source specification '%s'"%options.sources;
          sys.exit(1);
        if isinstance(ss,str):
          srcs.add(ss);
        else:
          srcs.update(ss);
    SRCS = sorted(srcs);
  
  print "Selected %d sources: %s"%(len(SRCS)," ".join(SRCS));

  XX = CORRS.index('xx');
  YY = CORRS.index('yy');
  # set of all source,antenna,corr combinations
  ALL = [ (src,ant,corr) for src in SRCS for ant in ANTS for corr in CORRS ];

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
          arr[i,j,k,:] = numpy.array([complex(fr.coeff,fi.coeff) for fr,fi in zip(fs_real,fs_imag)]);
    return arr;

  #
  # de is an NSPW x NSRC x NANT x NCORR x NTIME array of complex dEs
  # (NANT+1 actually, the last one is the mean-across-antennas value)
  #
  MEAN = len(ANTS);

  de = numpy.zeros((len(SPWS),len(SRCS),len(ANTS)+1,len(CORRS),NTIMES),dtype=complex);
  for spw,tabname in enumerate(args):
    print "Reading",tabname;
    pt = ParmTab(tabname);
    freq0[spw] = sum(pt.envelope_domain().freq)/2;
    de[spw,:,0:len(ANTS),:,:] = de1 = read_complex_parms(pt);

  # take mean along antenna axis
  de[:,:,MEAN,:,:] = de1 = numpy.mean(de[:,:,0:len(ANTS),:,:],2);

  # get phase component
  de_phase_rad = numpy.angle(de);
  de_phase = de_phase_rad*180/math.pi;
  de_ampl = abs(de);
  # average XX/YY dE phase (NSPWxNSRCxNANTxNTIME)
  deph = de_phase_rad.mean(3);
  # average spw/XX/YY dE phase (NSRCxNANTxNTIME)
  deph_mean = deph.mean(0);

  ### phase analysis
  # subtract mean phase over all antennas and times
  # dep0 is what's left, NSPW x NSRC x NANT x NCORR x NTIME
  dep0 = de_phase - (de_phase.mean(4).mean(2))[:,:,numpy.newaxis,:,numpy.newaxis];
  # cut outb the last ("mean") antenna
  dep0 = dep0[:,:,0:-1,:,:];

  # now fit phase offset and slope over array
  # dep0_b0/1 will be NSPWxNSRCxNCORRxNTIME
  sx = ANTX.sum()
  sx2 = (ANTX**2).sum();
  sy = dep0.sum(2)
  sxy = (dep0*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis,numpy.newaxis]).sum(2);
  n = len(ANTX);
  ## slope of linear fit
  dep0_b1 = (n*sxy - sy*sx)/(n*sx2-sx**2);
  ## offset of linear fit
  dep0_b0 = sy/n - dep0_b1*(sx/n);
  ## stddev w.r.t. linear fit
  #deph_fit_std = (deph - deph_b1*fq[:,numpy.newaxis,numpy.newaxis,numpy.newaxis] - deph_b0).std(0);

  # now also compute residuals w.r.t. the slope fits
  dep0_fitted = dep0_b0[:,:,numpy.newaxis,:,:] + \
    dep0_b1[:,:,numpy.newaxis,:,:]*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis,numpy.newaxis];
  dep0_res = dep0 - dep0_fitted;

  # same analysis for amplitudes
  dea0 = de_ampl - (de_ampl.mean(4).mean(2))[:,:,numpy.newaxis,:,numpy.newaxis];
  # cut out the last ("mean") antenna
  dea0 = dea0[:,:,0:-1,:,:];

  # now fit phase offset and slope over array
  # dep0_b0/1 will be NSPWxNSRCxNCORRxNTIME
  sy = dea0.sum(2)
  sxy = (dea0*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis,numpy.newaxis]).sum(2);
  n = len(ANTX);
  ## slope of linear fit
  dea0_b1 = (n*sxy - sy*sx)/(n*sx2-sx**2);
  ## offset of linear fit
  dea0_b0 = sy/n - dea0_b1*(sx/n);
  ## stddev w.r.t. linear fit
  #deph_fit_std = (deph - deph_b1*fq[:,numpy.newaxis,numpy.newaxis,numpy.newaxis] - deph_b0).std(0);

  # now also compute residuals w.r.t. the slope fits
  dea0_fitted = dea0_b0[:,:,numpy.newaxis,:,:] + \
    dea0_b1[:,:,numpy.newaxis,:,:]*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis,numpy.newaxis];
  dea0_res = dea0 - dea0_fitted;


  ## fitr offset and slope of mean (in freq) phases
  #dep0m = dep0.mean(0);   # NSRCxNANTxNCORRxNTIME
  #sy = dep0.sum(0)
  #sxy = (dep0*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis,numpy.newaxis]).sum(0);
  #n = len(ANTX);
  ## slope of linear fit
  #dep0_b1 = (n*sxy - sy*sx)/(n*sx2-sx**2);
  ## offset of linear fit
  #dep0_b0 = sy/n - dep0_b1*(sx/n);


  ### TEC analysis
  ## convert to dTEC corresponding to each spw's reference freq
  #K = (1.3445*2*math.pi*1e+9)   # radians*GHz/TECU
  #dtec = de_phase_rad*freq0[:,numpy.newaxis,numpy.newaxis,numpy.newaxis,numpy.newaxis]/K;
  ## get mean TEC across corrs and spws (NSRCxNANTxNTIME)
  #dtec_mean = dtec.mean(3).mean(0);
  ## stddev of dE phase - mean dE phase in freq (NSRCxNANTxNTIME)
  #deph_mean_std = (deph - deph_mean).std(0);
  ## stddev of dE phase - dE phase given by mean TEC (NSRCxNANTxNTIME)
  #deph_tec_std = (deph - dtec_mean*K/freq0[:,numpy.newaxis,numpy.newaxis,numpy.newaxis]).std(0);
  ## if slope fits better, then deph_fromtec will be consistently lower than deph_meanspw. This is then the
  ## confidence ratio (NSRCxNANTxNTIME)
  ## conf_ratio = deph_mean_std/deph_fromtec;

  ### try a linear fit in freq to the phases
  #fq = freq0*1e-9;
  #sx = fq.sum()
  #sx2 = (fq**2).sum();
  #sy = deph.sum(0)
  #sxy = (deph*fq[:,numpy.newaxis,numpy.newaxis,numpy.newaxis]).sum(0);
  #n = len(fq);
  ## slope of linear fit
  #deph_b1 = (n*sxy - sy*sx)/(n*sx2-sx**2);
  ## offset of linear fit
  #deph_b0 = sy/n - deph_b1*(sx/n);
  ## stddev w.r.t. linear fit
  #deph_fit_std = (deph - deph_b1*fq[:,numpy.newaxis,numpy.newaxis,numpy.newaxis] - deph_b0).std(0);

  # now normalize each dE by mean dE across all times and antennas (this takes out
  # corrections for due to imperfect LSMs)
  de /= numpy.mean(de1,3)[:,:,numpy.newaxis,:,numpy.newaxis];

  # add "mean" to antenna labels
  ANTS += [ "mean" ];

  # phase
  #de_phase = numpy.angle(de)*180/math.pi;

  # de_norm[spw,src,ant,time] is sqrt(|dExx|^2+|dEyy|^2), this is the "I" gain
  de_norm = numpy.sqrt((abs(de[:,:,:,XX,:])**2+abs(de[:,:,:,YY,:])**2)/2);
  de_norm[:,:,MEAN,:] = numpy.mean(de_norm[:,:,0:len(ANTS),:],2);
  # de_pol[spw,src,ant,time] is sqrt(|dExx|^2-|dEyy|^2), this is the instrumental Q polarization
  # (subtracted mean over antennas and times to take out LSM bias)
  de_pol  = (abs(de[:,:,:,XX,:])**2-abs(de[:,:,:,YY,:])**2)/2;
  de_pol[:,:,MEAN,:] = numpy.mean(de_pol[:,:,0:len(ANTS),:],2);
  de_pol -= numpy.mean(de_pol[:,:,MEAN,:],2)[:,:,numpy.newaxis,numpy.newaxis];
  # de_fpol[spw,src,ant,time] is de_norm/de_polz, i.e. fractional instrumental Q
  de_fpol = de_pol/(de_norm**2);
  # de_renorm[spw,src,ant,time] is de_norm divided by the mean de_norm[spw,src,:,time],
  # i.e. mean gain across all antennas
  de_renorm = de_norm/((de_norm[:,:,MEAN,:])[:,:,numpy.newaxis,:]);

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

  # figure sizes
  # more sources than antennas
  if len(SRCS) < (len(ANTS)+1):
    sz_as = ( 210,min(290,30*len(ANTS)));
    sz_sa = ( 290,min(210,30*len(SRCS)));
  else:
    sz_as = ( 290,min(210,30*len(ANTS)));
    sz_sa = ( 210,min(290,30*len(SRCS)));


  def make_figure (rows,cols,  # (irow,row) and (icol,col) list
                  datafunc,   # datafunc(irow,icol) returns plot data for plot i,j
                  mode=PLOT_SINGLE,xaxis=None,
                  hline=None, # plot horizontal line at Y position, None for none
                  ylock="row", # lock Y scale. "row" locks across rows, "col" locks across columns, True locks across whole plot
                  suptitle=None,      # title of plot
                  save=None,          # filename to save to
                  format=None,        # format: use options.output_type by default
                  figsize=sz_as       # figure width,height in mm
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
          plt.errorbar(range(len(y)) if xaxis is None else xaxis,y,yerr,fmt=None,capsize=1);
        if ylock:
          plt.set_ylim(ymin[irow,icol],ymax[irow,icol]);
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
    figsize=sz_as       # figure width,height in mm
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

  if options.per_corr:
    for corr in XX,YY:
      if options.phase:
        if options.per_band:
          make_figure(enumerate(ANTS),enumerate(SRCS),
            lambda iant,isrc:[de_phase[spw,isrc,iant,corr,:] for spw in SPWS],
            hline=0,ylock="row",mode=PLOT_MULTI,
            suptitle="dE%s phase (deg) per band"%CORRS[corr],
            save="dEphase_%s_perband"%CORRS[corr]);
        #make_figure(enumerate(ANTS),enumerate(SRCS),
          #lambda iant,isrc:(numpy.mean(de_phase[:,isrc,iant,corr,:],0),
                            #numpy.std(de_phase[:,isrc,iant,corr,:],0)),
          #hline=0,ylock="row",mode=PLOT_ERRORBARS,
          #suptitle="dE%s phase (deg) mean & stddev across all bands"%CORRS[corr],
          #save="dEphase_%s_mean"%CORRS[corr]);
        make_figure(enumerate(SRCS),enumerate(ANTS),
          lambda isrc,iant:(numpy.mean(de_phase[:,isrc,iant,corr,:],0),
                            numpy.std(de_phase[:,isrc,iant,corr,:],0)),
          hline=0,ylock="row",figsize=sz_sa,mode=PLOT_ERRORBARS,
          suptitle="dE%s phase (deg) mean & stddev across all bands"%CORRS[corr],
          save="dEphase_%s_mean"%CORRS[corr]);
      if options.ampl:
        make_figure(enumerate(SRCS),enumerate(ANTS),
          lambda isrc,iant:(numpy.mean(de_ampl[:,isrc,iant,corr,:],0),
                            numpy.std(de_ampl[:,isrc,iant,corr,:],0)),
          hline=1,ylock="row",figsize=sz_sa,mode=PLOT_ERRORBARS,
          suptitle="dE%s amplitude mean & stddev across all bands"%CORRS[corr],
          save="dEampl_%s_mean"%CORRS[corr]);
      #make_figure(enumerate(SRCS),enumerate(ANTS),
        #lambda isrc,iant:(numpy.mean(dtec[:,isrc,iant,corr,:],0),
                          #numpy.std(dtec[:,isrc,iant,corr,:],0)),
        #hline=0,ylock="row",figsize=sz_sa,mode=PLOT_ERRORBARS,
        #suptitle="dETEC(%s) mean & stddev across all bands"%CORRS[corr],
        #save="dETEC_%s_mean"%CORRS[corr]);
      #make_figure(enumerate(ANTS[:-1]),enumerate(SRCS),
        #lambda iant,isrc:(dep0_res[:,isrc,iant,corr,:].mean(0),
                          #dep0_res[:,isrc,iant,corr,:].std(0)),
        #hline=0,ylock="row",mode=PLOT_ERRORBARS,
        #suptitle="dE%s residual phase (w.r.t. slope fit) mean & stddev across all bands"%CORRS[corr],
        #save="dEphase_array_%s_res"%CORRS[corr]);
      if options.phase_slope:
        make_figure(enumerate(SRCS),enumerate(ANTS[:-1]),
          lambda isrc,iant:(dep0_res[:,isrc,iant,corr,:].mean(0),
                            dep0_res[:,isrc,iant,corr,:].std(0)),
          hline=0,ylock="row",mode=PLOT_ERRORBARS,figsize=sz_sa,
          suptitle="dE%s residual phase (w.r.t. slope fit) mean & stddev across all bands"%CORRS[corr],
          save="dEphase_array_%s_res"%CORRS[corr]);
      if options.ampl_slope:
      #make_figure(enumerate(ANTS[:-1]),enumerate(SRCS),
        #lambda iant,isrc:(dea0_res[:,isrc,iant,corr,:].mean(0),
                          #dea0_res[:,isrc,iant,corr,:].std(0)),
        #hline=0,ylock="row",mode=PLOT_ERRORBARS,
        #suptitle="dE%s residual ampl (w.r.t. slope fit) mean & stddev across all bands"%CORRS[corr],
        #save="dEampl_array_%s_res"%CORRS[corr]);
        make_figure(enumerate(SRCS),enumerate(ANTS[:-1]),
          lambda isrc,iant:(dea0_res[:,isrc,iant,corr,:].mean(0),
                            dea0_res[:,isrc,iant,corr,:].std(0)),
          hline=0,ylock="row",mode=PLOT_ERRORBARS,figsize=sz_sa,
          suptitle="dE%s residual ampl (w.r.t. slope fit) mean & stddev across all bands"%CORRS[corr],
          save="dEampl_array_%s_res"%CORRS[corr]);

  # dep0_b0/1 will be NSPWxNSRCxNCORRxNTIME
  if options.phase_slope:
    dep0s = [ dep0_b1[:,:,XX,:],dep0_b1[:,:,YY,:],dep0_b0[:,:,XX,:],dep0_b0[:,:,YY,:] ];
    make_figure(enumerate(("X slope deg/m","Y slope deg/m","X offset, deg","Y offset, deg")),enumerate(SRCS),
      lambda i,isrc:(numpy.mean(dep0s[i][:,isrc,:],0),
                    numpy.std(dep0s[i][:,isrc,:],0)),
      ylock=False,mode=PLOT_ERRORBARS,
      suptitle="Fitted differential phase slope over array",
      save="dEphase_array_slopes");

  if options.ampl_slope:
    dea0s = [ dea0_b1[:,:,XX,:],dea0_b1[:,:,YY,:],dea0_b0[:,:,XX,:],dea0_b0[:,:,YY,:] ];
    make_figure(enumerate(("X slope 1/m","Y slope 1/m","X offset","Y offset")),enumerate(SRCS),
      lambda i,isrc:(numpy.mean(dea0s[i][:,isrc,:],0),
                    numpy.std(dea0s[i][:,isrc,:],0)),
      ylock=False,mode=PLOT_ERRORBARS,
      suptitle="Fitted differential ampl slope over array",
      save="dEampl_array_slopes");

  #TIMELIST = [(h*2,"%dh"%h) for h in range(10)];

  #make_figure(TIMELIST,enumerate(SRCS),
    #lambda itime,isrc:(numpy.mean(dep0[:,isrc,:,XX,itime],0),
                      #numpy.mean(dep0[:,isrc,:,YY,itime],0)),
    #hline=0,ylock="row",mode=PLOT_MULTI,xaxis=ANTX,
    #suptitle="Differential XX+YY phase (deg) over array, evolution in time",
    #save="dEphase_array");

  #make_figure(TIMELIST,enumerate(SRCS),
    #lambda itime,isrc:(numpy.mean(dep0[:,isrc,:,XX,itime],0),
                      #numpy.std(dep0[:,isrc,:,XX,itime],0)),
    #hline=0,ylock="row",mode=PLOT_ERRORBARS,xaxis=ANTX,
    #suptitle="Differential XX phase (deg) over array, evolution in time",
    #save="dEphase_array_xx");

  #make_figure(TIMELIST,enumerate(SRCS),
    #lambda itime,isrc:(numpy.mean(dep0[:,isrc,:,YY,itime],0),
                      #numpy.std(dep0[:,isrc,:,YY,itime],0)),
    #hline=0,ylock="row",mode=PLOT_ERRORBARS,xaxis=ANTX,
    #suptitle="Differential YY phase (deg) over array, evolution in time",
    #save="dEphase_array_yy");

  #make_figure(TIMELIST,enumerate(SRCS),
    #lambda itime,isrc:(numpy.mean(dea0[:,isrc,:,XX,itime],0),
                      #numpy.mean(dea0[:,isrc,:,YY,itime],0)),
    #hline=0,ylock="row",mode=PLOT_MULTI,xaxis=ANTX,
    #suptitle="Differential XX+YY ampl (deg) over array, evolution in time",
    #save="dEampl_array");

  #make_figure(TIMELIST,enumerate(SRCS),
    #lambda itime,isrc:(numpy.mean(dea0[:,isrc,:,XX,itime],0),
                      #numpy.std(dea0[:,isrc,:,XX,itime],0)),
    #hline=0,ylock="row",mode=PLOT_ERRORBARS,xaxis=ANTX,
    #suptitle="Differential XX ampl (deg) over array, evolution in time",
    #save="dEampl_array_xx");

  #make_figure(TIMELIST,enumerate(SRCS),
    #lambda itime,isrc:(numpy.mean(dea0[:,isrc,:,YY,itime],0),
                      #numpy.std(dea0[:,isrc,:,YY,itime],0)),
    #hline=0,ylock="row",mode=PLOT_ERRORBARS,xaxis=ANTX,
    #suptitle="Differential YY ampl (deg) over array, evolution in time",
    #save="dEampl_array_yy");


  if options.pol:
    make_figure(enumerate(SRCS),enumerate(ANTS),
      lambda isrc,iant:(numpy.mean(de_fpol[:,isrc,iant,:],0),numpy.std(de_fpol[:,isrc,iant,:],0)),
      hline=0,ylock="row",mode=PLOT_ERRORBARS,
      suptitle="Fractional instrumental polarization, average across all bands",
      save="dEpol_mean");
    if options.per_band:
      make_figure(enumerate(SRCS),enumerate(ANTS),
        lambda isrc,iant:[de_fpol[spw,isrc,iant,:] for spw in SPWS],
        hline=0,ylock="row",mode=PLOT_MULTI,
        suptitle="Fractional instrumental polarization per band",
        save="dEpol_mean_sa_perband");

  if options.ampl:
    make_figure(enumerate(SRCS),enumerate(ANTS[:-1]),  # omit MEAN antenna since it's 1 by definition
      lambda isrc,iant:(numpy.mean(de_renorm[:,isrc,iant,:],0),numpy.std(de_renorm[:,isrc,iant,:],0)),
      hline=1,ylock="row",figsize=sz_sa,mode=PLOT_ERRORBARS,
      suptitle="Renormalized ||dE|| mean & stddev across all bands",
      save="dE_mean");
    if options.per_band:
      make_figure(enumerate(SRCS),enumerate(ANTS),
        lambda isrc,iant:[de_renorm[spw,isrc,iant,:] for spw in SPWS],
        hline=1,ylock="row",figsize=sz_sa,mode=PLOT_MULTI,
        suptitle="Renormalized ||dE|| per band",
        save="dE_mean_perband");

  #make_figure(enumerate(SRCS),enumerate(ANTS),
    #lambda isrc,iant:(numpy.mean(de_norm[:,isrc,iant,:],0),numpy.std(de_norm[:,isrc,iant,:],0)),
    #hline=1,ylock="row",figsize=sz_sa,mode=PLOT_ERRORBARS,
    #suptitle="||dE|| mean & stddev across all bands",
    #save="dE_mean");


  #make_figure(enumerate(ANTS),enumerate(SRCS),
    #lambda iant,isrc:[numpy.mean(abs(de[:,isrc,iant,corr,:]),0) for corr in XX,YY],
    #hline=1,ylock="row",mode=PLOT_MULTI,
    #suptitle="|dExx|,|dEyy| averaged across all bands",
    #save="dE_mean_xx_yy");

  if options.circle_ampl or options.circle_phase or options.circle_ampl_ant or options.circle_phase_ant:
    maxabsdelog = math.sqrt(0.5); # abs(delog).max();
    maxsize = 36
    minsize = 0
    lsrc = numpy.zeros(len(SRCS),float);
    msrc = numpy.zeros(len(SRCS),float);
    for i,name in enumerate(SRCS):
      src = model.findSource(name);
      lsrc[i],msrc[i] = src._lm_ncp;

  if options.circle_ampl:
    # de_renorm is NSPWxNSRCxNANTxNTIME
    # reshape as NSRCxNANTxNSPWxNTIME
    de_r1 = de_renorm.copy().transpose((1,2,0,3));
    demean = de_r1.mean(2).mean(2);
    # stddev is sqrt(s1^2+s2^2), where s1 is mean error, and s2 is error of mean
    destd  = numpy.sqrt((de_r1.std(2).mean(2)/math.sqrt(NTIMES-1))**2 + de_r1.mean(2).std(2)**2);
    delog = numpy.where(demean>=1,numpy.sqrt(demean-1),-numpy.sqrt(1-demean));

    antplots = [];

    for iant,ant in list(enumerate(ANTS))[:-1]:
      # make figure
      markers = [];
      labels = [];
      for name,log,mean,std in zip(SRCS,delog[:,iant],demean[:,iant],destd[:,iant]):
        color = "blue" if log >=0 else "red";
        mark = dict(marker='o',markersize=minsize + (maxsize-minsize)*(abs(log)/maxabsdelog),
                    markeredgecolor=color,markerfacecolor='None');
        sigmas = abs(1-mean)/std;
        if sigmas >= 3:
          mark['markerfacecolor'] = color;
        elif sigmas >= 2:
          mark['markeredgewidth'] = 3;
        elif sigmas >= 1:
          mark['markeredgewidth'] = 2;
        else:
          mark['markeredgewidth'] = 1;
        markers.append(mark);
        labels.append("%s %.2f+-%.2f"%(name,mean,std));

      antplots.append(make_skymap(lsrc,msrc,markers,labels=labels,
        figsize=(210,210),suptitle="Average ||dE||, RT%s"%ant,save="dE_ant%s"%ant));
    
    # try to make summary plot
    if options.output_type.upper() == "PNG":
      try:
        import PIL
        import PIL.Image
      except ImportError:
        print "Python Imaging Library (PIL) not found, can't make the antenna rogues' gallery";
        PIL = None;
      if PIL:
        ncol = 5;  # number of plots per column
        plotsize = 256; # column width, in pixels
        ncol = min(ncol,len(antplots));
        nrow,rem = divmod(len(antplots),ncol);
        if rem:
          nrow += 1;
        # make image
        img = PIL.Image.new('RGBA',(ncol*plotsize,nrow*plotsize));
        for i,plot in enumerate(antplots):
          subimg = PIL.Image.open(plot);
          subimg = subimg.resize((plotsize,plotsize),PIL.Image.ANTIALIAS);
          y0,x0 = divmod(i,ncol);
          img.paste(subimg,(x0*plotsize,y0*plotsize));
        img.save("dE_ant_gallery.png","PNG");

  if circle_ampl_ant is not None:
    # make animation for given antenna
    antname = ANTS[circle_ampl_ant];

    # de_renorm is NSPWxNSRCxNANTxNTIME
    demean = de_renorm[:,:,circle_ampl_ant,:].mean(0);
    destd  = de_renorm[:,:,circle_ampl_ant,:].std(0);
    delog = numpy.sqrt(demean);  # numpy.log10(demean);
    delog = numpy.where(demean>=1,numpy.sqrt(demean-1),-numpy.sqrt(1-demean));

    frames = [];
    for k in range(NTIMES):
      markers = [];
      labels = [];
      for name,log,mean,std in zip(SRCS,delog[:,k],demean[:,k],destd[:,k]):
        color = "blue" if log >=0 else "red";
        mark = dict(marker='o',markersize=minsize + (maxsize-minsize)*(abs(log)/maxabsdelog),
                    markeredgecolor=color,markerfacecolor='None');
        sigmas = abs(1-mean)/std;
        if sigmas >= 3:
          mark['markerfacecolor'] = color;
        elif sigmas >= 2:
          mark['markeredgewidth'] = 3;
        elif sigmas >= 1:
          mark['markeredgewidth'] = 2;
        else:
          mark['markeredgewidth'] = 1;
        markers.append(mark);
        labels.append("%s %.2f+-%.2f"%(name,mean,std));

      frames.append(make_skymap(lsrc,msrc,markers,labels=labels,
        figsize=(210,210),suptitle="Average ||dE||, antenna %s, time slice %d"%(antname,k),save="dE_ant%s_%03d"%(antname,k)));
      
    # try to make animation
    if options.output_type.upper() == "PNG":
      for f in frames:
        os.system("pngtopnm %s | ppmquant 256 | ppmtogif >%s.gif"%(f,f));
      os.system("gifsicle -d 15 -l0 --colors 256 %s >%s"%
                (" ".join([f+".gif" for f in frames + frames[-1::-1]]),
                  "dE_ant%s_anim.gif"%antname));
      os.system("rm -f %s"%" ".join([f+".gif" for f in frames]));


  if options.output_type.upper() == "X11":
    from pylab import plt
    plt.show();
