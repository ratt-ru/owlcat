#!/usr/bin/python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser,OptionGroup
  parser = OptionParser(usage="""%prog: [plots & options] <parmtable(s) or cache file>""",
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
                    help="makes plots of an East-West slope fit for dE-amplitudes");
  parser.add_option("--phase-slope",action="store_true",
                    help="makes plots of a East-West slope fit for dE-phases");
  parser.add_option("--phase-slope-rot-offset",metavar="DEG",type="float",
                    help="...including phase slopes corresponding to rotational offset (EXPERIMENTAL)");
  parser.add_option("--phase-slope-time-offset",metavar="FRAC",type="float",
                    help="...including phase slopes corresponding to time offset (EXPERIMENTAL)");
  parser.add_option("--phase-slope-res",action="store_true",
                    help="makes plots of residuals w.r.t. slope fit for dE-phases");
  parser.add_option("--phase-slope-dlm",action="store_true",
                    help="makes plots of dl/dm offsets corresponding to phase slopes (EXPERIMENTAL)");
  parser.add_option("--lsm",type="string",
                    help="LSM file (for source positions where needed). Default is first *.lsm.html file found.");
  parser.add_option("--update-lsm",action="store_true",
                    help="updates LSM file based on dE solutions (EXPERIMENTAL)");
  parser.add_option("--ms",type="string",
                    help="Measurement Set (for UVW coordinates, antenna positions and such)");
  parser.add_option("-s","--sources",type="string",
                    help="source subset (default is ':' meaning all, use --list to list)");
  parser.add_option("-c","--cache",metavar="FILENAME",type="string",
                    help="cache parms to file, which can be given instead of parmtables next time, for quicker startup");
  parser.add_option("--pe",metavar="FILENAME",type="string",
                    help="cache file from plot-pointing-errors.py script. If given, |dE| given by the pointing model will be plotted");

  group = OptionGroup(parser,"Plotting options");
  group.add_option("--radius",type="int",
                    help="radius of circle plots, in arcmin. Default is auto-sizing.");
  group.add_option("-G","--grid-circle",metavar="ARCMIN",type="int",action="append",
                    help="sets radius of grid circle(s) in circle plots. May be given multiple times. Default is 30 and 60.");
  group.add_option("--title-fontsize",metavar="POINTS",type="int",
                  help="Set plot title font size, 0 for no title. Default is %default.");
  group.add_option("--subtitle-fontsize",metavar="POINTS",type="int",
                  help="Set plot subtitle font size, 0 for none. Default is %default.");
  group.add_option("--label-fontsize",metavar="POINTS",type="int",
                  help="Set plot label font size, 0 for no labels. Default is %default.");
  group.add_option("--axis-fontsize",metavar="POINTS",type="int",
                  help="Set axis label font size, 0 for no axis labels. Default is %default.");
  group.add_option("--circle-title-fontsize",metavar="POINTS",type="int",
                  help="Set plot title font size in circle plots, 0 for no title. Default is %default.");
  group.add_option("--circle-label-fontsize",metavar="POINTS",type="int",
                  help="Set plot label font size in circle plots, 0 for no labels. Default is %default.");
  group.add_option("--circle-axis-fontsize",metavar="POINTS",type="int",
                  help="Set axis label font size in circle_plots, 0 for no axis labels. Default is %default.");
  group.add_option("--names-only",action="store_true",
                  help="Use only source names for labels in circle plots. Default is name, mean and std values.");
  group.add_option("--circle-minsize",metavar="POINTS",type="int",
                  help="Set minimum size of marker in circle plots. Default is %default.");
  group.add_option("--circle-maxsize",metavar="POINTS",type="int",
                  help="Set maximum size of marker in circle plots. Default is %default.");
  group.add_option("--circle-linewidth",metavar="SCALE",type="int",
                  help="Max linewidth used to draw open circles in circle plots. The width corresponds to 3 sigmas. Default is %default.");
  group.add_option("--circle-label-offset",metavar="FRAC",type="float",
                  help="Vertical offset of circle labels. Use 0 to center labels on source position. Default is %default.");
  group.add_option("--dlm-scale",metavar="N",type="float",
                  help="Scale for arrows indicating dl/dm offset, as N:1. Default is to auto-scale based on max length.");
  group.add_option("--borders",metavar="L,R,B,T",type="str",
                  help="Set left/right/bottom/top plot borders, in normalized page coordinates. Default is %default.");
  group.add_option("--circle-borders",metavar="L,R,B,T",type="str",
                  help="Set left/right/bottom/top circle-plot borders, in normalized page coordinates. Default is %default.");
  group.add_option("--subplot-wspace",metavar="W",type="float",
                  help="Set spacing between subplots, as fraction of subplot width. Default is automatic.");
  parser.add_option_group(group);

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
  group.add_option("-W","--width",type="int",
                    help="set explicit plot width, in mm. (Useful for .eps output)");
  group.add_option("-H","--height",type="int",
                    help="set explicit plot height, in mm. (Useful for .eps output)");
  group.add_option("--portrait",action="store_true",
                  help="Force portrait orientation. Default is to select orientation based on plot size");
  group.add_option("--landscape",action="store_true",
                  help="Force landscape orientation.");
  parser.add_option_group(group);

  parser.set_defaults(circle_ampl_ant="",circle_phase_ant="",sources="",cache=None,
    borders="0.05,0.99,0.05,0.95",circle_borders=".05,.9,.05,.9",subplot_wspace=0,
    title_fontsize=12,subtitle_fontsize=10,label_fontsize=8,axis_fontsize=5,
    circle_title_fontsize=12,circle_label_fontsize=8,circle_axis_fontsize=5,
    phase_slope_rot_offset=0,
    radius=0,dlm_scale=0,
    circle_maxsize=36,circle_minsize=0,circle_linewidth=5,circle_label_offset=0.02,
    output_prefix="",output_type="png",papertype='a4',width=0,height=0);

  (options,args) = parser.parse_args();

  if not options.circle_ampl and not options.circle_phase and \
      not options.circle_ampl_ant and not options.circle_phase_ant and \
      not options.ampl_slope and not options.phase_slope and not options.phase_slope_res and \
      not options.ampl and not options.phase and not options.pol and \
      not options.list:
    parser.error("No plots specified.");

  if not args:
    parser.error("No parmtables or cachefile specified.");

  if (options.width or options.height) and not (options.width and options.height):
    parser.error("Either both -W/--width and -H/--height must be specified, or neither.");

  try:
    borders = [ float(x) for x in options.borders.split(",") ];
  except:
    borders = [];
  if len(borders) != 4:
    parser.error("Bad --borders specification.");
  try:
    circle_borders = [ float(x) for x in options.circle_borders.split(",") ];
  except:
    circle_borders = [];
  if len(circle_borders) != 4:
    parser.error("Bad --circle-borders specification.");

  if not options.grid_circle:
    options.grid_circle = (30,60);

  import os.path
  import os
  import cPickle
  import sys
  # load other stuff
  import numpy
  from scipy import linalg
  from Owlcat import Coordinates

  import matplotlib
  import math
  if options.output_type.upper() != "X11":
    matplotlib.use('agg');
  else:
    matplotlib.use('qt4agg');
  import matplotlib.pyplot as pyplot
  from matplotlib.ticker import MaxNLocator

  DEG = math.pi/180;
  ARCMIN = math.pi/(180*60);

  #=== loads MS (that has been specified in options, else picks one in current dir)
  #=== returns tuple of ms,msant, where the first is the main table, and the second is the ANTENNA subtable
  _ms = _msant = None;
  def load_ms ():
    global _ms;
    global _msant;
    if _ms is None:
      import pyrap.tables
      # if MS file not specified, use first one found
      if not options.ms:
        mss = [ filename for filename in os.listdir(".") if filename.lower().endswith(".ms") ];
        if not mss:
          parser.error("No MSs found. Use --ms to specify one explicitly.");
        msname = mss[0];
      else:
        msname = options.ms;
      print "Using MS %s for coordinates information"%msname;
      _ms = pyrap.tables.table(msname);
      _msant = pyrap.tables.table(_ms.getkeyword('ANTENNA'));
    return _ms,_msant;

  #
  #========================== LOAD LSM FILE AND MS
  #
  if options.circle_ampl or options.circle_phase or options.circle_ampl_ant or options.circle_phase_ant \
      or options.phase_slope_rot_offset or options.phase_slope_time_offset or options.phase_slope_dlm \
      or options.update_lsm or (options.pe and options.ampl):
    # if LSM file not specified, use first one found
    if not options.lsm:
      lsms = [ filename for filename in os.listdir(".") if filename.endswith(".lsm.html") ];
      if not lsms:
        parser.error("No LSMs found. Use --lsm to specify one explicitly.");
      lsm_filename = lsms[0];
    else:
      lsm_filename = options.lsm;
    #
    # find tigger
    try:
      import Tigger;
    except ImportError:
      # make plot of average dE to model
      sys.path.append(os.getenv('HOME'));
      import Tigger
    print "Using LSM file %s"%lsm_filename;
    model = Tigger.load(lsm_filename);
    lsm_timestamp = os.path.getmtime(lsm_filename);

    # if MS file not specified, use first one found
    ms,msant = load_ms();
    nant = msant.nrows();
    # read UVWs
    uvw0 = ms.query('ANTENNA1==0 && ANTENNA2==%d'%(nant-1)).getcol('UVW');
    uvw = uvw0[30::60];
    # read phase center
    radec0 = pyrap.tables.table(ms.getkeyword('FIELD')).getcol('PHASE_DIR')[0][0];
    print "Phase center is at",radec0;
  else:
    model = None;

  #
  # ========================== read parmtables
  #
  mep_timestamp = os.path.getmtime(args[0]);

  if len(args) > 1 or not args[0].endswith('.cache'):
    from Owlcat.ParmTables import ParmTab

    # get antenna positions
    if options.ampl_slope or options.phase_slope or \
      options.phase_slope_rot_offset or options.phase_slope_time_offset or options.phase_slope_dlm:
      ms,msant = load_ms();
      names = list(msant.getcol('NAME'));
      pos = msant.getcol('POSITION');
      # get reference position
      x,y,z = pos0 = pos[0,:];
      # convert to unit vector pointing east
      dir_east = numpy.array([-y,x,0]);
      dir_east /= math.sqrt((dir_east**2).sum());
      # get coordinate of each antenna along dir_east: this is just the scalar product of
      # (pos[i,:]-pos[0,:]) by dir_east.
      antx = ((pos-pos0[numpy.newaxis,:])*dir_east[numpy.newaxis,:]).sum(1);
      # now convert to dict of antenna coordinates
      ANTX_dict = dict(zip(names,antx));
      # update with abbreviated antenna names, if they share a common prefix
      prefix = 1;
      while prefix < len(names[0]) and all([n[:prefix] == names[0][:prefix] for n in names[1:]]):
        prefix += 1;
      if prefix > 1:
        ANTX_dict.update([(n[prefix-1:],ax) for n,ax in zip(names,antx)]);
      # print "Antenna positions: "," ".join([ "%s:%f"%(n,ANTX_dict[n]) for n in sorted(ANTX_dict.keys()) ]);
    # else antenna positions not needed
    else:
      ANTX_dict = ANTX = None;

    # spectral windows
    SPWS = range(len(args));
    # set of all sources, antennas and correlations
    SRCS = set();
    ANTS = set();
    CORRS = set();

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

    # check that antenna positions are known
    if ANTX_dict:
      unknown_antennas = [ p for p in ANTS if p not in ANTX_dict ];
      if unknown_antennas:
        print "Don't have positions for antenna(s) %s"%",".join(unknown_antennas);
        print "Perhaps you should specify a proper MS with --ms?";
        sys.exit(1);

      # make vector of antenna positions
      ANTX = numpy.array([ ANTX_dict[p] for p in ANTS ]);
      # recenter at middle of the array
      ANTX -= ANTX[-1]/2;

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
      #mask = numpy.zeros((len(SRCS),len(ANTS),len(CORRS),NTIMES),dtype=bool);
      arr = numpy.zeros((len(SRCS),len(ANTS),len(CORRS),NTIMES),dtype=complex);
      #arr = numpy.ma.masked_array(mask,arr,fill_value=1.0);
      for i,src in enumerate(SRCS):
        for j,ant in enumerate(ANTS):
          for k,corr in enumerate(CORRS):
            fs_real = pt.funkset(':'.join([prefix,src,ant,corr,"r"])).get_slice();
            fs_imag = pt.funkset(':'.join([prefix,src,ant,corr,"i"])).get_slice();
            # table is allowed to contain either NO funklets for this src/ant/corr combination,
            # or the same number as all other tables.
            if len(fs_real) or len(fs_imag):
              if len(fs_real) != len(fs_imag) or len(fs_real) != NTIMES:
                print "Error: table contains %d real and %d imaginary funklets; %d expected"%(len(fs_real),len(fs_imag),NTIMES);
                sys.exit(1);
              arr[i,j,k,:] = numpy.array([complex(c00(fr),c00(fi)) for fr,fi in zip(fs_real,fs_imag)]);
            # no funklets -- use 1, and mask the array
            else:
              arr[i,j,k,:] = numpy.ones((NTIMES,),dtype=complex);
  #            arr.mask[i,j,k,:] = True;
      return arr;

    # stuff down the line does not handle masked arrays properly, so we don't use them
    #mask = numpy.zeros((len(SPWS),len(SRCS),len(ANTS)+1,len(CORRS),NTIMES),dtype=bool);
    de = numpy.zeros((len(SPWS),len(SRCS),len(ANTS)+1,len(CORRS),NTIMES),dtype=complex);
    #de = numpy.ma.masked_array(de,mask,fill_value=1.0);

    for spw,tabname in enumerate(args):
      print "Reading",tabname;
      pt = ParmTab(tabname);
      mep_timestamp = max(mep_timestamp,os.path.getmtime(tabname));
      freq0[spw] = sum(pt.envelope_domain().freq)/2;
      de[spw,:,0:len(ANTS),:,:] = read_complex_parms(pt);

    print "Read %d parmtables"%len(args);

    # write cache
    if options.cache:
      cachefile = options.cache+'.cache';
      cPickle.dump((freq0,de,SPWS,SRCS,ANTS,CORRS,NTIMES,ANTX),file(cachefile,'w'));
      print "Cached all structures to file",cachefile;

  #
  # ========================== read cache file
  #
  else:
    print "Reading cache file",args[0];
    freq0,de,SPWS,SRCS,ANTS,CORRS,NTIMES,ANTX = cPickle.load(file(args[0]));
    print "Read %s: %d spws, %d srcs, %d ants, %d corrs, %d times"%(args[0],
        len(SPWS),len(SRCS),len(ANTS),len(CORRS),NTIMES);

  # check options that specify antennas by name, to make sure antenna is known
  for ant in 'circle_ampl_ant','circle_phase_ant':
    antname = getattr(options,ant,None);
    if antname:
      try:
        globals()[ant] = ANTS.index(antname);
      except IndexError:
        parser.error("Antenna name '%s' not found."%antname);
    else:
      globals()[ant] = None;

  # read pointing errors cache
  if options.pe:
    pe_dlm,pe_bsz,pe_beam_sizes = cPickle.load(file(options.pe));

    # pe_dlm is 2 x NSPW x NANT x NTIME
    # pe_bsz is 2 x 2 x NSPW x NANT x NTIME
    if pe_dlm.shape != (2,len(SPWS),len(ANTS),NTIMES) or \
      pe_bsz.shape != (2,2,len(SPWS),len(ANTS),NTIMES):
      print "Shape of cached arrays in file %s does not match: %s, %s"%(options.pe,pe_dlm.shape,pe_bsz_shape);
      sys.exit(1);
    print "Read pointing errors and %d beam extent(s) from %s"%(pe_beam_sizes,options.pe);
    if pe_beam_sizes > 1:
      print "WARNING: more than 1 beam extent is not currently supported";

  # de is an NSPW x NSRC x NANT x NCORR x NTIME array of complex dEs
  # (NANT+1 actually, the last one is the mean-across-antennas value)
  #

  # some more constants
  MEAN = len(ANTS);         # index of "mean" antenna
  XX = CORRS.index('xx');
  YY = CORRS.index('yy');
  # set of all source,antenna,corr combinations
  ALL = [ (src,ant,corr) for src in SRCS for ant in ANTS for corr in CORRS ];

  # print per-source mean de_xx and de_yy amplitudes
  # (skip the last antenna since it's the MEAN one.)
  dex = abs(de[:,:,:-1,XX,:]).mean(2).mean(2).mean(0);
  dey = abs(de[:,:,:-1,YY,:]).mean(2).mean(2).mean(0);
  print "=== Mean |dExx|, |dEyy|:";
  for isrc,src in enumerate(SRCS):
    print "%8s=[%16.8f,%16.8f],"%(src,dex[isrc],dey[isrc]);

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
  # Then we take the mean across correlations. dep0 becomes NSPW x NSRC x NANT x NTIME.
  dep0 = dep0.mean(3);
  # cut out the last ("mean") antenna
  dep0 = dep0[:,:,0:-1,:];

  if ANTX is not None:
    # now fit phase offset and slope over array
    # dep0_b0/1 will be NSPWxNSRCxNCORRxNTIME
    sx = ANTX.sum()
    sx2 = (ANTX**2).sum();
    sy = dep0.sum(2)
    sxy = (dep0*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis]).sum(2);
    n = len(ANTX);
    ## slope of linear fit
    dep0_b1 = (n*sxy - sy*sx)/(n*sx2-sx**2);
    ## offset of linear fit
    dep0_b0 = sy/n - dep0_b1*(sx/n);
    ## stddev w.r.t. linear fit
    #deph_fit_std = (deph - deph_b1*fq[:,numpy.newaxis,numpy.newaxis,numpy.newaxis] - deph_b0).std(0);

    # now also compute residuals w.r.t. the slope fits, NSPW x NSRC x NANT x NTIME
    dep0_fitted = dep0_b0[:,:,numpy.newaxis,:] + \
      dep0_b1[:,:,numpy.newaxis,:]*ANTX[numpy.newaxis,numpy.newaxis,:,numpy.newaxis];
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

  # unnormalized average dE (for pointing error plots), NSPWxNSRCxNANTxNTIME
  de_unnorm = numpy.sqrt((abs(de[:,:,:,XX,:])**2+abs(de[:,:,:,YY,:])**2)/2);
  de_unnorm[:,:,MEAN,:] = numpy.mean(de_unnorm[:,:,0:len(ANTS),:],2);

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

  # read in source positions, if plotting mode requires it
  if model:
    maxabsdelog = math.sqrt(0.5); # abs(delog).max();
    maxsize = options.circle_maxsize;
    minsize = options.circle_minsize;
    lsrc = numpy.zeros(len(SRCS),float);
    msrc = numpy.zeros(len(SRCS),float);
    nsrc = numpy.zeros(len(SRCS),float);
    beam_lm_src = numpy.zeros((2,len(SRCS)),float);
    for i,name in enumerate(SRCS):
      src = model.findSource(name);
      lsrc[i],msrc[i],nsrc[i] = Coordinates.radec_to_lmn(src.pos.ra,src.pos.dec,*radec0);
      if hasattr(src,"beam_lm"):
        beam_lm_src[:,i] = src.beam_lm;
      else:
        beam_lm_src[:,i] = lsrc[i],msrc[i];

  # figure out optimal plot sizes
  if options.width or options.height:
    sz_as = sz_sa = (options.width,options.height);
  # figure out automatically depending on whether we have more sources or antennas
  else:
    # more sources than antennas
    if len(SRCS) < (len(ANTS)+1):
      sz_as = ( 210,min(290,30*len(ANTS)));
      sz_sa = ( 290,min(210,30*len(SRCS)));
    else:
      sz_as = ( 290,min(210,30*len(ANTS)));
      sz_sa = ( 210,min(290,30*len(SRCS)));

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

  #
  #============================ MAKE_FIGURE() plotting function
  #
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
        save = "%s_%s"%(options.output_prefix,save);
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
      plt.text(.5,0,col,fontsize=options.subtitle_fontsize,
        horizontalalignment='center',verticalalignment='bottom');
    iplot += 1;
    # now make plots
    for irow,row in rows:
      for icol,col in cols:
        iplot += 1;
        plt = fig.add_subplot(len(rows)+1,len(cols)+1,iplot);
        # plt.axis("off");
        if ylock and ylock != "col" and icol:
          plt.set_yticklabels([]);
        else:
          plt.yaxis.set_major_locator(MaxNLocator(4));
        plt.set_xticks([]);
        data = datafunc(irow,icol);
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
            if x == range(len(y)):
              plt.set_xlim(-1,len(y));
        elif mode is PLOT_ERRORBARS:
          y,yerr = data;
          plt.errorbar(range(len(y)) if xaxis is None else xaxis,y,yerr,fmt=None,capsize=1);
          if xaxis is None:
            plt.set_xlim(-1,len(y));
        if ylock:
          plt.set_ylim(ymin[irow,icol],ymax[irow,icol]);
        if hline is not None:
          plt.axhline(y=hline,color='black');
        for lab in plt.get_yticklabels():
          lab.set_fontsize(options.axis_fontsize);
      # add row label
      iplot += 1;
      plt = fig.add_subplot(len(rows)+1,len(cols)+1,iplot);
      plt.text(0,.5,row,fontsize='x-small',horizontalalignment='left',verticalalignment='center');
      plt.axis("off");
    # adjust layout
    fig.subplots_adjust(left=borders[0],right=borders[1],top=borders[3],bottom=borders[2],
      wspace=options.subplot_wspace or None);
    # plot if asked to
    if suptitle and options.title_fontsize:
      fig.suptitle(suptitle,fontsize=options.title_fontsize);
    if save:
      if options.portrait:
        orientation = 'portrait';
      elif options.landscape:
        orientation = 'landscape';
      else:
        orientation='portrait' if figsize[0]<figsize[1] else 'landscape';
      fig.savefig(save,papertype=options.papertype,
                  orientation=orientation);
      print "Wrote",save,"in",orientation,"orientation";
      fig = None;
      pyplot.close("all");
    return save;

  #
  #============================ MAKE_SKYMAP() plotting function
  #
  def make_skymap (ll,mm,
    markers,  # marker is a list of strings or dicts or tuples.
              # For each marker, if str: axes.plot(x,y,str) is called
              # if dict: axes.plot(x,y,**dict)
              # else tuple is func,arg,kwargs
              # and axes.func(*args,**kwargs) is called
    labels=None,
    radius=None,        # plot radius. If None, then it's set automatically
    suptitle=None,      # title of plot
    save=None,          # filename to save to
    format=None,        # format: use options.output_type by default
    figsize=(210,210),       # figure width,height in mm
    ):
    if save and (format or options.output_type.upper()) != 'X11':
      save = "%s.%s"%(save,format or options.output_type);
      if options.output_prefix:
        save = "%s_%s"%(options.output_prefix,save);
      # exit if figure already exists, and we're not refreshing
      if os.path.exists(save) and not options.refresh: # and os.path.getmtime(save) >= os.path.getmtime(__file__):
        print save,"exists, not redoing";
        return save;
    else:
      save = None;

    figsize = (figsize[0]/25.4,figsize[1]/25.4);
    fig = pyplot.figure(figsize=figsize,dpi=600);
    plt = fig.add_axes([circle_borders[0],circle_borders[2],circle_borders[1]-circle_borders[0],
                        circle_borders[3]-circle_borders[2]]);
    plt.axhline(y=0,color='black',linestyle=':');
    plt.axvline(x=0,color='black',linestyle=':');

    # if ll is None, use length of markers for everything
    if ll is None:
      ll = mm = numpy.zeros(len(markers),float);

    if not isinstance(markers,(list,tuple)):
      markers = [markers]*len(ll);

    if not labels:
      labels = [None]*len(ll);

    # set limits
    if options.radius:
      lim = options.radius;
    elif radius:
      lim = radius;
    else:
      lim = max(max(abs(ll)),max(abs(mm)))*1.05/ARCMIN;

    # plot markers
    for l,m,mark,label in zip(ll/ARCMIN,mm/ARCMIN,markers,labels):
      if isinstance(mark,str):
        plt.plot(l,m,mark);
      elif isinstance(mark,dict):
        plt.plot(l,m,**mark);
      elif isinstance(mark,tuple):
        func,args,kwargs = mark;
        getattr(plt,func)(*args,**kwargs);
      if label and options.circle_label_fontsize:
        plt.text(l,m-lim*options.circle_label_offset,
          label,fontsize=options.circle_label_fontsize,
          horizontalalignment='center',
          verticalalignment='top' if options.circle_label_offset else 'center');

    x = numpy.cos(numpy.arange(0,361)*math.pi/180);
    y = numpy.sin(numpy.arange(0,361)*math.pi/180);

    for R in options.grid_circle:
      plt.plot(x*R,y*R,linestyle=':',color='black');

    plt.set_xlim(lim,-lim);
    plt.set_ylim(-lim,lim);
    for lab in list(plt.get_xticklabels())+list(plt.get_yticklabels()):
      lab.set_fontsize(options.circle_axis_fontsize);

    if suptitle and options.title_fontsize:
      fig.suptitle(suptitle,fontsize=options.circle_title_fontsize);
    if save:
      if options.portrait:
        orientation = 'portrait';
      elif options.landscape:
        orientation = 'landscape';
      else:
        orientation='portrait' if figsize[0]<figsize[1] else 'landscape';
      fig.savefig(save,papertype=options.papertype,
                  orientation=orientation);
      print "Wrote",save,"in",orientation,"orientation";
      fig = None;
      pyplot.close("all");
    return save;

  #
  #============================ VARIOUS PER-CORRELATION PLOTS
  #
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

  #
  #============================ DLM OFFSETS
  #
  dlm_offsets = dphase_lmfit = None; # will be set if we fit them below
  if options.phase_slope_dlm:
    # compute best-fitting per-source offset
    # convert phase slope into phase offset on 0D baseline
    p0wl = dep0_b1*(ANTX[-1]-ANTX[0])*DEG;
    p0wl /= 2*math.pi*freq0[:,numpy.newaxis,numpy.newaxis]/299792458.
    # take mean across SPWs (since we've eliminated dependence on freq),
    # p0wl is now NSRCxNTIME
    p0wl = p0wl.mean(0);
    # uv needs to be NTIMEx2
    uv = uvw[:,:2];
    # dlm will be NSRCx2
    dlm = dlm_offsets = numpy.zeros((len(SRCS),2),float);
    print "=== Fitted dl,dm offsets:";
    for isrc,src in enumerate(SRCS):
      # flip sign of phase offset, since phase is -(ul+vm+w(n-1))
      x,res,rank,sing = linalg.lstsq(uv,-p0wl[isrc,:]);
      print "%8s=[%16.8g,%16.8g],     # %.2f\", %.2f\""%(src,x[0],x[1],x[0]/ARCMIN*60,x[1]/ARCMIN*60);
      dlm[isrc,:] = x;
    # shape is now 2xNSRC
    dlm = dlm.transpose();
    # now make fitted phase curves
    ulvmwn = -(uvw[:,:2,numpy.newaxis]*dlm[numpy.newaxis,:,:]).sum(1);
    dphase_lmfit = (ulvmwn*2*math.pi*freq0.mean()/299792458.)/(DEG*(ANTX[-1]-ANTX[0]));

    # make figure
    markers = [];
    labels = [];
    # set dlmscale so that longet arrow is .2 of max radius
    maxlm = math.sqrt((lsrc**2+msrc**2).max());
    maxdlm = math.sqrt((dlm_offsets**2).sum(0).max());
    dlmscale = options.dlm_scale or 0.2*maxlm/maxdlm;
    # make l/m and dl/dm array (in arcmin)
    l0m0 = numpy.array([lsrc,msrc]).copy()/ARCMIN;
    dldm = dlm_offsets.transpose().copy()*dlmscale/ARCMIN;
    radius = options.radius; max(math.sqrt((l0m0**2).sum(0).max()),math.sqrt(((l0m0+dldm)**2).sum(0).max()))*1.01;
    # loop over sources
    for isrc,name in enumerate(SRCS):
      l0,m0 = l0m0[:,isrc];
      r = (math.sqrt((dlm_offsets[isrc,:]**2).sum())/DEG)*3600;
      mark = ('arrow',list(l0m0[:,isrc])+list(dldm[:,isrc]),dict(
        width=1,head_width=2,head_length=2,
      ));
      markers.append(mark);
      markers.append(('text',(l0,m0,"%.2f\""%r),dict(
          fontsize=options.circle_label_fontsize,
          horizontalalignment='left' if dldm[0,isrc]>0 else 'right',
          verticalalignment='top' if dldm[1,isrc]>0 else 'bottom')));

    make_skymap(None,None,markers,labels=labels,radius=radius,
      figsize=(210,210),suptitle="Fitted l/m offsets",save="dE_lm_offsets");

  #
  #============================ PHASE SLOPE PLOTS
  #
  # dep0_b0/1 will be NSPWxNSRCxNTIME
  if options.phase_slope:
    # compute rotational offset curves, if needed
    dr = options.phase_slope_rot_offset;
    dt = options.phase_slope_time_offset;
    if dr or dt:
      # dr: rotate lm by fixed angle
      if dr:
        # compute l1,m1,n1: rotated source lmn
        dr *= DEG;
        cr,sr = math.cos(dr),math.sin(dr);
        l1 = lsrc*cr - msrc*sr;
        m1 = lsrc*sr + msrc*cr;
        n1 = numpy.sqrt(1-l1**2-m1**2);
        # compute dlmn: delta-lmn, shape is 3xNSRC
        dlmn = numpy.array([lsrc-l1,msrc-m1,nsrc-n1]);
        # work out the corresponding phase difference (ul+vm+wn term)
        # shape is NTIMEx3xNSRC, then sum second axis, giving NTIMExNSRC
        ulvmwn = (uvw[:,:,numpy.newaxis]*dlmn[numpy.newaxis,:,:]).sum(1);
      # dt: introduce a time lag into the UVWs (in fractions of a timeslot)
      if dt:
        duvw = (uvw0[31::60]-uvw)*dt;
        lmn = numpy.array([lsrc,msrc,nsrc-1]);
        # work out the corresponding phase difference (ul+vm+wn term)
        # shape is NTIMEx3xNSRC, then sum second axis, giving NTIMExNSRC
        ulvmwn = (duvw[:,:,numpy.newaxis]*lmn[numpy.newaxis,:,:]).sum(1);
      # now multiply (ul+vm+wn) by 2*pi/c
      dphase = ulvmwn*2*math.pi/299792458.
      # now, convert this into a phase offset per NTIMExNSRC
      dphase = dphase[:,:]*freq0.mean();
      # and to deg/m
      dphase = dphase/(DEG*(ANTX[-1]-ANTX[0]));
      # if we've fitted a dphase corresponding to dlm, add it here
      if dphase_lmfit is None:
        plotfunc = lambda i,isrc:(
            (None,numpy.mean(dep0_b1[:,isrc,:]*1000,0),numpy.std(dep0_b1[:,isrc,:]*1000,0)),
            dphase[:,isrc]*1000
          );
      else:
        plotfunc = lambda i,isrc:(
            (None,numpy.mean(dep0_b1[:,isrc,:]*1000,0),numpy.std(dep0_b1[:,isrc,:]*1000,0)),
            dphase_lmfit[:,isrc]*1000,dphase[:,isrc]*1000
          );
      # make the plot
      make_figure(enumerate(("",)),enumerate(SRCS),plotfunc,
        ylock=False,mode=PLOT_MULTI,
        suptitle="Fitted differential phase slope over array",
        save="dEphase_array_slopes");
    else:
      if dphase_lmfit is None:
        plotfunc = lambda i,isrc:(
            (None,numpy.mean(dep0_b1[:,isrc,:]*1000,0),numpy.std(dep0_b1[:,isrc,:],0)*1000),);
      else:
        plotfunc = lambda i,isrc:(
            (None,numpy.mean(dep0_b1[:,isrc,:]*1000,0),numpy.std(dep0_b1[:,isrc,:]*1000,0)),
            dphase_lmfit[:,isrc]*1000
          );
      make_figure(enumerate(("",)),enumerate(SRCS),plotfunc,
        ylock=False,mode=PLOT_MULTI,
        suptitle="Fitted differential phase slope over array",
        save="dEphase_array_slopes");

  #
  #============================ RESIDUALS W.R.T. PHASE SLOPE
  #
  if options.phase_slope_res:
    make_figure(enumerate(SRCS),enumerate(ANTS[:-1]),
      lambda isrc,iant:(dep0_res[:,isrc,iant,:].mean(0),
                        dep0_res[:,isrc,iant,:].std(0)),
      hline=0,ylock="row",mode=PLOT_ERRORBARS,figsize=sz_sa,
      suptitle="dE residual phase (w.r.t. slope fit) mean & stddev across all bands",
      save="dEphase_array_slopes_res");

  #
  #============================ AMPLITUDE SLOPE PLOTS
  #
  if options.ampl_slope:
    dea0s = [ dea0_b1[:,:,XX,:],dea0_b1[:,:,YY,:],dea0_b0[:,:,XX,:],dea0_b0[:,:,YY,:] ];
    make_figure(enumerate(("X slope 1/m","Y slope 1/m","X offset","Y offset")),enumerate(SRCS),
      lambda i,isrc:(numpy.mean(dea0s[i][:,isrc,:],0),
                    numpy.std(dea0s[i][:,isrc,:],0)),
      ylock=False,mode=PLOT_ERRORBARS,
      suptitle="Fitted differential ampl slope over array",
      save="dEampl_array_slopes");

  #
  #============================ POLARIZATION PLOT
  #
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


  #
  #============================ dE PHASE PLOTS
  #
  if options.phase:
    de_mean_phase = de_phase.mean(3);
    make_figure(enumerate(SRCS),enumerate(ANTS[:-1]),  # omit MEAN antenna since it's 1 by definition
      lambda isrc,iant:(numpy.mean(de_mean_phase[:,isrc,iant,:],0),
                        numpy.std(de_mean_phase[:,isrc,iant,:],0)),
      hline=0,ylock="row",figsize=sz_sa,mode=PLOT_ERRORBARS,
      suptitle="dE phase (deg) mean & stddev across all bands",
      save="dEphase_mean");

  #
  #============================ dE AMPLITUDE PLOTS
  #
  if options.ampl:
    if options.pe:
      de_pe = numpy.zeros(de_unnorm.shape,float);
      de_peshape = numpy.zeros(de_unnorm.shape,float);
      # pe_dlm is 2 x NSPW x NANT x NTIME
      # beam_lm_src is 2 x NSRC
      # compute theoretical gain [NSRC]
      r0 = numpy.sqrt((beam_lm_src**2).sum(0));
      # gain0 is NSRC x NSPW
      gain0 = r0[:,numpy.newaxis]*freq0[numpy.newaxis,:]*6.5e-8;
      gain0 = numpy.maximum(numpy.cos(gain0)**3,0.1);
      # lm is then the actual source offset: this is 2 x NSRC x NSPW x NANT x NTIME
#      lm = numpy.zeros((2,len(SRCS),len(ANTS),NTIMES),float);
      lm = beam_lm_src[:,:,numpy.newaxis,numpy.newaxis,numpy.newaxis] - \
        pe_dlm[:,numpy.newaxis,:,:,:];
      # convert to r=sqrt(l^2+m^2)*freq*C: NSRC x NSPW x NANT x NTIME
      pe_r = numpy.sqrt((lm**2).sum(0)); ## sqrt(l^2+m^2)
      # multiply by C and freq
      pe_r *= freq0[numpy.newaxis,:,numpy.newaxis,numpy.newaxis]*6.5e-8;
      # convert to gain
      beamgain_pe = numpy.maximum(numpy.cos(pe_r)**3,.1);
      beamgain_pe = beamgain_pe / gain0[:,:,numpy.newaxis,numpy.newaxis];
      # if we have shape, compute also beamgain due to shape
      if pe_beam_sizes > 0:
        beamgain_pe_shape = numpy.maximum(numpy.cos(pe_r/pe_bsz[0,0,numpy.newaxis,:,:,:])**3,.1);
        beamgain_pe_shape = beamgain_pe_shape / gain0[:,:,numpy.newaxis,numpy.newaxis]
        # make plotting func
        plotfunc = lambda isrc,iant:(
            (None,numpy.mean(de_unnorm[:,isrc,iant,:],0),numpy.std(de_unnorm[:,isrc,iant,:],0)),
            beamgain_pe[isrc,:,iant,:].mean(0),
            beamgain_pe_shape[isrc,:,iant,:].mean(0),
          );
      else:
        plotfunc = lambda isrc,iant:(
            (None,numpy.mean(de_unnorm[:,isrc,iant,:],0),numpy.std(de_unnorm[:,isrc,iant,:],0)),
            beamgain_pe[isrc,:,iant,:].mean(0),
          );
      make_figure(enumerate(SRCS),enumerate(ANTS[:-1]),
        plotfunc,
        hline=1,ylock="row",figsize=sz_sa,mode=PLOT_MULTI,
        suptitle="||dE|| mean & stddev across all bands, plus fitted ||dE|| due to pointing error",
        save="dE_mean");
    else:
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

  #
  #============================ ROGUES GALLERY AMPLITUDES
  #
  output_prefix = options.output_prefix+"_" if options.output_prefix else "";
  if options.circle_ampl:
    # de_renorm is NSPWxNSRCxNANTxNTIME
    # reshape as NSRCxNANTxNSPWxNTIME
    de_r1 = de_renorm.copy().transpose((1,2,0,3));
    demean = de_r1.mean(2).mean(2);
    # stddev is sqrt(s1^2+s2^2), where s1 is mean error, and s2 is error of mean
    destd  = de_r1.std(2).mean(2)/math.sqrt(NTIMES-1);
#    destd  = numpy.sqrt((de_r1.std(2).mean(2)/math.sqrt(NTIMES-1))**2 + de_r1.mean(2).std(2)**2);
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
        else:
          sigmas = max(sigmas,1);
          mark['markeredgewidth'] = options.circle_linewidth*((sigmas-1)/2.)+1;
#        elif sigmas >= 1:
#          mark['markeredgewidth'] = 2*options.circle_linewidth;
#        else:
#          mark['markeredgewidth'] = 1;
        markers.append(mark);
        if options.names_only:
          labels.append(name);
        else:
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
        outname = "%sdE_ant_gallery.png"%output_prefix;
        img.save(outname,"PNG");
        print "Wrote",outname;


  #
  #============================ ROGUES GALLERY ANIMATION
  #
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
        color = "#2b60d3" if log >=0 else "red";
        mark = dict(marker='o',markersize=minsize + (maxsize-minsize)*(abs(log)/maxabsdelog),
                    markeredgecolor=color,markerfacecolor=color);
        markers.append(mark);
        if options.names_only:
          labels.append(name);
        else:
          labels.append("%s %.2f+-%.2f"%(name,mean,std));

      frames.append(make_skymap(lsrc,msrc,markers,labels=labels,
        figsize=(210,210),suptitle="Average ||dE||, antenna %s, time slice %d"%(antname,k),save="%sdE_ant%s_%03d"%(output_prefix,antname,k)));

    # try to make animation
    if options.output_type.upper() == "PNG":
      for f in frames:
        os.system("pngtopnm %s | ppmquant 256 | ppmtogif >%s.gif"%(f,f));
      outname = "%sdE_ant%s_anim.gif"%(output_prefix,antname);
      os.system("gifsicle -d 15 -l0 --colors 256 %s >%s"%
                (" ".join([f+".gif" for f in frames + frames[-1::-1]]),
                  outname));
      os.system("rm -f %s"%" ".join([f+".gif" for f in frames]));
      print "Wrote",outname;

  #
  #============================ UPDATE LSM
  #
  if options.update_lsm:
    if lsm_timestamp >= mep_timestamp:
      if raw_input("Your LSM file appears to be newer than the dE solutions. Really update (y/n)? "
          ).lower()[:1] != "y":
        print "LSM update cancelled.";
        sys.exit(1);
    # update sources
    print "=== Updating model sources";
    A = (dex**2+dey**2)/2;
    B = (dex**2-dey**2)/2;
    for isrc,name in enumerate(SRCS):
      for src in model.sources:
        if src.name.startswith(name):
          if hasattr(src.flux,'Q'):
            I,Q = src.flux.I,src.flux.Q;
            src.flux.I = A[isrc]*I+B[isrc]*Q;
            src.flux.Q = A[isrc]*Q+B[isrc]*I;
            print "%8s I=%f Q=%f --> I=%f Q=%f"%(src.name,I,Q,src.flux.I,src.flux.Q);
          if dlm_offsets is not None:
            ra,dec = src.pos.ra,src.pos.dec;
            l,m,n = Coordinates.radec_to_lmn(ra,dec,*radec0);
            l += dlm_offsets[isrc,0];
            m += dlm_offsets[isrc,1];
            src.pos.ra,src.pos.dec = Coordinates.lm_to_radec(l,m,*radec0);
            print "%8s position %.8f,%.8f --> %.8f %.8f"%(src.name,ra,dec,src.pos.ra,src.pos.dec);
    newname = "updated-"+lsm_filename;
    model.save(newname);
    print "Wrote updated sky model",newname;

  if options.output_type.upper() == "X11":
    from pylab import plt
    plt.show();
