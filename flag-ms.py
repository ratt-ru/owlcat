#!/usr/bin/python
# -*- coding: utf-8 -*-

import os.path
import math
import sys
import re
import traceback

flagger = parser = ms = msname = None;

def error (message):
  print "%s: %s"%(os.path.basename(sys.argv[0]),message);
  sys.exit(1);

def get_ms ():
  global ms;
  global msname;
  if not ms:
    ms = Owlcat.table(msname);
  return ms;

def parse_subset_options (options):
  global parser;
  global flagger;
  global msname;
  global ms;
  subset = {};
  from Owlcat import Parsing

  # DDID and FIELD_ID
  if options.ddid is not None:
    try:
      subset['ddid'] = map(int,options.ddid.split(","));
      print "  ===> DATA_DESC_ID:",subset['ddid'];
    except:
      parser.error("Invalid -D/--ddid option");
  if options.field is not None:
    try:
      subset['fieldid'] = map(int,options.field.split(","));
      print "  ===> FIELD_ID:",subset['fieldid'];
    except:
      parser.error("Invalid -F/--field option");
  # taql
  taqls = [];
  if options.taql:
    taqls.append(options.taql);
    print "  ===> TaQL query:",options.taql;
  # channels
  if options.channels:
    subset['channels'] = map(Parsing.parse_slice,options.channels.split(",")); 
    print "  ===> channels:",subset['channels'];
  # corr list
  if options.corrs is not None:
    try:
      subset['corrs'] = map(int,options.corrs.split(','));
      print "  ===> correlations:",subset['corrs'];
    except:
      parser.error("Invalid -X/--corrs option");
  # station list
  if options.stations is not None:
    try:
      subset['antennas'] = map(int,options.stations.split(','));
      print "  ===> stations:",subset['antennas'];
    except:
      parser.error("Invalid -S/--stations option");
  # IFR set
  if options.ifrs is not None:
    import Meow.IfrSet
    ifrset = Meow.IfrSet.from_ms(get_ms());
    # print help and exit
    if options.ifrs == "help":
      # print help string, but trim away RTF tags
      print re.sub("<[^>]+>","",ifrset.subset_doc).replace("&lt;","<").replace("&gt;",">");
      sys.exit(0);
    try:
      ifrset = ifrset.subset(options.ifrs);
      print "  ===> ifrs:"," ".join([ifrset.ifr_label(ip,iq) for (ip,p),(iq,q) in ifrset.ifr_index()]);
    except:
      parser.error("Invalid -I/--ifrs option");
    taqls.append(ifrset.taql_string());
  # clipping
  subset['data_column'] = options.data_column;
  if options.above is not None:
    subset['data_above'] = options.above;
    print "  ===> select |%s|>%f"%(options.data_column,options.above);
  if options.below is not None:
    subset['data_below'] = options.below;
    print "  ===> select |%s|<%f"%(options.data_column,options.below);
  if options.fm_above is not None:
    subset['data_fm_above'] = options.fm_above;
    print "  ===> select mean|%s|>%f"%(options.data_column,options.fm_above);
  if options.fm_below is not None:
    subset['data_fm_below'] = options.fm_below;
    print "  ===> select mean|%s|<%f"%(options.data_column,options.fm_below);
  # join taql queries
  if taqls:
    subset['taql'] = "( " + " ) && ( ".join(taqls) + " )";
  # fill flag args
  for opt in 'data_flagmask','flagmask','flagmask_all','flagmask_none':
    subset[opt] = getattr(options,opt);
  return subset;


if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser,OptionGroup
  parser = OptionParser(usage="""%prog: [actions] [options] MS""",
      description="Manipulates flags (bitflags and legacy FLAG/FLAG_ROW columns) in the MS. "
      "Use the selection options to narrow down a subset of the data, and use the action options "
      "to change flags within that subset. Without any action options, statistics on the current "
      "selection are printed -- this is useful as a preview of your intended action."
  );

  group = OptionGroup(parser,"Selection by subset");
  group.add_option("-L","--channels",type="string",
                    help="channel selection: single number or start:end[:step] to select channels start through end-1, "
                    "or start~end[:step] to select channels start through end, with an optional stepping.");
  group.add_option("-X","--corrs",type="string",
                    help="correlation selection. Use comma-separated list of correlation indices.");
  group.add_option("-S","--stations",type="string",
                    help="station (=antenna) selection. Use comma-separated list of station indices."),
  group.add_option("-I","--ifrs",type="string",
                    help="interferometer selection. Use \"-I help\" to get help on selecting ifrs.");
  group.add_option("-D","--ddid",type="string",
                    help="DATA_DESC_ID selection. Single number, or comma-separated list.");
  group.add_option("-F","--field",type="string",
                    help="FIELD_ID selection. Single number, or comma-separated list.");
  group.add_option("-Q","--taql",dest="taql",type="str",
                    help="additional TaQL selection to restrict subset.");
  parser.add_option_group(group);

  group = OptionGroup(parser,"Selection by data value");
  group.add_option("--above",metavar="X",type="float",
                    help="select on abs(data)>X");
  group.add_option("--below",metavar="X",type="float",
                    help="select on abs(data)<X");
  group.add_option("--fm-above",metavar="X",type="float",
                    help="select on mean(abs(data))>X, where mean is over frequencies");
  group.add_option("--fm-below",metavar="X",type="float",
                    help="select on mean(abs(data))<X, where mean is over frequencies");
  group.add_option("--data-column",metavar="COLUMN",type="string",
                    help="data column for --above or --below options. Default is %default.");
  group.add_option("--data-flagmask",metavar="FLAGS",type="string",
                    help="flags to apply to data column (when e.g. computing mean). Default is %default. See below for an "
                    "poverview of what FLAGS can be set to.");
  parser.add_option_group(group);

  group = OptionGroup(parser,"Selection by current flags");
  group.add_option("-Y","--flagged-any",metavar="FLAGS",dest='flagmask',type="string",
                    help="selects if any of the specified flags are raised. For this and all other options taking "
                    "a FLAGS argument, FLAGS can be a flagset name or an integer bitmask "
                    "(if bitflags are in use -- see also the -l/--list option). Prefix the bitmask by '0x' to use hex. "
                    "Append a '+L' to include legacy boolean FLAG/FLAG_ROW columns. Use 'all' for "
                    "all bitflags, and 'ALL' for all bitflags plus legacy flags (equivalent to 'all+L'). FLAGS may "
                    "also be a comma-separated list of any of the above terms.");
  group.add_option("-A","--flagged-all",metavar="FLAGS",dest='flagmask_all',type="string",
                    help="selects if all of the specified flags are raised");
  group.add_option("-N","--flagged-none",metavar="FLAGS",dest='flagmask_none',type="string",
                    help="selects if none of the specified flags are raised");
  parser.add_option_group(group);

  group = OptionGroup(parser,"Actions to take on selection (may be combined)");
  group.add_option("-f","--flag",metavar="FLAGS",type="string",
                  help="raise the specified FLAGS");
  group.add_option("-u","--unflag",metavar="FLAGS",type="string",
                  help="clear the specified flags");
  group.add_option("-g","--fill-legacy",metavar="FLAGS",type="string",
                  help="fills legacy FLAG/FLAG_ROW columns using the specified FLAGS. When -f/--flag or -u/--unflag "
                  "or -r/--remove is used, legacy flags are implicitly reset using all bitflags: use '-g -' "
                  "to skip this step. You may also use this option on its own to reset legacy flags (within the "
                  "specified data subset) using some bitmask. Use '-g 0' to clear legacy flags.");
  group.add_option("-c","--create",action="store_true",
                  help="for -f/--flag option only: if a named flagset doesn't exist, creates "
                  "it. Without this option, an error is reported.");
  parser.add_option_group(group);

  group = OptionGroup(parser,"Other options");
  group.add_option("-l","--list",action="store_true", 
                  help="lists various info about the MS, including its flagsets.");
  group.add_option("-r","--remove",metavar="FLAGSET(s)",type="string",
                  help="unflags and removes named flagset(s). You can use a comma-separated list.");
  group.add_option("-v","--verbose",metavar="LEVEL",type="int",
                    help="verbosity level for messages. Higher is more verbose, default is 0.");
  parser.add_option_group(group);

  parser.set_defaults(data_column="CORRECTED_DATA",data_flagmask="ALL",
    flagged_any=None,flaged_all=None,flagged_none=None,
    flag=None,unflag=None,fill_legacy=None,
    verbose=0);

  # parse args
  (options,args) = parser.parse_args();
  if len(args) != 1:
    parser.error("Incorrect number of arguments. Use '-h' for help.");
  msname  = args[0];

  # import stuff
  import numpy
  import numpy.ma
  import Owlcat.Flagger
  from Owlcat.Flagger import Flagger

  # create flagger object
  flagger = Flagger(msname,verbose=options.verbose);

  #
  # -l/--list: list MS info
  #
  if options.list:
    ms = get_ms();
    ants = Owlcat.table(ms.getkeyword('ANTENNA')).getcol('NAME');
    ddid_tab = Owlcat.table(ms.getkeyword('DATA_DESCRIPTION'));
    spwids = ddid_tab.getcol('SPECTRAL_WINDOW_ID');
    polids = ddid_tab.getcol('POLARIZATION_ID');
    corrs = Owlcat.table(ms.getkeyword('POLARIZATION')).getcol('CORR_TYPE');
    spw_tab = Owlcat.table(ms.getkeyword('SPECTRAL_WINDOW'));
    ref_freq = spw_tab.getcol('REF_FREQUENCY');
    nchan = spw_tab.getcol('NUM_CHAN');
    fields = Owlcat.table(ms.getkeyword('FIELD')).getcol('NAME');

    print "===> MS is %s"%msname;
    print "  %d antennas: %s"%(len(ants)," ".join(ants));
    print "  %d DATA_DESC_ID(s): "%len(spwids);
    for i,(spw,pol) in enumerate(zip(spwids,polids)):
      print "    %d: %.3f MHz, %d chans x %d correlations"%(i,ref_freq[spw]*1e-6,nchan[spw],len(corrs[pol,:]));
    print "  %d field(s): %s"%(len(fields),", ".join(["%d: %s"%ff for ff in enumerate(fields)]));
    if not flagger.has_bitflags:
      print "No BITFLAG/BITFLAG_ROW columns in this MS. Use the 'addbitflagcol' command to add them.";
    else:
      names = flagger.flagsets.names();
      if names:
        print "  %d flagset(s): "%len(names);
        for name in names:
          mask = flagger.flagsets.flagmask(name);
          print "    '%s': %d (0x%02X)"%(name,mask,mask);
      else:
        print "  No flagsets.";
    print "";
    if options.flag or options.unflag or options.fill_legacy or options.remove:
      print "-l/--list was in effect, so all other options were ignored.";
    sys.exit(0);

  # --flag/--unflag/--remove implies '-g all' by default, '-g -' skips the fill-legacy step
  if options.flag or options.unflag or options.remove:
    if options.fill_legacy is None:
      options.fill_legacy = 'all';
    elif options.fill_legacy == '-':
      options.fill_legacy = None;

  # if no other actions supplied, enable stats
  if not (options.flag or options.unflag or options.fill_legacy):
    statonly = True;
  else:
    statonly = False;

  # convert all the various FLAGS to flagmasks (or Nones)
  for opt in 'data_flagmask','flagmask','flagmask_all','flagmask_none','flag','unflag','fill_legacy':
    value = getattr(options,opt);
    try:
      flagmask = flagger.lookup_flagmask(value,create=(opt=='flag' and options.create));
    except Exception,exc:
      msg = str(exc);
      if opt == 'flag' and not options.create:
        msg += "\nPerhaps you forgot the -c/--create option?";
      error(msg);
    setattr(options,opt,flagmask);

  # clear the legacy flag itself from fill_legacy, otherwise it can have no effect
  if options.fill_legacy is not None:
    options.fill_legacy &= ~Flagger.LEGACY;

  #
  # -r/--remove: remove flagsets 
  #
  if options.remove is not None:
    if options.flag or options.unflag:
      error("Can't combine -r/--remove with --f/--flag or -u/--unflag.");
    # get list of names and flagmask to remove
    remove_names = options.remove.split(",");
    names_found = [];
    names_not_found = [];
    flagmask = 0;
    for name in remove_names:
      try:
        flagmask |= flagger.flagsets.flagmask(name);
        names_found.append(name);
      except:
        names_not_found.append(name);
    if not names_found:
      if names_not_found:
        error("No such flagset(s): %s"%",".join(names_not_found));
      error("No flagsets to remove");
    # unflag flagsets
    if names_not_found:
      print "===> WARNING: flagset(s) not found, ignoring: %s"%",".join(names_not_found);
    print "===> removing flagset(s) %s"%",".join(names_found);
    print "===> and clearing corresponding flagmask %s"%Flagger.flagmaskstr(flagmask); 
    if options.fill_legacy is not None:
      print "===> and filling FLAG/FLAG_ROW using flagmask %s"%Flagger.flagmaskstr(options.fill_legacy);
    flagger.xflag(unflag=flagmask,fill_legacy=options.fill_legacy);
    flagger.flagsets.remove_flagset(*names_found);
    sys.exit(0);

  # parse subset options
  subset = parse_subset_options(options);

  # at this stage all remaining options are handled the same way
  flagstr = unflagstr = legacystr = None;
  if options.flag is not None:
    flagstr = flagger.flagmaskstr(options.flag);
    print "===> flagging with flagmask %s"%flagstr;
  if options.unflag is not None:
    unflagstr = flagger.flagmaskstr(options.unflag);
    print "===> unflagging with flagmask %s"%unflagstr; 
  if options.fill_legacy is not None:
    legacystr = flagger.flagmaskstr(options.fill_legacy);
    print "===> filling legacy flags with flagmask %s"%legacystr;

  # do the job
  totrows,sel_nrow,sel_nvis,nvis_A,nvis_B,nvis_C = \
    flagger.xflag(flag=options.flag,unflag=options.unflag,fill_legacy=options.fill_legacy,**subset);

  # print stats
  if statonly:
    print "===> No actions were performed. Showing the result of your selection:"
  else:
    print "===> Flagging stats:";
  rpc = 100.0/totrows;
  print "===>   MS size:               %8d rows"%totrows;
  print "===>   Data selection:        %8d rows, %8d visibilities (%.3g%% of MS rows)"%(sel_nrow,sel_nvis,sel_nrow*rpc);
  if legacystr:
    print "===>     (over which legacy flags were filled using flagmask %s)"%legacystr;

  percent = 100.0/sel_nvis;
  if options.channels or options.corrs:
    print "===>   Chan/corr slicing reduces this to     %8d visibilities (%.3g%% of selection)"%(nvis_A,nvis_A*percent);
  if not (options.flagmask is None and options.flagmask_all is None and options.flagmask_none is None):  
    print "===>   Flag selection reduces this to        %8d visibilities (%.3g%% of selection)"%(nvis_B,nvis_B*percent);
  if options.above is not None or options.below is not None or options.fm_above is not None or options.fm_below is not None:
    print "===>   Data clipping reduces this to         %8d visibilities (%.3g%% of selection)"%(nvis_C,nvis_C*percent);
  if unflagstr:
    print "===>     (which were unflagged using flagmask %s)"%unflagstr;
  if flagstr:
    print "===>     (which were flagged using flagmask %s)"%flagstr;

  flagger.close();