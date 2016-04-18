#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#% $Id$
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
import math
import sys
import re
import gzip
import traceback
import cPickle
import Owlcat.Tables

flagger = parser = ms = msname = None;

def error (message):
  sys.stderr.write("%s: %s\n"%(os.path.basename(sys.argv[0]),message));
  sys.exit(1);

def get_ms (readonly=True):
  global ms;
  global msname;
  if not ms:
    ms = Owlcat.table(msname,readonly=readonly);
  return ms;

def shape_str (label,arr):
  return "%s: %s"%(label,"x".join(map(str,arr.shape)) if arr is not None else "None");

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
      if not ifrset.ifrs():
        return None;
    except:
      parser.error("Invalid -I/--ifrs option");
    taqls.append(ifrset.taql_string());
  # clipping
  subset['data_column'] = options.data_column;
  if options.nan:
    subset['data_nan'] = True;
    print "  ===> select %s = NAN or INF"%options.data_column;
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
  group.add_option("-T","--timeslots",type="string",
                    help="timeslot selection: single number or start:end to select timeslots start through end-1, "
                    "or start~end to select timeslots start through end.");
  group.add_option("-M","--timeslot-multiplier",type="int",default=1,
                    help="multiplies the timeslot numbers given to -T by the given factor. Default is 1.");
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
  group.add_option("--nan",action="store_true",
                    help="select on invalid data (NaN or infinite)");
  group.add_option("--fm-above",metavar="X",type="float",
                    help="select on mean(abs(data))>X, where mean is over frequencies");
  group.add_option("--fm-below",metavar="X",type="float",
                    help="select on mean(abs(data))<X, where mean is over frequencies");
  group.add_option("-C","--data-column",metavar="COLUMN",type="string",
                    help="data column for --above/--below/--nan options. Default is %default.");
  group.add_option("--data-flagmask",metavar="FLAGS",type="string",
                    help="flags to apply to data column (when e.g. computing mean). Default is %default. See below for "
                    "details on specifying flags.");
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
  group.add_option("-x","--extend-all-corr",action="store_true",
                   help="apply selection to all correlations if at least one is selected");
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
  group.add_option("-s","--stats",action="store_true",
                  help="prints per-flagset flagging stats.");
  group.add_option("-r","--remove",metavar="FLAGSET(s)",type="string",
                  help="unflags and removes named flagset(s). You can use a comma-separated list.");
  group.add_option("--export",type="string",metavar="FILENAME",
                  help="exports all flags to flag file. FILENAME may end with .gz to produce a gzip-compressed file. If any flagging actions are specified, these will be done before the export." );
  group.add_option("--import",type="string",dest="_import",metavar="FILENAME",
                  help="imports flags from flag file. If any flagging actions are specified, these will be done after the import.");
  group.add_option("-v","--verbose",metavar="LEVEL",type="int",
                    help="verbosity level for messages. Higher is more verbose, default is 0.");
  group.add_option("--timestamps",action="store_true",
                  help="adds timestamps to verbosity messages.");
  group.add_option("-z","--chunk-size",metavar="NROWS",type="int",default=200000,
                    help="Number of rows to process at once. Default is %default. Set to higher values if you have RAM to spare.");
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

  import Owlcat

  # import flags from file, if so specified
  # these are the columns that are imported and exported
  FLAGCOLS = "FLAG","FLAG_ROW","BITFLAG","BITFLAG_ROW";
  if options._import:
    try:
      dump = Owlcat.Tables.TableDump(options._import,compress=True);
      ms = get_ms(readonly=False);
      print "Importing flags from %s:"%options._import;
      dump.load(ms,verbose=True);
      dump.close();
      ms.close();
    except:
      traceback.print_exc();
      error("Error importing flags from %s, exiting"%options._import);
    print "Flags imported OK.";

  # if no other actions supplied, enable stats (unless flags were imported, in which case just exit)
  if not (options.flag or options.unflag or options.fill_legacy):
    if options._import:
      sys.exit(0);
    statonly = True;
  else:
    statonly = False;

  import numpy
  import numpy.ma
  import Owlcat.Flagger
  from Owlcat.Flagger import Flagger

  # now, skip most of the actions below if we're in statonly mode and exporting
  if not (statonly and options.export):
    # create flagger object
    flagger = Flagger(msname,verbose=options.verbose,timestamps=options.timestamps,chunksize=options.chunk_size);

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

    # if no other actions supplied, enable stats (unless flags were imported, in which case just exit)
    if not (options.flag or options.unflag or options.fill_legacy):
      if options._import:
        sys.exit(0);
      statonly = not options.export;
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
      # get set of flagsets to remove
      remove_flagsets = set(options.remove.split(","))
      # get set of all flagsets
      all_flagsets = set(flagger.flagsets.names())
      # warn if any named flagsets were not found
      if remove_flagsets - all_flagsets:
        print "===> WARNING: flagset(s) %s not found, ignoring"%",".join(remove_flagsets - all_flagsets)
      # build flagmask of remaining flagsets
      retain = all_flagsets-remove_flagsets
      if not retain:
        #if names_not_found:
        # error("No such flagset(s): %s"%",".join(names_not_found));
        print "===> WARNING: no flagsets to remove, exiting"
        sys.exit(0)
      flagmask = 0
      for name in retain:
        flagmask |= flagger.flagsets.flagmask(name)
      print "===> removing flagset(s) %s"%",".join(all_flagsets-retain);
      print "===> and clearing flagmask %s"%Flagger.flagmaskstr(~flagmask);
      if options.fill_legacy is not None:
        print "===> and filling FLAG/FLAG_ROW using flagmask %s"%Flagger.flagmaskstr(options.fill_legacy);
      flagger.xflag(unflag=~flagmask,fill_legacy=options.fill_legacy);
      flagger.flagsets.remove_flagset(*list(all_flagsets-retain));
      sys.exit(0);

    # parse subset options
    subset = parse_subset_options(options);
    if not subset:
      print "===> ended up with empty subset, exiting";
      sys.exit(0);

    # convert timeslots to reltime option, if specified
    if options.timeslots:
      from Owlcat import Parsing
      tslice = Parsing.parse_slice(options.timeslots,options.timeslot_multiplier);
      times = sorted(set(get_ms().getcol('TIME')));
      time0 = times[0] if tslice.start is None else times[tslice.start];
      time1 = times[-1] if tslice.stop is None else times[tslice.stop-1];
      time0 -= times[0];
      time1 -= times[0];
      subset['reltime'] = time0,time1;
      print "  ===> select timeslots %s (reltime %g~%g s)"%(tslice,time0,time1);

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


    # if --stats in effect, loop over all flagsets and print stats
    if options.stats:
      print "===> --stats in effect, showing per-flagset statistics"
      printed_header = False;
      for flagset in list(flagger.flagsets.names())+["+L"]:
        # get stats for this flagset
        subset['flagmask'] = flagger.lookup_flagmask(flagset);
        totrows,sel_nrow,sel_nvis,nvis_A,nvis_B,nvis_C = flagger.xflag(**subset);
        percent = 100.0/sel_nvis if sel_nvis else 0;
        # print them
        if flagset is "+L":
          label = "legacy FLAG/FLAG_ROW";
        else:
          label =  "Flagset %s (0x%02X)"%(flagset,flagger.flagsets.flagmask(flagset));
        if not printed_header:
          printed_header = True;
          rpc = 100.0/totrows if totrows else 0;
          print "===>   MS size:               %8d rows"%totrows;
          print "===>   Data/time selection:   %8d rows, %10d visibilities (%.3g%% of MS rows)"%(sel_nrow,sel_nvis,sel_nrow*rpc);
          if options.channels or options.corrs:
            print "===>   Chan/corr slicing reduces this to    %12d visibilities (%.3g%% of selection)"%(nvis_A,nvis_A*percent);
        print "===>   %-29s includes %10d visibilities (%.3g%% of selection)"%(label,nvis_B,nvis_B*percent);
      sys.exit(0);

    # else not stats mode, do the actual flagging job
    totrows,sel_nrow,sel_nvis,nvis_A,nvis_B,nvis_C = \
      flagger.xflag(flag=options.flag,unflag=options.unflag,fill_legacy=options.fill_legacy,
        flag_allcorr=options.extend_all_corr,
        **subset);
      
    # print stats
    if statonly:
      print "===> No actions were performed. Showing the result of your selection:"
    else:
      print "===> Flagging stats:";
    rpc = 100.0/totrows if totrows else 0;
    print "===>   MS size:               %8d rows"%totrows;
    print "===>   Data/time selection:   %8d rows, %10d visibilities (%.3g%% of MS rows)"%(sel_nrow,sel_nvis,sel_nrow*rpc);
    if legacystr:
      print "===>     (over which legacy flags were filled using flagmask %s)"%legacystr;

    percent = 100.0/sel_nvis if sel_nvis else 0;
    if options.channels or options.corrs:
      print "===>   Chan/corr slicing reduces this to     %10d visibilities (%.3g%% of selection)"%(nvis_A,nvis_A*percent);
    if not (options.flagmask is None and options.flagmask_all is None and options.flagmask_none is None):
      print "===>   Flag selection reduces this to        %10d visibilities (%.3g%% of selection)"%(nvis_B,nvis_B*percent);
    if options.nan or options.above is not None or options.below is not None or \
        options.fm_above is not None or options.fm_below is not None:
      print "===>   Data selection reduces this to         %10d visibilities (%.3g%% of selection)"%(nvis_C,nvis_C*percent);
    if unflagstr:
      print "===>     (which were unflagged using flagmask %s)"%unflagstr;
    if flagstr:
      print "===>     (which were flagged using flagmask %s)"%flagstr;

    flagger.close();

  # export flags from file, if so specified
  if options.export:
    try:
      dump = Owlcat.Tables.TableDump(options.export,write=True,compress=True);
      ms = get_ms();
      colnames = set(ms.colnames());
      print "Exporting flags to %s, %d rows:"%(options.export,ms.nrows());
      for colname in FLAGCOLS:
        if colname in colnames:
          dump.dump_column(ms,colname,verbose=True);
      ms.close();
      dump.close();
    except:
      traceback.print_exc();
      error("Error exporting flags to %s"%options.export);
    print "Flags exported OK.";
