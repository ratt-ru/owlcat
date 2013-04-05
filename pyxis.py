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
import sys
import os.path
import re
import time
import subprocess
import imp

import Pyxis
import PyxisImpl.Internals

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser
  parser = OptionParser(usage="""%prog: [options] a=value command_x b=value c=value command_y ...""",
    description="Runs a sequence of reduction commands");
  parser.add_option("--log",type="string",metavar="FILENAME",
                    help="log output to file");
  parser.add_option("-l",dest="default_log",action="store_true",
                    help="equivalent to --log pyxis.log");

  (options,args) = parser.parse_args();

  if options.log:
    v.LOG = options.log;
  elif options.default_log:
    v.LOG = "pyxis.log";
    
  # if options
  # sort arguments into commands and MSs
  mslist = [];
  commands = []
  for arg in args:
    if not PyxisImpl.Internals._re_assign.match(arg) and re.match(".*\\.(MS|ms)$",arg):
      assign("MS",arg);
      mslist.append(arg);
    else:
      commands.append(arg);
      
  # MS list from command line overrides defaults
  if mslist:
    info("setting MS list with",mslist);
    globals().pop("MS_List_Template",None);
  
  # run commands
  if not commands:
    info("no commands to execute");
  else:
    run(*commands);
  