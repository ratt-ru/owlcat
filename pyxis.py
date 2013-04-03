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
from Pyxis.Context import *


if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser
  parser = OptionParser(usage="""%prog: [options] a=value command_x b=value c=value command_y ...""",
    description="Runs a sequence of reduction commands");
  parser.add_option("--log",type="string",
                    help="name of log file");
  parser.add_option("-f","--force",dest="force",action="store_true",
                    help="proceed without confirmation, and overwrite output MS if it already exists");
  parser.add_option("-s","--renumber-spws",dest="renumber",action="store_true",
                    help="treat each MS as a separate spectral window");

  (options,args) = parser.parse_args();

  # sort arguments into commands and MSs
  mslist = [];
  commands = []
  for arg in args:
    if re.match(".*\\.(MS|ms)$",arg):
      mslist.append(arg);
    else:
      commands.append(arg);
      
  # load config files -- all variable assignments go into the Context scope
  Pyxis.initconf();
  
  # MS list from command line overrides defaults
  if mslist:
    info("PYRATE: overriding MS list with",mslist);
    Pyxis.Context.MS_List = mslist;
    Pyxis.Context.pop("MS_List_Template",None);
  
  # run commands
  if not commands:
    Pyxis._info("PYRATE: no commands to execute");
  else:
    Pyxis.run(*commands);
  