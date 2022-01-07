#!/bin/bash
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

if [ "$2" == "" ]; then
  echo "Usage: $0 msname ifr_subset_string"
  exit 2
fi

export PYTHONPATH=$PYTHONPATH:$MEQTREES_CATTERY_PATH:~/Cattery
#python -c "import Meow.IfrSet;ss=Meow.IfrSet.from_ms('$1').subset('$2'); print '%d//%s'%(len(ss.ifr_index()),ss.taql_string());" 
python3 -c "import Meow.IfrSet;ss=Meow.IfrSet.from_ms('$1').subset('$2'); print('%d//%s'%(len(ss.ifr_index()),ss.taql_string()));" 2>/dev/null | grep ANTENNA1==
