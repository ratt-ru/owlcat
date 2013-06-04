# -*- coding: utf-8 -*-
#
#% $Id: __init__.py 6998 2009-06-26 08:51:15Z cwilliams $ 
#
#
# Copyright (C) 2002-2007
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

import __main__
setattr(__main__,"_meow_verbosity",0);

## ugly hack to get around UGLY FSCKING ARROGNAT (misspelling fully intentional) pyfits-2.3 bug
## (just in case somebody imports pyfits somewhere)
import Kittens.utils
Kittens.utils.import_pyfits();

# import table class
try:
  from pyrap_tables import table,tablecopy,tableexists,tabledelete,addImagingColumns
except:
  try:
    from pyrap.tables import table,tablecopy,tableexists,tabledelete,addImagingColumns
  except:
    print "Failed to import pyrap_tables or pyrap.tables. Please install the pyrap "
    "package from http://code.google.com/p/pyrap/, or from a MeqTrees binary repository "
    "(see http://www.astron.nl/meqwiki/Downloading)"
    raise;

# pyrap likes to display a lot of these, so shut them off
import warnings
warnings.simplefilter('ignore',DeprecationWarning);

# list of optional packages which will be added to the include path
_Packages = [ "Cattery" ];

# list of locations where packages will be searched for
_PackageLocations = [ "~","~/Frameworks/",
  "/usr/local/MeqTrees","/usr/local/lib/MeqTrees","/usr/lib/MeqTrees","/usr/lib64/MeqTrees","/usr/lib32/MeqTrees"
  "/usr/local/meqtrees","/usr/local/lib/meqtrees","/usr/lib/meqtrees","/usr/lib64/meqtrees","/usr/lib32/meqtrees"
  ];

# mapping of package: path. Filled in as we find packages
_packages = {};

import sys
import os
import os.path

def packages ():
  """Returns mapping of available packages to their paths""";
  return _packages;
  # print "Using %s, set the %s_PATH environment variable to override this."%(path,package.upper());

def _tryPackageDir (path,package):
  """Tests if path refers to a valid directory, adds it to system include path if so.
  Marks package as having this path.""";
  if os.path.isdir(path):
    sys.path.insert(0,path);
    # check for version info
    try:
      version = ' '.join(file(os.path.join(path,'version_info')));
    except:
      version = 'no version info';
    # insert into packages
    global _packages;
    _packages[package] = path,version;
    return True;
  return False;

def _setPackagePath (package):
  """Finds the given package, by first looking in $MEQTREES_PACKAGE_PATH, then checking for 
  subdirectories of the standard _PackageLocations list.""";
  # check for explicit MEQTREES_PACKAGE_PATH first
  varname = 'MEQTREES_%s_PATH'%package.upper()
  path = os.environ.get(varname,None);
  if path:
    if not _tryPackageDir(path,package):
      print "Warning: your %s environment variable is set to"%varname;
      print "%s, but this is not a valid directory."%path;
      print "The %s package will not be available."%package;
    return;
  # else look in standard places
  for path in _PackageLocations:
    path = os.path.expanduser(path);
    if _tryPackageDir(os.path.join(path,package),package):
      return;
  # none found
  print "Warning: No %s package found."%package;
  print "If you have %s in a non-standard location, please set the %s environment"%(package,varname);
  print "variable to point to it."

for pkg in _Packages:
  _setPackagePath(pkg);
  
def find_exec (execname):
  import os
  path = os.environ.get('PATH') or os.defpath;
  for dirname in path.split(os.pathsep):
    fname = os.path.join(dirname,execname);
    if os.access(fname,os.R_OK|os.X_OK):
      return fname;
  return None;


