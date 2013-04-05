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
import re
import os.path
import numpy
## ugly hack to get around UGLY FSCKING ARROGNAT (misspelling fully intentional) pyfits-2.3 bug
import Kittens.utils
pyfits = Kittens.utils.import_pyfits();
import scipy.ndimage.measurements
import math

SANITIZE_DEFAULT = 12345e-7689

if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  from optparse import OptionParser
  parser = OptionParser(usage="""%prog: [options] <image names...>""");
  parser.add_option("-o","--output",dest="output",type="string",
                    help="name of output FITS file");
  parser.add_option("-r","--replace",action="store_true",
                    help="replace (first) input file by output. Implies '--force'.");
  parser.add_option("-f","--force",dest="force",action="store_true",
                    help="overwrite output file even if it exists");
  parser.add_option("-S","--sanitize",type="float",metavar="VALUE",
                    help="sanitize FITS files by replacing NANs and INFs with VALUE");
  parser.add_option("-N","--nonneg",action="store_true",
                    help="replace negative values by 0");
  parser.add_option("-m","--mean",dest="mean",action="store_true",
                    help="take mean of input images");
  parser.add_option("-d","--diff",dest="diff",action="store_true",
                    help="take difference of 2 input images");
  parser.add_option("-t","--transfer",action="store_true",
                    help="transfer data from image 2 into image 1, preserving the FITS header of image 1");
  parser.add_option("-z","--zoom",dest="zoom",type="int",metavar="NPIX",                    
                    help="zoom into central region of NPIX x NPIX size");
  parser.add_option("-R","--rescale",dest="rescale",type="float",
                    help="rescale image values");
  parser.add_option("-E","--edit-header",metavar="KEY=VALUE",type="string",action="append",
                    help="replace header KEY with VALUE. Use KEY=VALUE for floats and KEY='VALUE' for strings.");
  parser.add_option("-H","--header",action="store_true",help="print header(s) of input image(s)");
  parser.add_option("-s","--stats",action="store_true",help="print stats on images and exit. No output images will be written.");

  parser.set_defaults(output="",mean=False,zoom=0,rescale=1,edit_header=[]);

  (options,imagenames) = parser.parse_args();

  if not imagenames:
    parser.error("No images specified. Use '-h' for help.");

  # print "%d input image(s): %s"%(len(imagenames),", ".join(imagenames));
  images = [ pyfits.open(img) for img in imagenames ];
  updated = False;
  
  if options.header:
    for filename,img in zip(imagenames,images):
      if len(imagenames)>1:
        print "======== FITS header for",filename;
      for hdrline in img[0].header.ascard:
        print hdrline; 
  
  if options.replace or len(imagenames)<2:
    if options.output:
      parser.error("Cannot combine -r/--replace with -o/--output");
    outname = imagenames[0];
    options.force = True;
    autoname = False;
  else:
    outname = options.output;
    autoname = not outname;
    if autoname:
      outname = re.split('[_]',imagenames[0],1)[-1];

  for keyval in options.edit_header:
    key,val = keyval.split("=");
    if val[0] == "'" and val[-1] == "'":
      images[0][0].header[key] = val[1:-1:];
    elif val[-1] == 'd' or key.startswith('NAXIS'):
      images[0][0].header[key] = int(val[:-1] if val[-1]=='d' else val);
    else:
      images[0][0].header[key] = float(val);
    print "Setting header %s=%s"%(key,val);
    updated = True;

  if options.sanitize is not None:
    print "Sanitizing: replacing INF/NAN with",options.sanitize;
    for img in images:
      d = img[0].data;
      d[numpy.isnan(d)+numpy.isinf(d)] = options.sanitize;
    # if using stats, do not generate output
    if not options.stats:
      updated = True;

  if options.nonneg:
    print "Replacing negative value by 0";
    for img,name in zip(images,imagenames)[:1]:
      d = img[0].data;
      wh = d<0;
      d[wh] = 0;
      print "Image %s: replaced %d points"%(name,wh.sum());
    updated = True;

  if options.transfer:
    if len(images) != 2:
      parser.error("The --transfer option requires exactly two input images.");
    if autoname:
      outname = "xfer_" + outname;
    print "Transferring %s into coordinate system of %s"%(imagenames[1],imagenames[0]);
    images[0][0].data[...] = images[1][0].data;
    updated = True;
  elif options.diff:
    if len(images) != 2:
      parser.error("The --diff option requires exactly two input images.");
    if autoname:
      outname = "diff_" + outname;
    print "Computing difference";
    data = images[0][0].data;
    data -= images[1][0].data;
    updated = True;
  elif options.mean:
    if autoname:
      outname = "mean%d_"%len(images) + outname;
    print "Computing mean";
    data = images[0][0].data;
    for img in images[1:]:
      data += img[0].data;
    data /= len(images);
    images = [ images[0] ];
    updated = True;

  if options.zoom:
    z = options.zoom;
    if autoname:
      outname = "zoom%d_"%z + outname;
    if len(images) > 1:
      "Too many input images specified for this operation, at most 1 expected";
      sys.exit(2);
    data = images[0][0].data;
    nx = data.shape[2];
    ny = data.shape[3];
    zdata = data[:,:,(nx-z)/2:(nx+z)/2,(ny-z)/2:(ny+z)/2];
    print "Making zoomed image of shape","x".join(map(str,zdata.shape));
    images = [ pyfits.PrimaryHDU(zdata) ];
    updated = True;

  if options.rescale != 1:
    if autoname and not updated:
      outname = "rescale_" + outname;
    if len(images) > 1:
      "Too many input images specified for this operation, at most 1 expected";
      sys.exit(2);
    print "Applying scaling factor of %f to image values"%options.rescale;
    images[0][0].data *= options.rescale;
    updated = True;

  if updated:
    imagenames[0] = outname;

  if options.stats:
    for ff,filename in zip(images,imagenames):
      data = ff[0].data;
      min,max,dum1,dum2 = scipy.ndimage.measurements.extrema(data);
      sum = data.sum();
      mean = sum/data.size;
      std = math.sqrt(((data-mean)**2).mean());
      print "%s: min %g, max %g, sum %g, np %d, mean %g, std %g"%(filename,min,max,sum,data.size,mean,std); 
    sys.exit(0);
      
  if updated:
    print "Writing output image",outname;
    if os.path.exists(outname) and not options.force:
      print "Output image exists, rerun with the -f switch to overwrite.";
      sys.exit(1);
    images[0].writeto(outname,clobber=True);
  elif not options.header:
    print "No operations specified. Use --help for help."
