#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import pyfits
import re
import os.path
import numpy

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
  parser.add_option("-m","--mean",dest="mean",action="store_true",
                    help="take mean of input images");
  parser.add_option("-d","--diff",dest="diff",action="store_true",
                    help="take difference of 2 input images");
  parser.add_option("-z","--zoom",dest="zoom",type="int",metavar="NPIX",
                    help="zoom into central region of NPIX x NPIX size");
  parser.add_option("-s","--scale",dest="scale",type="float",
                    help="rescale image values");

  parser.set_defaults(output="",mean=False,zoom=0,scale=1);

  (options,imagenames) = parser.parse_args();

  if not imagenames:
    parser.error("No images specified. Use '-h' for help.");

  print "%d input image(s): %s"%(len(imagenames),", ".join(imagenames));
  images = [ pyfits.open(img) for img in imagenames ];
  updated = False;

  if options.replace:
    if options.output:
      parser.error("Cannot combine -r/--replace with -o/--output");
    outname = imagenames[0];
    options.force = True;
  else:
    outname = options.output;
    autoname = not outname;
    if autoname:
      outname = re.split('[_]',imagenames[0],1)[-1];

  if options.sanitize is not None:
    print "Sanitizing: replacing INF/NAN with",options.sanitize;
    for img in images:
      d = img[0].data;
      d[numpy.isnan(d)+numpy.isinf(d)] = options.sanitize;
    updated = True;

  if options.diff and options.mean:
    parser.error("Cannot do both --mean and --diff at once.");
  
  if options.diff:
    if len(images) != 2:
      parser.error("The --diff option requires exactly two input images.");
    if autoname:
      outname = "diff_" + outname;
    print "Computing difference";
    data = images[0][0].data;
    data -= images[1][0].data;
    updated = True;

  if options.mean:
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
  
  if options.scale != 1:
    if autoname and not updated:
      outname = "rescale_" + outname;
    if len(images) > 1:
      "Too many input images specified for this operation, at most 1 expected";
      sys.exit(2);
    print "Applying scaling factor of %f to image values"%options.scale;
    images[0][0].data *= options.scale;
    updated = True;

  if updated:
    print "Writing output image",outname;
    if os.path.exists(outname) and not options.force:
      print "Output image exists, rerun with the -f switch to overwrite.";
      sys.exit(1);
    images[0].writeto(outname,clobber=True);
  else:
    print "No operations specified. Use --help for help."
