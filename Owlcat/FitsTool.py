# -*- coding: utf-8 -*-

#
# % $Id$
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
from __future__ import print_function
import sys
import re
import os.path
import numpy
from astropy.io import fits as pyfits
import scipy.ndimage.measurements
import math
#from astLib.astWCS import WCS
from astropy.wcs import WCS
import glob

SANITIZE_DEFAULT = 12345e-7689


def stack_planes(fitslist, outname='combined.fits', axis=0, ctype=None, keep_old=False, fits=False):
    """ Stack a list of fits files along a given axiis.
       
       fitslist: list of fits file to combine
       outname: output file name
       axis: axis along which to combine the files
       fits: If True will axis FITS ordering axes
       ctype: Axis label in the fits header (if given, axis will be ignored)
       keep_old: Keep component files after combining?
    """

    hdu = pyfits.open(fitslist[0])[0]
    hdr = hdu.header
    naxis = hdr['NAXIS']

    # find axis via CTYPE key
    if ctype is not None:
        for i in range(1, naxis + 1):
            if hdr['CTYPE%d' % i].lower().startswith(ctype.lower()):
                axis = naxis - i  # fits to numpy convention
    elif fits:
        axis = naxis - axis

    fits_ind = abs(axis - naxis)
    crval = hdr['CRVAL%d' % fits_ind]

    imslice = [slice(None)] * naxis
    _sorted = sorted([pyfits.open(fits) for fits in fitslist],
                     key=lambda a: a[0].header['CRVAL%d' % (naxis - axis)])

    # define structure of new FITS file
    nn = [hd[0].header['NAXIS%d' % (naxis - axis)] for hd in _sorted]
    shape = list(hdu.data.shape)
    shape[axis] = sum(nn)
    data = numpy.zeros(shape, dtype=float)

    for i, hdu0 in enumerate(_sorted):
        h = hdu0[0].header
        d = hdu0[0].data
        imslice[axis] = list(range(sum(nn[:i]), sum(nn[:i + 1])))
        data[imslice] = d
        if crval > h['CRVAL%d' % fits_ind]:
            crval = h['CRVAL%d' % fits_ind]
        try:
            hdr['BMAJ%d' % (i + 1)] = h['BMAJ']
            hdr['BMIN%d' % (i + 1)] = h['BMIN']
            hdr['BPA%d' % (i + 1)] = h['BPA']
        except KeyError:
            pass

    # update header
    hdr['CRVAL%d' % fits_ind] = crval
    hdr['CRPIX%d' % fits_ind] = 1

    pyfits.writeto(outname, data.astype(numpy.float32), hdr, overwrite=True)
    print(("Successfully stacked images. Output image is %s" % outname))

    # remove old files
    if not keep_old:
        for fits in fitslist:
            os.system('rm -f %s' % fits)


def unstack_planes(fitsname, each_chunk, axis=None, ctype=None, prefix=None, fits=False, keep_old=True):
    """ 
        Unstack FITS image along a given axis into multiple 
        images each having each_chunk planes along that axis 
    """

    prefix = prefix or fitsname[:-5]  # take everthing but .FITS/.fits
    hdu = pyfits.open(fitsname)
    hdr = hdu[0].header
    data = hdu[0].data.copy()
    naxis = hdr["NAXIS"]

    if axis is None and ctype is None:
        raise SystemExit('Please specify either axis or ctype')
    # find axis via CTYPE key

    if ctype:
        for i in range(1, naxis + 1):
            if hdr['CTYPE%d' % i].lower().startswith(ctype.lower()):
                axis = naxis - i  # fits to numpy indexing
    elif fits:
        axis = naxis - axis

    crval = hdr['CRVAL%d' % (naxis - axis)]
    cdelt = hdr['CDELT%d' % (naxis - axis)]
    crpix = hdr['CRPIX%d' % (naxis - axis)]
    # shift crval to crpix=1
    crval = crval - (crpix - 1) * cdelt

    nstacks = hdr['NAXIS%d' % (naxis - axis)]
    nchunks = nstacks // each_chunk
    print(("The FITS file %s has %s stacks along this axis. Breaking it up to %d images" % (fitsname, nstacks, nchunks)))

    outfiles = []
    for i in range(0, nchunks):
        _slice = [slice(None)] * naxis
        _slice[axis] = list(range(i * each_chunk, (i + 1) * each_chunk if i + 1 != nchunks else nstacks))
        hdu[0].data = data[_slice].astype(numpy.float32)
        hdu[0].header['CRVAL%d' % (naxis - axis)] = crval + i * cdelt * each_chunk
        hdu[0].header['CRPIX%d' % (naxis - axis)] = 1
        outfile = '%s-%04d.fits' % (prefix, i)
        outfiles.append(outfile)
        print(("Making chunk %d : %s. File is %s" % (i, repr(_slice[axis]), outfile)))
        hdu.writeto(outfile, overwrite=True)
    hdu.close()

    if keep_old == False:
        os.system('rm -rf %s' % fitsname)

    return outfiles


def reorder(fitsname, order=[], outfile=None):
    """ 
    Re-order FITS image axes. 

    Example:
        If your FITS files has axes (NAXIS{1,2,3,4})-> RA, DEC, STOKES, FREQ
        You can re-order it so that it has RA, DEC, FREQ, STOKES
        by running: reorder(fitsname, order=[1,2,4,3], "my_reodered_image.fits")
    """

    # Get to know input FITS image
    hdu = pyfits.open(fitsname)
    hdr0 = hdu[0].header
    ndim = hdr0["NAXIS"]
    order0 = list(range(1, ndim + 1))

    def fits2py(a):
        if isinstance(a, (list, tuple)):
            return list(ndim - numpy.array(a))
        else:
            return ndim - a

    if order == order0:
        print("The FITS image already has the order you requested")
        return
    elif len(order) != ndim:
        raise ValueError("The dimensions of the FITS image do not match your request.")

    # Ok, Lock and Load
    data = hdu[0].data.astype(numpy.float32)
    # First re-order data
    hdu[0].data = numpy.transpose(data, list((ndim - numpy.array(order))[list(range(ndim - 1, -1, -1))]))

    hdr = hdr0.copy()
    mendatory = "CTYPE CRVAL CDELT CRPIX".split()
    optional = "CUNIT CROTA".split()

    for key in mendatory + optional:
        for old, new in zip(order0, order):
            if old == new:
                continue
            try:
                val = hdr0["%s%d" % (key, old)]
                comment = hdr0.comments["%s%d" % (key, old)]
                idx = hdr.index("%s%d" % (key, old))
                # Seems like can't replace non empy val with an empty string. So delete first

                del hdr["%s%d" % (key, new)]
                hdr.insert(idx, ("%s%d" % (key, new), val, comment))

            except KeyError:
                if key in mendatory:
                    raise KeyError("ARBORTNG: FITS file doesn't have the '%s' key in the header" % key)
                else:
                    pass

    hdu[0].header = hdr
    outfile = outfile or "reodered_" + fitsname
    print(("Successfully re-ordered axes in FITS image. Output image is at %s" % outfile))
    hdu.writeto(outfile, overwrite=True)


def main():
    # setup some standard command-line option parsing
    #
    from optparse import OptionParser
    parser = OptionParser(usage="""%prog: [options] <image names...>""")
    parser.add_option("-o", "--output", dest="output", type="string",
                      help="name of output FITS file")
    parser.add_option("-r", "--replace", action="store_true",
                      help="replace (first) input file by output. Implies '--force'.")
    parser.add_option("-f", "--force", dest="force", action="store_true",
                      help="overwrite output file even if it exists")
    parser.add_option("-S", "--sanitize", type="float", metavar="VALUE",
                      help="sanitize FITS files by replacing NANs and INFs with VALUE")
    parser.add_option("-Z", "--zero-to-nan", action='store_true',
                      help="Replace zeros with NaNs")
    parser.add_option("-N", "--nonneg", action="store_true",
                      help="replace negative values by 0")
    parser.add_option("-m", "--mean", dest="mean", action="store_true",
                      help="take mean of input images")
    parser.add_option("-d", "--diff", dest="diff", action="store_true",
                      help="take difference of 2 input images")
    parser.add_option("--div", dest="div", action="store_true",
                      help="take ratio of 2 input images (e.g. as in primary beam correction)")
    parser.add_option("--div-min", type="float",
                      help="min value to use for second image, when taking a ratio")
    parser.add_option("--div-min-mask", type="float", metavar="MIN", default=0,
                      help="mask value to apply to output where second image < MIN")
    parser.add_option("--sum", dest="sum", action="store_true",
                      help="sum of input images")
    parser.add_option("--prod", dest="prod", action="store_true",
                      help="product of input images")
    parser.add_option("-t", "--transfer", action="store_true",
                      help="transfer data from image 2 into image 1, preserving the FITS header of image 1")
    parser.add_option("-z", "--zoom", dest="zoom", type="int", metavar="NPIX",
                      help="zoom into central region of NPIX x NPIX size")
    parser.add_option("-p", "--paste", action="store_true",
                      help="paste image(s) into first image, using the nearest WCS pixel positions. See also --empty-canvas")
    parser.add_option("--empty-canvas", action="store_true",
                      help="null first image, then paste image(s) into it")
    parser.add_option("--null-region", type="str", metavar="NXxNY", help="null central region of specified size")
    parser.add_option("-M", "--set-max", metavar="MAX", type="float",
                      help="set values over MAX to MAX")
    parser.add_option("-R", "--rescale", dest="rescale", type="float",
                      help="rescale image values")
    parser.add_option("-E", "--edit-header", metavar="KEY=VALUE", type="string", action="append",
                      help="replace header KEY with VALUE. Use KEY=VALUE for floats and KEY='VALUE' for strings.")
    parser.add_option("-D", "--delete-header", metavar="KEY", type="string", action="append",
                      help="Remove KEY from header")
    parser.add_option("--stack", metavar="outfile:axis",
                      help="Stack a list of FITS images along a given axis. This axis may given as an integer"
                           "(as it appears in the NAXIS keyword), or as a string (as it appears in the CTYPE keyword)")
    parser.add_option("--reorder",
                      help="Required order. List of comma seperated indeces")
    parser.add_option("--add-axis", metavar="CTYPE:CRVAL:CRPIX:CDELT[:CUNIT:CROTA]", action="append", default=[],
                      help="Add axis to a FITS image. The AXIS will be described by CTYPE:CRVAL:CRPIX:CDELT[:CUNIT:CROTA]. "
                           "The keywords in brackets are optinal, while those not in brackets are mendatory. "
                           "This axis will be the last dimension. Maybe specified more than once.")
    parser.add_option("--unstack", metavar="prefix:axis:each_chunk",
                      help="Unstack a FITS image into smaller chunks each having [each_chunk] planes along a given axis. "
                           "This axis may given as an integer (as it appears in the NAXIS keyword), or as a string "
                           "(as it appears in the CTYPE keyword)")
    parser.add_option("--delete-files", action="store_true",
                      help="Delete original file(s) after stacking/unstacking using --stack/--unstack")
    parser.add_option("-H", "--header", action="store_true", help="print header(s) of input image(s)")
    parser.add_option("-s", "--stats", action="store_true",
                      help="print stats on images and exit. No output images will be written.")
    parser.add_option("-F", "--file_pattern",
                      help="Speicfy input images via a pattern string, e.g, \"prefix*June2016.fits\". NB: The qouatation marks are important.")

    parser.set_defaults(output="", mean=False, zoom=0, rescale=1, edit_header=[], delete_header=[])

    (options, imagenames) = parser.parse_args()

    if options.file_pattern:
        imagenames = glob.glob(options.file_pattern)

    if not imagenames:
        parser.error("No images specified. Use '-h' for help.")

    # print "%d input image(s): %s"%(len(imagenames),", ".join(imagenames))
    images = [pyfits.open(img) for img in imagenames]
    updated = False

    # Stack FITS images
    if options.stack:
        if len(imagenames) < 1:
            parser.error("Need more than one image to stack")
        stack_args = options.stack.split(":")
        if len(stack_args) != 2:
            parser.error("Two --stack options are required. See ./fitstool.py -h")

        outfile, axis = stack_args

        _string = True
        try:
            axis = int(axis)
            _string = False
        except ValueError:
            _string = True

        stack_planes(imagenames, ctype=axis if _string else None, keep_old=not options.delete_files,
                     axis=None if _string else axis, outname=outfile, fits=True)

    # Unstack FITS image
    if options.unstack:
        image = imagenames[0]
        if len(imagenames) < 1:
            parser.error("Need more than one image to stack")
        unstack_args = options.unstack.split(":")
        if len(unstack_args) != 3:
            parser.error("Two --unstack options are required. See ./fitstool.py -h")

        prefix, axis, each_chunk = unstack_args

        _string = True
        try:
            axis = int(axis)
            _string = False
        except ValueError:
            _string = True

        each_chunk = int(each_chunk)

        unstack_planes(image, each_chunk, ctype=axis if _string else None,
                       axis=None if _string else axis, prefix=prefix, fits=True,
                       keep_old=not options.delete_files)
        sys.exit(0)

    if options.add_axis:
        for axis in options.add_axis:
            hdu = pyfits.open(imagenames[0])
            hdr = hdu[0].header
            ndim = hdr["NAXIS"]
            hdu[0].data = hdu[0].data[numpy.newaxis, ...]
            _mendatory = "CTYPE CRVAL CDELT CRPIX".split()
            _optional = "CUNIT CROTA".split()
            L = len(_mendatory + _optional)
            values = axis.split(":")
            Lv = len(values)

            if len(_mendatory) > len(values):
                parser.error("Something with the way specified --add-axis. See %prog -h for help")

            for i, value in enumerate(values):
                try:
                    value = float(value)
                except ValueError:
                    if isinstance(value, str):
                        value = value.upper()

                hdu[0].header["%s%d" % ((_mendatory + _optional)[i], ndim + 1)] = value

            hdu.writeto(imagenames[0], overwrite=True)
            print(("Successfully added axis %s to %s" % (values[0], imagenames[0])))

    if options.reorder:
        for image in imagenames:
            order = list(map(int, options.reorder.split(",")))
            reorder(image, order=order, outfile=image)

    if options.header:
        for filename, img in zip(imagenames, images):
            img.verify('silentfix')
            if len(imagenames) > 1:
                print("======== FITS header for", filename)
            for hdrline in img[0].header.cards:
                print(hdrline)

    if options.replace:
        if options.output:
            parser.error("Cannot combine -r/--replace with -o/--output")
        outname = options.output or imagenames[0]
        options.force = True
        autoname = False
    else:
        outname = options.output
        autoname = not outname
        if autoname:
            outname = re.split('[_]', imagenames[0], 1)[-1]

    for keyval in options.edit_header:
        key, val = keyval.split("=")
        q = ''
        if val[0] == "'" and val[-1] == "'":
            images[0][0].header[key] = val[1:-1:]
            q = '"'
        elif val[-1] == 'd' or key.startswith('NAXIS'):
            images[0][0].header[key] = int(val[:-1] if val[-1] == 'd' else val)
        else:
            try:
                images[0][0].header[key] = float(val)
            except:
                images[0][0].header[key] = val
                q = '"'
        print("Setting header %s=%s%s%s" % (key, q, val, q))
        updated = True

    for key in options.delete_header:
        try:
            del images[0][0].header[key]
        except KeyError:
            raise "Key '%s' not found in FITS header" % key
        print("Deleting key '%s' header" % key)
        updated = True

    if options.sanitize is not None:
        print("Sanitizing: replacing INF/NAN with", options.sanitize)
        for img in images:
            d = img[0].data
            d[numpy.isnan(d) + numpy.isinf(d)] = options.sanitize
        # if using stats, do not generate output
        if not options.stats:
            updated = True

    if options.zero_to_nan:
        print("Replacing zeros with NaNs")
        for img in images:
            d = img[0].data
            d[d == 0] = numpy.nan
        # if using stats, do not generate output
        if not options.stats:
            updated = True

    if options.nonneg:
        print("Replacing negative value by 0")
        for img, name in zip(images, imagenames)[:1]:
            d = img[0].data
            wh = d < 0
            d[wh] = 0
            print("Image %s: replaced %d points" % (name, wh.sum()))
        updated = True

    def compute_in_out_slice(N, N0, R, R0):
        """
        given an input axis of size N, and an output axis of size N0, and reference pixels of I and I0
        respectively, computes two slice objects such that
                A0[slice_out] = A1[slice_in]
        would do the correct assignment (with I mapping to I0, and the overlapping regions transferred)
        """
        i, j = 0, N     # input slice
        i0 = R0 - R
        j0 = i0 + N     # output slice
        if i0 < 0:
            i = -i0
            i0 = 0
        if j0 > N0:
            j = N - (j0 - N0)
            j0 = N0
        if i >= j or i0 >= j0:
            return None, None
        return slice(i0, j0), slice(i, j)

    if options.transfer:
        if len(images) != 2:
            parser.error("The --transfer option requires exactly two input images.")
        if autoname:
            outname = "xfer_" + outname
        print("Transferring %s into coordinate system of %s" % (imagenames[1], imagenames[0]))
        images[0][0].data = images[1][0].data
        updated = True
    elif options.diff:
        if len(images) != 2:
            parser.error("The --diff option requires exactly two input images.")
        if autoname:
            outname = "diff_" + outname
        print("Computing difference")
        data = images[0][0].data
        data -= images[1][0].data
        updated = True
    elif options.sum:
        if autoname:
            outname = "sum_" + outname
        print("Computing sum")
        data = images[0][0].data
        for d in images[1:]:
            data += d[0].data
        updated = True
    elif options.div:
        if len(images) != 2:
            parser.error("The --div option requires exactly two input images.")
        if autoname:
            outname = "div_" + outname
        print("Computing ratio")
        a, b = images[0][0].data, images[1][0].data
        if options.div_min is not None:
            b[b<=options.div_min] = options.div_min
        a /= b
        if options.div_min_mask is not None:
            a[b<=options.div_min] = options.div_min_mask
        updated = True
    elif options.prod:
        if autoname:
            outname = "prod_" + outname
        print("Computing product")
        data = images[0][0].data
        for d in images[1:]:
            data *= d[0].data
        updated = True
    elif options.mean:
        if autoname:
            outname = "mean%d_" % len(images) + outname
        print("Computing mean")
        data = images[0][0].data
        for img in images[1:]:
            data += img[0].data
        data /= len(images)
        images = [images[0]]
        updated = True
    elif options.paste:
        if autoname:
            outname = "paste%d_" % len(images) + outname
        print("Pasting images into canvas given by {}".format(imagenames[0]))
        print("WARNING: the reference pixel of the pasted image will be aligned with the WCS of the canvas image,")
        print("         but no WCS reprojection will be attempted. Results may be grossly incorrect if your WCSs")
        print("         are sufficiently different. Use at own risk.")
        # get original image size and WCS
        data0 = images[0][0].data
        ny0, nx0 = data0.shape[-2:]
        hdr0 = images[0][0].header
        wcs0 = WCS(hdr0, mode="pyfits")
        if options.empty_canvas:
            print("  canvas image will be cleared initially")
            data0[:] = 0

        for image, imagename in zip(images[1:], imagenames[1:]):
            hdr = image[0].header
            data = image[0].data
            ny, nx = data.shape[-2:]
            rx, ry = int(hdr["CRPIX1"]-1), int(hdr["CRPIX2"]-1)
            # find coordinate of this image's reference pixel in the canvas WCS
            rx0, ry0 = wcs0.wcs2pix(hdr["CRVAL1"], hdr["CRVAL2"])  # get WCS of center pixel
            print("image {} rpix {},{} is at {},{} in canvas image".format(imagename, rx, ry, rx0, ry0))
            rx0 = int(round(rx0))
            ry0 = int(round(ry0))
            xout, xin = compute_in_out_slice(nx, nx0, rx, rx0)
            yout, yin = compute_in_out_slice(ny, ny0, ry, ry0)
            if xout is not None and yout is not None:
                print("  canvas X {} will be assigned from {}".format(xout, xin))
                print("  canvas Y {} will be assigned from {}".format(yout, yin))
                data0[..., yout, xout] += data[..., yin, xin]
                updated = True
            else:
                print("  image {} footprint is outside of canvas".format(imagename))

    if options.null_region:
        if len(images) > 1:
            print("Too many input images specified for this operation, at most 1 expected")
            sys.exit(2)
        data = images[0][0].data
        nx = data.shape[-1]
        ny = data.shape[-2]
        xc, yc = nx//2, ny//2
        rx, ry = map(int,options.null_region.split("x", 1))
        rx, ry = min(rx, nx), min(ry, ny)
        x0, y0 = xc - rx//2, yc - ry//2
        data[..., y0:y0+ry, x0:x0+rx] = 0
        print("Setting {}:{} {}:{} to zero".format(x0, x0+rx, y0, y0+ry))
        updated = True

    if options.zoom:
        z = options.zoom
        if autoname:
            outname = "zoom%d_" % z + outname
        if len(images) > 1:
            print("Too many input images specified for this operation, at most 1 expected")
            sys.exit(2)
        data = images[0][0].data
        nx = data.shape[-1]
        ny = data.shape[-2]
        # make output array of shape z x z
        zdata_shape = list(data.shape)
        zdata_shape[-1] = zdata_shape[-2] = z
        zdata = numpy.zeros(zdata_shape, dtype=data.dtype)
        # make input/output slices
        rx, ry = nx//2, ny//2
        rx0 = ry0 = z//2
        xout, xin = compute_in_out_slice(nx, z, rx, rx0)
        yout, yin = compute_in_out_slice(ny, z, ry, ry0)

        zdata[..., yout, xout] = data[..., yin, xin]

        # update header
        hdr = images[0][0].header
        # hdr1 = hdr.copy()
        # hdr1["NAXIS"] = 2
        # # newer astropy don't like our radio axes 3 and 4 :((
        # for nax in "3", "4":
        #     for key in "NAXIS", "CTYPE", "CRVAL", "CRPIX", "CDELT", "CROTA":
        #         hdr1.remove(key + nax, ignore_missing=True)
        wcs = WCS(hdr) #, mode="pyfits")
        pixcoord = numpy.zeros((1, hdr["NAXIS"]), float)
        pixcoord[0, 0] = nx//2
        pixcoord[0, 1] = ny//2
        world = wcs.wcs_pix2world(pixcoord, 0)  # get WCS of center pixel
        hdr["CRVAL1"] = world[0,0]
        hdr["CRVAL2"] = world[0,1]
        hdr["CRPIX1"] = rx0 + 1
        hdr["CRPIX2"] = ry0 + 1

        print("Making zoomed image of shape", "x".join(map(str, zdata.shape)))

        images = [pyfits.PrimaryHDU(zdata, hdr)]
        updated = True

    if options.set_max is not None:
        print("Enforcing maximum value of {}".format(options.set_max))
        images[0][0].data[images[0][0].data > options.set_max] = options.set_max
        updated = True

    if options.rescale != 1:
        if autoname and not updated:
            outname = "rescale_" + outname
        if len(images) > 1:
            "Too many input images specified for this operation, at most 1 expected"
            sys.exit(2)
        print("Applying scaling factor of %f to image values" % options.rescale)
        images[0][0].data *= options.rescale
        updated = True

    if updated:
        imagenames[0] = outname

    if options.stats:
        for ff, filename in zip(images, imagenames):
            data = ff[0].data
            min1, max1, dum1, dum2 = scipy.ndimage.measurements.extrema(data)
            sum = data.sum()
            mean = sum / data.size
            std = math.sqrt(((data - mean) ** 2).mean())
            print("%s: min %g, max %g, sum %g, np %d, mean %g, std %g" % (filename, min1, max1, sum, data.size, mean, std))
        sys.exit(0)

    if updated:
        print("Writing output image", outname)
        if os.path.exists(outname) and not options.force:
            print("Output image exists, rerun with the -f switch to overwrite.")
            sys.exit(1)
        images[0].writeto(outname, overwrite=True)
    elif not (options.header or options.stack or options.add_axis or options.reorder):
        print("No operations specified. Use --help for help.")
