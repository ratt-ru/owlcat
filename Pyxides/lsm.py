from Pyxis.ModSupport import *

import imager,std

import pyfits
import Tigger

register_pyxis_module(superglobals="OUTFILE");

tigger_restore  = x("tigger-restore")
tigger_convert  = x("tigger-convert")
tigger_tag      = x("tigger-tag")

v.define("LSM","lsm.lsm.html",
  """current local sky model""");

define("LSM_TDL_Template","tiggerlsm.filename=$LSM",
  """TDL option for selecting current lsm""");
define("LSMREF","",
  """reference LSM (for transferring tags, etc.)""");
  
define('PYBDSM_OUTPUT_Template',"${OUTFILE}_pybdsm.lsm.html",
  """output LSM file for pybdsm search""");    
define('PYBDSM_POLARIZED',0,
  """set to True to run pybdsm in polarized mode""");
_pybdsm = x.pybdsm;

define('CLUSTER_DIST',0,
  """source clustering distance. If 0, then CLUSTER_DIST_BEAMS is used instead.""");
define('CLUSTER_DIST_BEAMS',3,
  """source clustering distance, in terms of number of PSFs (measured as (BMAJ+BMIN)/2). If BMAJ/BMIN is not defined,
  falls back to 60".""");
define('MIN_EXTENT',0,
  """minimum Gaussian source extent; sources smaller than this will be converted to point sources""");

def pybdsm_search (image="${imager.RESTORED_IMAGE}",output="$PYBDSM_OUTPUT",pol='$PYBDSM_POLARIZED',threshold=None,pbexp=None,**kw):
  """Runs pybdsm on the specified 'image', converts the results into a Tigger model and writes it to 'output'.
  Use 'threshold' to specify a non-default threshold (thresh_isl and thresh_pix).
  Use 'pol' to force non-default polarized mode.
  Use 'pbexp' to supply a primary beam expression (passed to tigger-convert), in which case the output model will contain
  intrinsic fluxes.
  """
  image,output,pol = interpolate_locals("image output pol");
  # setup parameters
  script = II("${output:FILE}.pybdsm");
  gaul = II("${output:FILE}.gaul");
  if threshold:
    kw['thresh_isl'] = kw['thresh_pix'] = threshold;
  kw['polarisation_do'] = is_true(pol);
  # run pybdsm
  from lofar import bdsm
  img = bdsm.process_image(image,**kw);
  img.write_catalog(outfile=gaul,format='ascii',catalog_type='gaul',clobber=True);
  # set clustering parameter from beam size
  cluster = CLUSTER_DIST;
  if not cluster:
    hdr = pyfits.open(image)[0].header;
    # BMAJ/BMIN is in degrees -- convert to seconds, or fall back to 60" if not set
    cluster = 1800*(hdr.get('BMAJ',0)+hdr.get('BMIN',0))*CLUSTER_DIST_BEAMS or 60;
  # convert catalog
  if pbexp:
    args = [ "--primary-beam",pbexp,"--app-to-int" ]
  else:
    args = []
  tigger_convert(gaul,output,"-t","ASCII","--format",
      "name Isl_id Source_id Wave_id ra_d E_RA dec_d E_DEC i E_Total_flux Peak_flux E_Peak_flux Xposn E_Xposn Yposn E_Yposn Maj E_Maj Min E_Min PA E_PA emaj_d E_DC_Maj emin_d E_DC_Min pa_d E_DC_PA Isl_Total_flux E_Isl_Total_flux Isl_rms Isl_mean Resid_Isl_rms Resid_Isl_mean S_Code"
     +
    ("q E_Total_Q u E_Total_U v E_Total_V Linear_Pol_frac Elow_Linear_Pol_frac Ehigh_Linear_Pol_frac "+
     "Circ_Pol_Frac Elow_Circ_Pol_Frac Ehigh_Circ_Pol_Frac Total_Pol_Frac Elow_Total_Pol_Frac Ehigh_Total_Pol_Frac Linear_Pol_Ang E_Linear_Pol_Ang"
    if pol else ""),
    "-f","--rename",
    "--cluster-dist",cluster,
    "--min-extent",MIN_EXTENT,
    split_args=False,
    *args);
    
document_globals(pybdsm_search,"PYBDSM_* imager.RESTORED_IMAGE");



def transfer_tags (fromlsm="$LSMREF",lsm="$LSM",output="$LSM",tags="dE",tolerance=60*ARCSEC):
  """Transfers tags from a reference LSM to the given LSM. That is, for every tag
  in the given list, finds all sources with those tags in 'fromlsm', then applies 
  these tags to all nearby sources in 'lsm' (within a radius of 'tolerance'). 
  Saves the result to an LSM file given by 'output'.
  """
  fromlsm,lsm,output,tags = interpolate_locals("fromlsm lsm output tags");
  # now, set dE tags on sources
  import Tigger
  refmodel = Tigger.load(fromlsm);
  model = Tigger.load(lsm);
  tagset = frozenset(tags.split());
  info("Transferring tags %s from %s to %s"%(",".join(tagset),fromlsm,lsm));
  # for each dE-tagged source in the reference model, find all nearby sources
  # in our LSM, and tag them
  for src0 in refmodel.getSourceSubset(",".join(["="+x for x in tagset])):
    for src in model.getSourcesNear(src0.pos.ra,src0.pos.dec,tolerance=tolerance):
      for tag in tagset:
        tagval = src0.getTag(tag,None);
        if tagval is not None:
          if src.getTag(tag,None) != tagval:
            src.setTag(tag,tagval);
            info("setting tag %s=%s on source %s (from reference source %s)"%(tag,tagval,src.name,src0.name))
  model.save(output);

  

CC_RESCALE = 1.  
CC_IMAGE_Template = "${OUTFILE}_ccmodel.fits"

def add_ccs (lsm="$LSM",filename="${imager.MODEL_IMAGE}",
             cc_image="$CC_IMAGE",srcname="ccmodel",output="$LSM",zeroneg=True,scale=None,pad=1):
  """Adds clean components from the specified FITS image 'filename' to the sky model given by 'lsm'.
  Saves the result to an LSM file given by 'output'.
  The CC image is copied to 'cc_image', optionally rescaled by 'scale', and optionally has negative pixels reset to zero (if zeroneg=True).
  'srcname' gives the name of the resulting LSM component.
  'pad' gives the padding attribute of the LSM component, use e.g. 2 if CC image has significant signal towards the edges.
  """;
  lsm,filename,cc_image,srcname,output = interpolate_locals("lsm filename cc_image srcname output");
  # rescale image
  ff = pyfits.open(filename);
  ff[0].data *= (scale if scale is not None else CC_RESCALE);
  if zeroneg:
    ff[0].data[ff[0].data<0] = 0;
  ff.writeto(cc_image,clobber=True);
  info("adding clean components from $filename ($cc_image), resulting in model $output");
  tigger_convert(lsm,output,"-f","--add-brick","$srcname:$cc_image:%f"%pad);

document_globals(add_ccs,"MODEL_CC_*");  


def pointify (lsm="$LSM",output="$LSM",name=""):
  """Replaces names sources with point sources""";
  lsm,output,name = interpolate_locals("lsm output name");
  model = Tigger.load(lsm);
  src = model.findSource(name);
  info("Setting source $name in model $lsm to point source, saving to $output");
  src.shape = None;
  model.save(output);

document_globals(add_ccs,"MODEL_CC_*");  

  
  