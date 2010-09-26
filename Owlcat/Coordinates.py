# -*- coding: utf-8 -*-

from math import *

def radec_to_lmn (ra,dec,ra0,dec0):
  """Returns l,m,n corresponding to direction ra,dec w.r.t. direction ra0,dec0""";
## our old formula, perhaps unjustly suspected by me
## See purrlog for 3C147_spw0, entries of Nov 21.
## Doesn't this break down at the pole (l always 0)?
  l = cos(dec) * sin(ra-ra0);
  m = sin(dec) * cos(dec0) - cos(dec) * sin(dec0) * cos(ra-ra0);
## Sarod's formula from LSM.common_utils. Doesn't seem to work right!
## (that's because it's for NCP lm coordinates used in NEWSTAR sky models)
#  l = sin(ra-ra0)*math.cos(dec);
#   sind0 = sin(dec0);
#   if sind0 != 0:
#     m = -(cos(ra-ra0)*cos(dec)-cos(dec0))/math.sin(dec0);
#   else:
#     m = 0
  n = sqrt(1-l*l-m*m);
  return l,m,n;

def lm_to_radec (l,m,ra0,dec0):
  """Returns ra,dec corresponding to l,m w.r.t. direction ra0,dec0""";
  # see formula at http://en.wikipedia.org/wiki/Orthographic_projection_(cartography)
  rho = sqrt(l**2+m**2);
  cc = asin(rho);

  ra = ra0 + atan2( l*sin(cc),rho*cos(dec0)*cos(cc)-m*sin(dec0)*sin(cc) );
  dec = asin( cos(cc)*sin(dec0) + m*sin(cc)*cos(dec0)/rho );

  return ra,dec;
