Dummy PR
Installing Owlcat
=================


Prerequisites
-------------

* Assorted python packages: pyfits, numpy, matplotlib, casacore, pyrap.

 There are debian packages available for casacore and pyrap at:

 https://launchpad.net/~ska-sa/+archive/main

* MeqTrees Cattery, which you can get here:

  https://github.com/ska-sa/meqtrees-cattery

  For this software there are also Debian packages available in the repo above.


Getting Owlcat
--------------

  Download the tarball or do a git checkout from:

  https://github.com/ska-sa/owlcat


Installing Olwcat
-----------------

to intall the package run:

  $ python setup.py install

Or add the source to your PYTHONPATH.


Running Owlcat
--------------

Currently, the most useful script in Owlcat is plot-ms.py. See 

http://www.astron.nl/meqwiki-data/users/oms/Owlcat-plotms-tutorial.purrlog/

for examples. If you need quick deterministic flagging, try flag-ms.py.


Questions or problems?
----------------------

Just e-mail Oleg Smirnov <osmirnov@gmail.com> or open an issue on Github.


Development
-----------

[![Build Status](https://travis-ci.org/ska-sa/owlcat.svg?branch=master)](https://travis-ci.org/ska-sa/owlcat)
