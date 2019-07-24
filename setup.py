#!/usr/bin/env python

import os
from setuptools import setup

install_requires = [
      'astropy',
      'numpy',
      'matplotlib',
      'python-casacore',
      'meqtrees_cattery',
      'scipy',
      'astlib',
      'astro-kittens',
      'six',
      'future',
]

setup(name='owlcat',
      version='1.5.2',
      description='miscellaneous utility scripts for manipulating radio interferometry data',
      author='Oleg Smirnov',
      author_email='Oleg Smirnov <osmirnov@gmail.com>',
      url='https://github.com/ska-sa/owlcat',
      packages=['Owlcat'],
      install_requires=install_requires,
      scripts=['Owlcat/bin/' + i for i in os.listdir('Owlcat/bin')],
      data_files=[('Owlcat/bin/', ['Owlcat/bin/commands.list'])],
)
