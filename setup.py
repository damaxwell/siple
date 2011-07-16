#!/usr/bin/env python

"""Distutils based setup script for siple.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command (you'll probably need root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands, e.g.:

    python setup.py clean -> will clean all trash (*.pyc and stuff)

To get a full list of avaiable commands, read the output of:

from distutils.core import setup
"""


from distutils.core import setup


setup(
    name='siple',
    description='siple: a small inverse problems library',
    version='0.1',
    packages=['siple', 'siple.gradient', 'siple.linalg', 'siple.opt'],
    license='GNU Public License v2',
    author="David Maxwell",
    author_email="damaxwell@alaska.edu",
    long_description=open('README').read(),
)
