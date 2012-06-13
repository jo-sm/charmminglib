#!/usr/bin/env python
"""
Installation script for standalone usage of pychmlib.

:Author: fcp
:Date: 10/28/2010
"""


from distutils.core import setup
import sys


# Make sure python version is 2.6 or 2.7
ver = sys.version_info
if not (ver[0] == 2 and ver[1] in [6,7]):
    print '\npychm requires python version 2.6.x or 2.7.x to function properly\n'
    sys.exit(1)


setup(
    name = 'pychm',
    version = '0.2',
    description = 'Python backend to pychm',
    author = 'Frank C. Pickard IV and Tim Miller',
    author_email = 'pickard81@gmail.com, btamiller@gmail.com',
    url = 'www.pychm.googlecode.com/svn/',
    packages = [
        'pychm',
        'pychm.analysis',
        'pychm.cg', 'pychm.cg.analysis',
        'pychm.const',
        'pychm.future', 'pychm.future.io', 'pychm.future.io.charmm', 'pychm.future.lib', 'pychm.future.scripts',
        'pychm.io',
        'pychm.lib',
        'pychm.scripts'
        ],
    package_dir = {'pychm': 'src/pychm'},
    package_data = {'pychm': ['data/*']},
    requires = ['numpy'],
    license = 'Public Domain',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Framework :: Django',
        'Intended Audience :: Science/Research',
        'License :: Public Domain',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python'
        'Topic :: Scientific/Engineering :: Molecular Modeling'
        ]
    )
