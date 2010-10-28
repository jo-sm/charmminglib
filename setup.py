#!/usr/bin/env python
"""
DOCME
"""
# fcp
# 10/28/2010


from distutils.core import setup
import sys


# Make sure python version is 2.6 or 2.7
ver = sys.version_info
if not (ver[0] == 2 and ver[1] in [6,7]):
    print '\ncharmming requires python version 2.6.x or 2.7.x to function properly\n'
    sys.exit(1)


setup(
    name = 'charmming',
    version = '0.1',
    description = 'Python backend to CHARMMing',
    author = 'Frank C. Pickard IV',
    author_email = 'pickard81@gmail.com',
    url = 'www.charmming.googlecode.com/svn/',
    packages = [
        'charmming', 'charmming.analysis', 'charmming.analysis.cg',
        'charmming.cg', 'charmming.const', 'charmming.io', 'charmming.lib',
        'charmming.scripts'
        ],
    package_dir = {'charmming': 'src/charmming'},
    package_data = {'charmming': ['data/*']},
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
