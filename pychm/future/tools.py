"""
This is the module to put commonly used utility functions, classes and
generators. These objects should have very broad applicability.

This module should *NEVER* import from other pychm modules.
"""


from contextlib import contextmanager
from os.path import expanduser, abspath
import logging


logger = logging.getLogger('pychm.tools')


def rwprop(func):
    """A decorator similar to the built in `property`, however this one allows
    setting fset, fdel and doc in addition to fget without namespace pollution.
    """
    return property(**func())


class mydict(dict):
    """The same as a regular class:`dict`, except instead of raising exec:`KeyError`
    `None` is returned.
    """
    def __getitem__(self, key):
        try:
            return super(mydict, self).__getitem__(key)
        except KeyError:
            return None


def myfloat(k):
    """Casts a variable to a float, but doesn't barf exceptions on
    uninitialized (`None`) values.
    """
    try:
        return float(k)
    except TypeError:
        if k is None:
            return None
        else:
            raise


def myint(k):
    """Casts a variable to an int, but doesn't barf exceptions on
    uninitialized (`None`) values.
    """
    try:
        return int(k)
    except TypeError:
        if k is None:
            return None
        else:
            raise


def _myexpandpath(path):
    if '~' in path:
        return expanduser(path)
    return abspath(path)


def paragraphs(iterable, splitter):
    """
    Cut a text stream up into 'paragraphs,' where partitions are
    determined by a :mod:`list` named `splitter`.

    >>> iterable = paragraphs(iterable, ['taco', 'beans'])
    """
    if isinstance(splitter, basestring):
        splitter = (splitter,)
    else:
        splitter = tuple(splitter)
    tmp = []
    for line in iterable:
        if line.startswith(splitter):
            if tmp:
                yield tmp
            tmp = [line]
        else:
            tmp.append(line)
    if tmp:
        yield tmp


def _myopenzip(fname, mode='r', buffering=8192, ftype=None):
    import zipfile
    import tarfile
    import gzip
    import bz2


    if not isinstance(fname, basestring):
        raise TypeError("Invalid fname")
    if not isinstance(mode, basestring):
        raise TypeError("Invalid mode")
    if not isinstance(buffering, int):
        raise TypeError("Invalid buffering")
    if ftype not in [None, 'zip', 'tar', 'gz', 'gzip', 'bz', 'bz2',
                    'bzip', 'bzip2', 'auto']:
        raise ValueError("Invalid ftype: %r" % ftype)
    #
    if ftype is None:
        logger.info("opening file: %s, without compression" % fname)
        fp = open(fname, mode=mode, buffering=buffering)
    elif ftype == 'zip':
        logger.info("opening file: %s, with zip compression" % fname)
        fp = zipfile.ZipFile(fname, mode=mode)
    elif ftype == 'tar':
        logger.info("opening file: %s, with tar compression" % fname)
        fp = tarfile.open(fname, mode=mode, bufsize=buffering)
    elif ftype in ['gz', 'gzip']:
        logger.info("opening file: %s, with gzip compression" % fname)
        fp = gzip.open(fname, mode=mode)
    elif ftype in ['bz', 'bz2', 'bzip', 'bzip2']:
        logger.info("opening file: %s, with bzip compression" % fname)
        fp = bz2.BZ2File(fname, mode=mode, buffering=buffering)
    elif ftype == 'auto':
        if zipfile.is_zipfile(fname):
            logger.info("opening file: %s, with zip compression" % fname)
            fp = zipfile.ZipFile(fname, mode=mode)
        elif tarfile.is_tarfile(fname):
            logger.info("opening file: %s, with tar compression" % fname)
            fp = tarfile.open(fname, mode=mode, bufsize=buffering)
        elif fname.endswith(('.gz', '.gzip')):
            logger.info("opening file: %s, with gzip compression" % fname)
            fp = gzip.open(fname, mode=mode)
        elif fname.endswith(('.bz', '.bz2', '.bzip', '.bzip2')):
            logger.info("opening file: %s, with bzip compression" % fname)
            fp = bz2.BZ2File(fname, mode=mode, buffering=buffering)
        else:
            logger.info("opening file: %s, without compression" % fname)
            fp = open(fname, mode=mode, buffering=buffering)
    else:
        raise ValueError("??? Invalid ftype: %r" % ftype)
    return fp
