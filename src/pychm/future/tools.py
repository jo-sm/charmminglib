"""
This is the module to put commonly used utility functions, classes and
generators. These objects should have very broad applicability.

This module should *NEVER* import from other pychm modules.
"""

from os.path import expanduser, abspath

def rwprop(func):
    """A decorator similar to the built in `property`, however this one allows
    setting fset, fdel and doc in addition to fget without namespace pollution.
    """
    return property(**func())


class _mydict(dict):
    """The same as a regular class:`dict`, except instead of raising exec:`KeyError`
    `None` is returned.
    """
    def __getitem__(self, key):
        try:
            return super(_mydict, self).__getitem__(key)
        except KeyError:
            return None


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
