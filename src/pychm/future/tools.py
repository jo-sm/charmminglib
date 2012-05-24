"""
This is the module to put commonly used utility functions, classes and
generators. These objects should have very broad applicability.

This module should *NOT* ever make any imports from other pychm modules.
"""


def rwprop(func):
    """A decorator similar to the built in `property`, however this one allows
    setting fset, fdel and doc in addition to fget without namespace pollution.
    """
    return property(**func())


class _mydict(dict):
    def __getitem__(self, key):
        try:
            return super(_mydict, self).__getitem__(key)
        except KeyError:
            return None


def paragraphs(iterable, splitter):
    """
    Cut a text stream up into 'paragraphs,' where partitions are
    determined by a :mod:`list` named `splitter`.

    >>> iterable = paragraphs(iterable, ['taco', 'beans'])
    """
    assert isinstance(splitter, (tuple, list))
    splitter = tuple(splitter)
    paragraph = []
    for line in iterable:
        if line.startswith(splitter):
            if paragraph:
                yield paragraph
            paragraph = [line]
        else:
            paragraph.append(line)
    if paragraph:
        yield paragraph
