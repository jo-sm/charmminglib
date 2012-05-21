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
