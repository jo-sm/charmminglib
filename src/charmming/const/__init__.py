"""
Miscellaneous constants are contained herein.  Biological, chemical and physical
constants are contained in the appropriate submodule.  All :class:`string`
are lowercase.

**Attributes:**
    | ``alphanum`` - A :class:`string` from `a` to `9`.
    | ``alpha`` - A :class:`string` from `a` to `z`.
    | ``alpha2num`` - A :class:`dict` mapping ``alpha`` to the 0-based index.
    | ``alphanum2num`` - A :class:`dict` mapping ``alphanum`` to the 0-based index.
"""


import charmming.const.bio
import charmming.const.units


__all__ = ['bio', 'units']


alphanum = 'abcdefghijklmnopqrstuvwxyz0123456789'


alpha = 'abcdefghijklmnopqrstuvwxyz'


alpha2num = {
    '?': -1,
    'a': 0,
    'b': 1,
    'c': 2,
    'd': 3,
    'e': 4,
    'f': 5,
    'g': 6,
    'h': 7,
    'i': 8,
    'j': 9,
    'k': 10,
    'l': 11,
    'm': 12,
    'n': 13,
    'o': 14,
    'p': 15,
    'q': 16,
    'r': 17,
    's': 18,
    't': 19,
    'u': 20,
    'v': 21,
    'w': 22,
    'x': 23,
    'y': 24,
    'z': 25
    }


alphanum2num = {
     '?': -1,
     'a': 0,
     'b': 1,
     'c': 2,
     'd': 3,
     'e': 4,
     'f': 5,
     'g': 6,
     'h': 7,
     'i': 8,
     'j': 9,
     'k': 10,
     'l': 11,
     'm': 12,
     'n': 13,
     'o': 14,
     'p': 15,
     'q': 16,
     'r': 17,
     's': 18,
     't': 19,
     'u': 20,
     'v': 21,
     'w': 22,
     'x': 23,
     'y': 24,
     'z': 25,
     '0': 26,
     '1': 27,
     '2': 28,
     '3': 29,
     '4': 30,
     '5': 31,
     '6': 32,
     '7': 33,
     '8': 34,
     '9': 35
    }
