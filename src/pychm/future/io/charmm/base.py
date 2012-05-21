"""
Insertion point for altering charmm specific I/O classes.
"""

from __future__ import division

__all__ = []

from io import DEFAULT_BUFFER_SIZE

from pychm.future.io.base import FortFile, TextFile


class CharmmBin(FortFile):
    """
    """
    def __init__(self, fname, mode='rb', buffering=DEFAULT_BUFFER_SIZE,
                endian='@', rec_head_prec='i'):
        super(CharmmBin, self).__init__(fname=fname, mode=mode,
            buffering=buffering, endian=endian, rec_head_prec=rec_head_prec)


class CharmmCard(TextFile):
    pass
