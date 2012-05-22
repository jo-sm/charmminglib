"""
Insertion point for altering charmm specific I/O classes.
"""

from __future__ import division

__all__ = []

from collections import deque

from pychm.future.io.base import FortFile, TextFile


class CharmmBin(FortFile):
    """
    """
    def __init__(self, fname, mode='rb', buffering=None,
                endian='@', rec_head_prec='i'):
        super(CharmmBin, self).__init__(fname=fname, mode=mode,
            buffering=buffering, endian=endian, rec_head_prec=rec_head_prec)


class CharmmCard(TextFile):
    """
    """
    comment_char = '!'
    continue_char = '-'

    def __init__(self, fname, mode='r', buffering=None):
        super(CharmmCard, self).__init__(fname=fname, mode=mode,
                                        buffering=buffering)
        self.title = None
        self.body = None
        self.deque = None

    def iter_normalize_card(self):
        """Returns an iterator over the card file's lines, with comments,
        exterior whitespace and blanklines removed. Also, lines that are
        broken over multiple physical lines in the file are combined.
        """
        cc = self.comment_char
        self.seek(0, 0)
        # filter comments, whitespace, blanklines
        iterable = ( line.strip() for line in self )
        iterable = ( line for line in iterable if line )
        iterable = ( line.lower() for line in iterable if not line.startswith(cc) )
        iterable = ( line.split(cc)[0] for line in iterable )
        iterable = ( line.strip() for line in iterable )
        # account for line continuation chars
        tmp = []
        for line in iterable:
            if line.endswith(self.continue_char):
                tmp.append(line[:-1].strip())
            else:
                if tmp:
                    tmp.append(line)
                    yield ' '.join(tmp)
                    tmp = []
                else:
                    yield line
        # flush
        if tmp:
            yield ' '.join(tmp)

    def parse(self):
        self.deque = deque(self.iter_normalize_card())
        title = []
        while 1:
            if self.deque[0].startswith('*'):
                title.append(self.deque.popleft()[1:].strip())
            else:
                break
        self.title = [ line.strip() for line in title if line ]

    def pack_title(self):
        try:
            tmp = []
            for line in self.title:
                tmp.append('* %s' % line)
            tmp.append('*\n')
            return '\n'.join(tmp)
        except TypeError:
            return "* A blank title.\n*\n"
