"""
Insertion point for altering charmm specific I/O classes.
"""

from __future__ import division

__all__ = []

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

    def get_title(self):
        cur_pos = self.tell()
        self.seek(0, 0)
        title = []
        for line in self:
            line = line.strip()
            if not line:
                continue
            if line == '*':
                self.seek(cur_pos, 0)
                return title
            elif line.startswith('*'):
                title.append(line)
            else:
                self.seek(cur_pos, 0)
                return None

    def pack_title(self):
        try:
            tmp = []
            for line in self.title:
                if not line.startswith('*'):
                    line = '* ' + line
                tmp.append(line)
            if tmp[-1].strip() != '*':
                tmp.append('*')
            return '\n'.join(tmp)
        except TypeError:
            return "* A blank title.\n*\n"

    def seek_top(self):
        if self.title:
            self.readlines(len(self.title)+1)
        else:
            self.seek(0, 0)

    def iter_normalize_card_data(self):
        """Returns an iterator over the card file's lines, with comments,
        exterior whitespace and blanklines removed. Also, lines that are
        broken over multiple physical lines in the file are combined.
        """
        cc = self.comment_char
        self.seek_top()
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
                tmp.append(line)
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
