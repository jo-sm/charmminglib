"""
Insertion point for altering charmm specific I/O classes.
"""

from __future__ import division

__all__ = []

from collections import deque
import logging

from pychm.future.io.base import FortFile, TextFile
from pychm.future.io.charmm import COMMENT_CHAR as CC
from pychm.future.io.charmm import CONTINUE_CHAR as KC


logger = logging.getLogger('pychm.io.charmm.base')


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
    def __init__(self, fname, mode='r', buffering=None):
        super(CharmmCard, self).__init__(fname=fname, mode=mode,
                                        buffering=buffering)
        self.title = None
        self.version = None
        self.body = None
        self.deque = None

    def iter_normalize_card(self, comments=False):
        """Returns an iterator over the card file's lines, with comments,
        exterior whitespace and blanklines removed. Also, lines that are
        broken over multiple physical lines in the file are combined.
        """
        if comments:
            return self._iter_normalize_card_yes_comments()
        else:
            return self._iter_normalize_card_no_comments()

    def parse(self):
        self.deque = deque(self.iter_normalize_card())
        self.title = self._parse_title()
        self.version = self._parse_version()

    def _parse_title(self):
        logger.debug('reading title')
        tmp = []
        while self.deque[0].startswith('*'):
            tmp.append(self.deque.popleft()[1:].strip())
        return [ line for line in tmp if line ]

    def _parse_version(self):
        logger.debug('reading version info')
        test_line = self.deque[0]
        try:
            ver, subver = map(int, test_line.split())
            self.deque.popleft()
            return (ver, subver)
        except ValueError:
            try:
                ver = int(test_line)
                self.deque.popleft()
                return (ver, 0)
            except ValueError:
                return None

    def pack_title(self):
        try:
            tmp = []
            for line in self.title:
                tmp.append('* %s' % line)
            tmp.append('*')
            return '\n'.join(tmp)
        except TypeError:
            return "* A blank title.\n*"

    def pack_version(self):
        try:
            return "%5d%5d" % self.version
        except TypeError:
            return "    0    0"

    def _iter_normalize_card_no_comments(self):
        logger.debug('preprocessing body text, inline comments stripped')
        self.seek(0, 0)
        # filter comments, whitespace, blanklines
        iterable = ( line.strip() for line in self )
        iterable = ( line for line in iterable if line )
        iterable = ( line.lower() for line in iterable if not line.startswith(CC) )
        iterable = ( line.split(CC)[0] for line in iterable )
        iterable = ( line.strip() for line in iterable )
        # account for line continuation chars
        tmp = []
        for line in iterable:
            if line.endswith(KC):
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

    def _iter_normalize_card_yes_comments(self):
        logger.debug('preprocessing body text, keeping inline comments')
        self.seek(0, 0)
        # filters
        iterable = ( line.strip() for line in self )
        iterable = ( line for line in iterable if line )
        iterable = ( line for line in iterable if not line.startswith(CC) )
        #
        tmp_line = []
        tmp_comm = []
        for line in iterable:
            if CC in line:
                line, comment = line.split(CC, 1)
                line = line.strip()
                comment = comment.strip()
            else:
                comment = ""

            if line.endswith(KC):
                tmp_line.append(line[:-1].strip())
                tmp_comm.append(comment)
                continue

            tmp_line.append(line)
            tmp_comm.append(comment)
            line = ' '.join(tmp_line).strip()
            comment = ' '.join(tmp_comm).strip()
            if comment:
                line += ' %s %s' % (CC, comment)
            yield line.lower()
            tmp_line = []
            tmp_comm = []
        # flush
        if tmp_line or tmp_comm:
            line = ' '.join(tmp_line).strip()
            comment = ' '.join(tmp_comm).strip()
            if comment:
                line += ' %s %s' % (CC, comment)
            yield line.lower()

