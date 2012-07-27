"""
"""

from __future__ import division
__author__ = ("Frank C. Pickard <frank.pickard@nih.gov>")
__all__ = ["open_fort"]

from io import DEFAULT_BUFFER_SIZE
import os
import struct
import warnings
import logging

import numpy as np

from pychm.future.tools import rwprop


logger = logging.getLogger('pychm.io.base')


def get_buffer_size(arg=None):
    """Function to calculate buffer size if it is not explicitly specified."""
    if arg is None or arg < 0 or arg==1:
        try:
            return os.stat(os.curdir).st_blksize
        except (os.error, AttributeError):
            return DEFAULT_BUFFER_SIZE
    else:
        return arg


def open_fort(fname, mode='r', buffering=None):
    """The public function responsible for mediating access to the Fort file-
    like objects. Opens a file and returns a stream. If the file cannot be
    opened, an IOError is raised.

    Parameters
    ----------
    fname : a string representing the path to a fortran binary file
    mode : an optional string that specifies the mode in which the file is
        opened. It defaults to 'r' which means open for reading in binary
        mode. Other values are shown in the table below:
    ========= =================================================================
    Character Meaning
    --------- -----------------------------------------------------------------
    'r'       open for reading (default)
    'w'       open for writing, truncating the file first
    'a'       open for writing, appending to the end of the file if it exists
    'x'       open for writing, noclobber, raises a value error if file exists
    '+'       open a disk for updating (reading and writing)
    'w+'     open for writing, random access, truncates the file to 0 bytes
    'r+'     open for reading and writing, random access, without truncation
    ========= =================================================================
    buffering: an optional integer used to set the buffering policy. Passing 0
        switches buffering off. Passing negative values, or 1 sets buffer
        to default size. Passing any other positive integer sets the buffer
        size in bytes.
    """
    if not isinstance(fname, basestring):
        raise TypeError("Invalid fname: %r" % fname)
    if not isinstance(mode, basestring):
        raise TypeError("Invalid mode: %r" % mode)
    if buffering is not None and not isinstance(buffering, int):
        raise TypeError("Invalid buffering: %r" % buffering)
    # parse modes
    modes = set(mode)
    if modes - set("arw+xb") or len(mode) > len(modes):
        raise ValueError("invalid mode: %r" % mode)
    reading = "r" in modes
    writing = "w" in modes or "x" in modes
    appending = "a" in modes
    updating = "+" in modes
    if reading + writing + appending > 1:
        raise ValueError("Invalid mode: %r" % mode)
    if not (reading or writing or appending):
        raise ValueError("must have exactly one read/write/append mode")
    if "x" in modes and os.path.isfile(fname):
        raise ValueError("you may not set `mode=x` for existing files")
    # instantiate!
    return FortFile(fname=fname, mode=
                    (reading and "r" or "") +
                    (writing and "w" or "") +
                    (appending and "a" or "") +
                    (updating and "+" or ""),
                    buffering=buffering)


class File(object):
    """
    """
    def __init__(self, fname, mode='r', buffering=None):
        super(File, self).__init__()
        self._buffer_size = get_buffer_size(buffering)
        self.fp = open(name=fname, mode=mode, buffering=self._buffer_size)
        logger.debug("opening file: %s" % self.fp.name)
        logger.debug("mode set to: %s" % self.fp.mode)
        logger.debug("buffering set to: %d" % self._buffer_size)

    # Wrapper to python file API ##############################################
    def close(self):
        self.flush()
        self.fp.close()

    def flush(self):
        self.fp.flush()
        os.fsync(self.fileno())

    def fileno(self):
        return self.fp.fileno()

    def isatty(self):
        return self.fp.isatty()

    def next(self):
        tmp = self.readline()
        if tmp:
            return tmp
        else:
            raise StopIteration

    def read(self, n=-1):
        return self.fp.read(n)

    def readline(self):
        return self.fp.readline()

    def readlines(self, n=-1):
        if n < 0:
            return [ line for line in self ]
        tmp = []
        for i in xrange(n):
            tmp.append(self.next())
        return tmp

    def seek(self, offset, whence=0):
        self.fp.seek(offset, whence)

    def tell(self):
        return self.fp.tell()

    def truncate(self, offset):
        self.fp.truncate(offset)

    def write(self, s):
        self.fp.write(s)

    def writelines(self, iterable):
        self.fp.writelines(iterable)

    # More python file API definitions ########################################
    @property
    def closed(self):
        return self.fp.closed

    @property
    def encoding(self):
        return self.fp.encoding

    @property
    def errors(self):
        return self.fp.errors

    @property
    def mode(self):
        return self.fp.mode

    @property
    def name(self):
        return self.fp.name

    @property
    def newlines(self):
        return self.fp.newlines

    @rwprop
    def softspace():
        def fget(self):
            return self.fp.softspace
        def fset(self, v):
            self.fp.softspace = v
        doc = "flag indicating the a space needs to be printed; used by print"
        return locals()

    @property
    def buffer_size(self):
        return self._buffer_size

    # Private Utility Methods #################################################
    def _checkClosed(self):
        if self.closed:
            raise ValueError("I/O operation on closed file.")

    # Special Methods #########################################################
    # #### Iterator Protocol
    def __iter__(self):
        self._checkClosed()
        return self

    def __repr__(self):
        try:
            return "%s(%d, %s)" % (self.__class__.__name__, self.fileno(), self.mode)
        except ValueError:
            return "%s(CLOSED, %s)" % (self.__class__.__name__, self.mode)

    # #### Context Manager Protocol
    def __enter__(self):
        self._checkClosed()
        return self

    def __exit__(self, *args):
        self.close()


class BinFile(File):
    """
    """
    def __init__(self, fname, mode='rb', buffering=None, endian='@'):
        if 'b' not in mode:
            mode += 'b'
        super(BinFile, self).__init__(fname=fname, mode=mode, buffering=buffering)
        self.ENDIAN = endian
        logger.debug("ENDIAN set to %s" % self._endian)

    def _read_exactly(self, nbytes):
        """Read exactly nbytes, raise an error otherwise."""
        tmp = self.read(nbytes)
        if len(tmp) == nbytes:
            return tmp
        else:
            raise IOError("Did not read enough data, wanted %r bytes, received %r bytes." % (nbytes, len(tmp)))

    # These methods have indeterminant meaning for bin files ##################
    def readline(self):
        raise NotImplementedError

    def writelines(self, iterable):
        raise NotImplementedError

    # Byte Bookeeping #########################################################
    @rwprop
    def ENDIAN():
        def fget(self):
            return self._endian
        def fset(self, value):
            if value is None or value in '<>@=!':
                if value is None or value in '@=':
                    self._endian = '='
                elif value in '>!':
                    self._endian = '>'
                else:
                    self._endian = '<'
            else:
                raise ValueError("Incorrect 'endian' specifier: %r" % value)
        doc = "Possible endian values are '<', '>', '@', '=' or '!'"
        return locals()

    _byte_dict = {
        None: None,
        'h':2,
        'i':4,
        'l':4,
        'q':8,
        'f':4,
        'd':8}

    _np_byte_dict = {
        None: None,
        'h': np.int16,
        'i': np.int32,
        'l': np.int32,
        'q': np.int64,
        'f': np.float32,
        'd': np.float64}


class FortFile(BinFile):
    def __init__(self, fname, mode='rb', buffering=None,
                endian='@', rec_head_prec='i'):
        """The `func:open_fort` should be used in lieu of this constructor.
        """
        super(FortFile, self).__init__(fname=fname, mode=mode,
                                    buffering=buffering, endian=endian)
        self.REC_HEAD_PREC = rec_head_prec
        logger.debug("record marker precision set to %s" % self._rec_head_prec)

    def _read_chk(self):
        return struct.unpack(self.ENDIAN + self.REC_HEAD_PREC,
                            self._read_exactly(self.REC_HEAD_BLEN))[0]

    def _write_chk(self, nbytes):
        self.write(struct.pack(self.ENDIAN + self.REC_HEAD_PREC, nbytes))

    def readline(self):
        """Read a single fortran record."""
        return self.read_record()

    def read_record(self):
        """Read a single fortran record."""
        l = self._read_chk()
        data = self._read_exactly(l)
        chk = self._read_chk()
        if chk != l:
            raise IOError("Record markers do not match.")
        return data

    def write_record(self, s):
        """Write a single fortran record, containing bytestring `s`."""
        l = len(s)
        self._write_chk(l)
        self.write(s)
        self._write_chk(l)
        self.flush()

    def writelines(self, iterable):
        """Write an iterable of bytestrings to multiple fortran records."""
        for i in iterable:
            self.write_record(i)

    def read_array(self, prec):
        """Read a fortran record, decoded as an array of floats or integers,
        as specified by prec.
        ========= =========================
        Character Meaning
        --------- -------------------------
        'f'       float, single precision
        'd'       float, double precision
        'h'       integer, half precision
        'i'       integer, single precision
        'l'       integer, single precision
        'q'       integer, double precision
        ========= =========================
        """
        if prec not in 'fdhilq':
            raise ValueError("Invalid precision specified: %r" % prec)
        data = self.read_record()
        num = len(data) / self._byte_dict[prec]
        return np.array(struct.unpack("%s%d%s" % (self.ENDIAN, num, prec), data),
                        dtype=np.dtype(prec))

    def read_reals(self, prec='f'):
        """Read a fortran record, decoded as an array of floats with single
        (prec='f') or double (prec='d') precision."""
        if prec not in 'fd':
            raise ValueError("Invalid precision specified: %r" % prec)
        return self.read_array(prec)

    def read_ints(self, prec='i'):
        """Read a fortran record, decoded as an array of ints with half
        (prec='h'), single (prec='i' or 'l') or double (prec='q') precision."""
        if prec not in 'hilq':
            raise ValueError("Invalid precision specified: %r" % prec)
        return self.read_array(prec)

    def write_array(self, arraylike, prec):
        """Write a fortran record, encoded as an array of floats or integers,
        as specified by prec.
        ========= =========================
        Character Meaning
        --------- -------------------------
        'f'       float, single precision
        'd'       float, double precision
        'h'       integer, half precision
        'i'       integer, single precision
        'l'       integer, single precision
        'q'       integer, double precision
        ========= =========================
        """
        if prec not in 'fdhilq':
            raise ValueError("Invalid precision specified: %r" % prec)
        l = len(arraylike) * self._byte_dict[prec]
        self._write_chk(l)
        self.write(struct.pack("%s%d%s" % (self.ENDIAN, l, prec), *arraylike))
        self._write_chk(l)

    def write_reals(self, arraylike, prec='f'):
        """Write a fortran record, encoded as an array of floats with single
        (prec='f') or double (prec='d') precision."""
        if prec not in 'fd':
            raise ValueError("Invalid precision specified: %r" % prec)
        self.write_array(arraylike, prec)

    def write_ints(self, arraylike, prec='i'):
        """Write a fortran record, encoded as an array of ints with half
        (pref='h'), single (prec='i' or 'l') or double (prec='q') precision."""
        if prec not in 'hilq':
            raise ValueError("Invalid precision specified: %r" % prec)
        self.write_array(arraylike, prec)

    # Byte Bookeeping #########################################################
    @rwprop
    def REC_HEAD_PREC():
        def fget(self):
            return self._rec_head_prec
        def fset(self, value):
            if value in 'hilq':
                self._rec_head_prec = value
            else:
                raise ValueError("Incorrect precision specified: %r" % value)
        doc = "Possible precision values are 'h', 'i', 'l' or 'q'"
        return locals()

    @property
    def REC_HEAD_BLEN(self):
        return self._byte_dict[self._rec_head_prec]


class TextFile(File):
    """
    """
    def __init__(self, fname, mode='r', buffering=None):
        if 'b' in mode:
            raise ValueError("Invalid mode, opening TextFile using 'b'.")
        super(TextFile, self).__init__(fname=fname, mode=mode,
                                    buffering=buffering)
