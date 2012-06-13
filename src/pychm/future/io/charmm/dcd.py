"""DOCME

Documentation mostly done!

:TODO:
    :meth:`DCDFile.info` -- new method to pretty print header information
    :meth:`write_nparray` -- test, add error checking, right now it blindly
    writes
"""

from __future__ import division

__author__ = ("Frank C. Pickard <frank.pickard@nih.gov>")
__all__ = ["open_dcd"]

from array import array
import os
import struct
import warnings
import pdb

import numpy as np

from pychm.future.io.base import FortFile
from pychm.future.io.charmm.base import CharmmBin
from pychm.future.tools import rwprop


def open_dcd(fname, mode='r', buffering=None):
    """The public function responsible for mediating access to DCD file- like
    objects. Opens a file and returns a stream. If the file cannot be opened,
    an :exec:`IOError` is raised.

    Parameters
    ----------
    fname: a string representing the path to a CHARMM DCD file
    mode: an optional string that specifies the mode in which the file is
    opened. It defaults to 'r' which means open for reading in binary mode.
    Other values are shown in the table below:
    ========= =================================================================
    Character Meaning
    --------- -----------------------------------------------------------------
    'r'       open for reading (default)
    'w'       open for writing, truncating the file first
    'a'       open for writing, appending to the end of the file if it exists
    'x'       open for writing, noclobber, raises a value error if file exists
    '+'       open a disk for updating (reading and writing)
    'w+'      open for writing, random access, truncates the file to 0 bytes
    'r+'      open for reading and writing, random access, without truncation
    ========= =================================================================
    buffering: an optional integer used to set the buffering policy. Passing 0
    switches buffering off. Passing negative values, or 1 sets buffer to
    default size. Passing any other positive integer sets the buffer size in
    bytes.
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
    tmp = DCDFile(fname, mode=
                (reading and "r" or "") +
                (writing and "w" or "") +
                (appending and "a" or "") +
                (updating and "+" or ""),
                buffering=buffering)
    if reading or appending:
        tmp.read_header()
        tmp.seek_frame(0, whence=0)
    return tmp


class DCDFile(CharmmBin):
    """This class has no public constructor, please use :func:`open_dcd`
    instead."""
    def __init__(self, fname, mode='rb', buffering=None, endian='@',
                rec_head_prec='i', c_array_prec='i', xyz_prec='f',
                xtl_prec=None, d4_prec=None, q_prec=None):
        super(DCDFile, self).__init__(fname=fname, mode=mode,
            buffering=buffering, endian=endian, rec_head_prec=rec_head_prec)
        self.C_ARRAY_PREC = c_array_prec
        self.XTL_PREC = xtl_prec
        self.XYZ_PREC = xyz_prec
        self.D4_PREC = d4_prec
        self.Q_PREC = q_prec
        # header metadata
        self.dcdtype = None
        self.nframes = None
        self.npriv = None
        self.nsavc = None
        self.nsteps = None
        self.nsavv = None
        self.ndegfree = None
        self.nfix = None
        self.del_t = None
        # self.has_xtl # this is a property set from self.XTL_PREC
        # self.has_d4 # this is a property set from self.D4_PREC
        # self.has_q # this is a property set from self.Q_PREC
        self.inconsistent = None
        self.validated = None
        self.charmm_ver = None
        self.title = None
        self.natoms = None
        self.header_size = None
        # build frame data structure
        self._frame_dt = None
        self._frame_size = None
        # etc
        self._leftovers = []

    # Header reading, writing , importing, exporting ##########################
    def read_header(self):
        with FortFile(fname=self.name, mode='rb', buffering=self.buffer_size,
                    endian=self.ENDIAN, rec_head_prec=self.REC_HEAD_PREC) as fp:
            rec0 = fp.read_record()
            rec1 = fp.read_record()
            rec2 = fp.read_record()
        # read rec0
        self.dcdtype = ''.join(struct.unpack('cccc', rec0[:4])).lower()
        c_array = array(self.C_ARRAY_PREC, rec0[4:])
        self.nframes = c_array[0]
        self.npriv = c_array[1]
        self.nsavc = c_array[2]
        self.nsteps = c_array[3]
        self.nsavv = c_array[4]
        self.ndegfree = c_array[7]
        self.nfix = c_array[8]
        self.del_t = c_array[9]
        if c_array[10] == 1:
            self._xtl_prec = 'd'
        if c_array[11] == 1:
            self._d4_prec = 'f'
        if c_array[12] == 1:
            self._q_prec = 'f'
        self.inconsistent = bool(c_array[13])
        self.validated = c_array[18] == 1
        self.charmm_ver = c_array[19]
        # read rec1
        self.title = ''.join(list(array('c', rec1)))
        # read rec2
        if len(rec2) == self.C_ARRAY_BLEN:
            self.natoms = struct.unpack(self.C_ARRAY_PREC, rec2)[0]
        else:
            raise IOError("Error when parsing natoms record of dcd header")
        # check first four chars in rec0 for gibberish
        if self.dcdtype not in ['cord', 'veld']:
            warnings.warn("unexpected dcdtype: %s" %self.dcdtype)
        # calc header_size
        self.header_size = len(rec0) + len(rec1) + len(rec2) + 24
        # compile frame info
        self._frame_dt = self.compile_npdt()
        self._frame_size = self._frame_dt.itemsize

    def write_header(self):
        """The class' current values for attributes defined by a CHARMM DCD
        header are read from the class, mapped to a binary format, and written
        to disk.
        """
        # length of header less the title
        hlen_minus_title = 6 * self.REC_HEAD_BLEN + 4 + 21 * self.C_ARRAY_BLEN
        if '+' in self.mode:
            # TODO if len(self.title) is too long, truncate it
            # TODO if len(self.title) is too short, pad with whitespace
            return
        # build rec0
        rec0 = list(self.dcdtype[:4].upper())
        rec0.append(int(self.nframes))      # c_array 1
        rec0.append(int(self.npriv))
        rec0.append(int(self.nsavc))
        rec0.append(int(self.nsteps))
        rec0.append(int(self.nsavv))        # c_array 5
        rec0.append(0)
        rec0.append(0)
        rec0.append(int(self.ndegfree))
        rec0.append(int(self.nfix))
        rec0.append(int(self.del_t))        # c_array 10
        rec0.append(int(self.has_xtl))
        rec0.append(int(self.has_d4))
        rec0.append(int(self.has_q))
        rec0.append(int(self.inconsistent))
        rec0.append(0)                      # c_array 15
        rec0.append(0)
        rec0.append(0)
        rec0.append(0)
        rec0.append(int(self.validated))
        rec0.append(int(self.charmm_ver))   # c_array 20
        rec0_formatting = '4c 20%s' % self.C_ARRAY_PREC
        rec0 = struct.pack(rec0_formatting, *rec0)
        # build rec1
        rec1 = list(self.title)
        rec1_formatting = '%dc' % len(self.title)
        rec1 = struct.pack(rec1_formatting, *rec1)
        # build rec2
        rec2 = int(self.natoms)
        rec2_formatting = self.C_ARRAY_PREC
        rec2 = struct.pack(rec2_formatting, rec2)
        # write
        self.fp.seek(0, 0)
        self.write_record(rec0)
        self.write_record(rec1)
        self.write_record(rec2)
        # verify and set frame vars
        self.read_header()

    def import_header(self, arg):
        tmp = self.export_header()
        try:
            for k, v in arg.iteritems():
                setattr(self, k, v)
        except:
            warnings.warn("Error importing values, reverting to previous")
            for k, v in tmp.iteritems():
                setattr(self, k, v)

    def export_header(self):
        tmp = {
            'C_ARRAY_PREC': self.C_ARRAY_PREC,
            'XTL_PREC': self.XTL_PREC,
            'XYZ_PREC': self.XYZ_PREC,
            'D4_PREC': self.D4_PREC,
            'Q_PREC': self.Q_PREC,
            'dcdtype': self.dcdtype,
            'nframes': self.nframes,
            'npriv': self.npriv,
            'nsavc': self.nsavc,
            'nsteps': self.nsteps,
            'nsavv': self.nsavv,
            'ndegfree': self.ndegfree,
            'nfix': self.nfix,
            'del_t': self.del_t,
            'inconsistent': self.inconsistent,
            'validated': self.validated,
            'charmm_ver': self.charmm_ver,
            'title': self.title,
            'natoms': self.natoms,
            'header_size': self.header_size
            }
        return tmp

    def pprint_header(self):
        pass

    # Highest level methods ###################################################
    def iter_frame(self, begin=0, end=None, stride=None):
        """Creates a generator over the frames of the DCDFile. Uses the
        (hopefully) familiar start:stop:stride syntax of :class:`slice`s.
        This generator returns raw binary data representing one MD frame.
        """
        if end is not None and not isinstance(end, int):
            raise TypeError("invalid type for end: %r" % end)
        if stride is not None and not isinstance(stride, int):
            raise TypeError("invalid type for stride: %r" % end)
        if stride is not None and stride < 1:
            raise ValueError("invalid value for stride: %r" % end)
        #
        self.seek_frame(begin, whence=0)
        if end is None:
            if stride is None:
                while 1:
                    yield self.next()
            else:
                while 1:
                    yield self.next()
                    self.seek_frame(stride-1, whence=1)
        else:
            if begin > end:
                warnings.warn("begin > end, this will be an empty iterator")
            if stride is None:
                while 1:
                    if begin > end:
                        raise StopIteration
                    yield self.next()
                    begin += 1
            else:
                while 1:
                    if begin > end:
                        raise StopIteration
                    yield self.next()
                    self.seek_frame(stride-1, whence=1)
                    begin += stride

    def iter_nparray(self, begin=0, end=None, stride=None):
        """Creates a generator over the frames of the DCDFile. Uses the
        (hopefully) familiar start:stop:stride syntax of :class:`slice`s.  This
        generator returns a :class:`numpy.ndarray` object representing one MD
        frame.
        """
        if end is not None and not isinstance(end, int):
            raise TypeError("invalid type for end: %r" % end)
        if stride is not None and not isinstance(stride, int):
            raise TypeError("invalid type for stride: %r" % end)
        if stride is not None and stride < 1:
            raise ValueError("invalid value for stride: %r" % end)
        #
        self.seek_frame(begin, whence=0)
        if end is None:
            if stride is None:
                try:
                    while 1:
                        yield np.fromfile(self.fp, dtype=self.frame_dt, count=1)
                except MemoryError:
                    raise StopIteration
            else:
                try:
                    while 1:
                        yield np.fromfile(self.fp, dtype=self.frame_dt, count=1)
                        self.seek_frame(stride-1, whence=1)
                except MemoryError:
                    raise StopIteration
        else:
            if begin > end:
                warnings.warn("begin > end, this will be an empty iterator")
            if stride is None:
                try:
                    while 1:
                        if begin > end:
                            raise StopIteration
                        yield np.fromfile(self.fp, dtype=self.frame_dt, count=1)
                        begin += 1
                except MemoryError:
                    raise StopIteration
            else:
                try:
                    while 1:
                        if begin > end:
                            raise StopIteration
                        yield np.fromfile(self.fp, dtype=self.frame_dt, count=1)
                        self.seek_frame(stride-1, whence=1)
                        begin += stride
                except MemoryError:
                    raise StopIteration

    def get_massive_dump(self):
        """Returns a large 3-D :class:`numpy.ndarray` containing all of the
        frame data for the entire trajectory.
        """
        self.seek_frame(0, whence=0)
        tmp = np.fromfile(self.fp, dtype=self.frame_dt, count=-1)
        self.readline()
        return tmp

    def read_frame(self):
        """Reads and returns a full frame in binary format, an empty binary
        string (if there are no more frames to be read) or reads a partial
        frame into the :attr:`leftovers` attribute.
        """
        return self.readline()

    def read_nparray(self):
        """Reads and returns a full frame formatted as a
        :class:`numpy.ndarray`.
        """
        return np.fromfile(self.fp, dtype=self.frame_dt, count=1)

    def write_nparray(self, nparray, order='C'):
        """Takes an appropriately structed :class:`numpy.ndarray`, converts it
        to binary format and writes it to disk.

        *Untested*

        :TODO: Error checking, right now it blindly writes
        """
        if not isinstance(nparray, np.ndarray):
            raise TypeError("invalid nparray")
        self.fp.write(nparray.tostring(order))

    @property
    def leftovers(self):
        return self._leftovers

    # Frame (meta)data ########################################################
    def compile_npdt(self):
        """Uses the availible precision specifications to build a
        :class:`numpy.dtype` object. This is then used in turn to convert
        between :class:`numpy.ndarray` objects and binary formatting for
        writing to disk.
        """
        tmp = []
        if self.has_xtl:
            tmp.append(('xtl', self.ENDIAN + self.XTL_PREC + self.XTL_BLEN, 6))
        tmp.append(('x', '%s%s%d' % (self.ENDIAN, self.XYZ_PREC, self.XYZ_BLEN),
                    self.natoms))
        tmp.append(('y', '%s%s%d' % (self.ENDIAN, self.XYZ_PREC, self.XYZ_BLEN),
                    self.natoms))
        tmp.append(('z', '%s%s%d' % (self.ENDIAN, self.XYZ_PREC, self.XYZ_BLEN),
                    self.natoms))
        if self.has_d4:
            tmp.append(('d4', '%s%s%d' % (self.ENDIAN, self.D4_PREC, self.D4_BLEN),
                        self.natoms))
        if self.has_q:
            tmp.append(('q', '%s%s%d' % (self.ENDIAN, self.Q_PREC, self.Q_BLEN),
                        self.natoms))
        # add record markers
        tmp2 = []
        if True:
            for i, t in enumerate(tmp):
                tmp2.append(('pad%d' % (i*2), '%s%s%d' % (self.ENDIAN,
                            self.REC_HEAD_PREC, self.REC_HEAD_BLEN), (1,)))
                tmp2.append(t)
                tmp2.append(('pad%d' % (i*2+1), '%s%s%d' % (self.ENDIAN,
                            self.REC_HEAD_PREC, self.REC_HEAD_BLEN), (1,)))
            return np.dtype(tmp2)
        return np.dtype(tmp)

    @property
    def frame_dt(self):
        return self._frame_dt

    @property
    def frame_size(self):
        return self._frame_size

    # Wrapper to python file API ##############################################
    def readline(self):
        """Reads and returns a full frame in binary format, an empty binary
        string (if there are no more frames to be read) or reads a partial
        frame into the :attr:`leftovers` attribute.
        """
        tmp = self.fp.read(self.frame_size)
        if len(tmp) == self.frame_size:
            return tmp
        elif len(tmp) == 0:
            return b''
        else:
            warnings.warn("Did not read enough data, wanted %r bytes, received %r bytes." % (self.frame_size, len(tmp)))
            warnings.warn("Dumping fractional frame to leftovers.")
            self._leftovers.append(tmp)

    def seek_frame(self, offset, whence=0):
        """Seek the dcd file using a byframe `offset` instead of a
        bytewise `offset`. If you want to seek the file in a bytewise
        manner, use the underlying file pointer directly.

        ====== =============================================================
        whence meaning
        ------ -------------------------------------------------------------
        0       seek in relation to the first frame (skips header) [default]
        -1      seek to the begining of the file, ignores `offset`
        1       seek in relation to current position
        2       seek in relation to last frame
        ====== =============================================================

        To seek to the very first frame in the file, use...
        >>> taco = open_dcd("some_file")
        >>> taco.seek_frame(0, 0)

        To seek to the 8109th byte of the file, use...
        >>> taco = open_dcd("some_file")
        >>> taco.fp.seek(8109, 0)

        """
        if whence == -1:
            self.fp.seek(0, 0)
        elif whence == 0:
            self.fp.seek(offset * self.frame_size + self.header_size, whence)
        else:
            self.fp.seek(offset * self.frame_size, whence)

    def tell_frame(self):
        """Returns the current position of the filepointer as a function
        of the frame. If an integer is returned the filepointer is between
        frames, if a float is returned the filepointer is in a frame.
        """
        tmp = (self.fp.tell() - self.header_size) / self.frame_size
        if tmp == int(tmp):
            return int(tmp)
        else:
            return tmp

    def truncate_frame(self, offset):
        """Truncates the file after the specified frame. If you want to
        truncate in a bytewise manner, the underlying filepointer should be
        used.  Offset must be an integer.
        """
        if offset == int(offset):
            self.fp.truncate(offset * self.frame_size + self.header_size)
        else:
            raise ValueError("offset must be an integer")

    # Byte Bookeeping #########################################################
    @rwprop
    def C_ARRAY_PREC():
        def fget(self):
            return self._c_array_prec
        def fset(self, value):
            if value is None or value in 'hilq':
                self._c_array_prec = value
            else:
                raise ValueError("Incorrect precision specified: %r" % value)
        doc = "Possible precision values are 'h', 'i', 'l' or 'q'"
        return locals()

    @property
    def C_ARRAY_BLEN(self):
        return self._byte_dict[self._c_array_prec]

    @rwprop
    def XTL_PREC():
        def fget(self):
            return self._xtl_prec
        def fset(self, value):
            if value is None or value in 'fd':
                self._xtl_prec = value
            else:
                raise ValueError("Incorrect precision specified: %r" % value)
        doc = "Possible precision values are 'f' or 'd'"
        return locals()

    @property
    def XTL_BLEN(self):
        return self._byte_dict[self._xtl_prec]

    @property
    def has_xtl(self):
        return not self._xtl_prec is None

    @rwprop
    def XYZ_PREC():
        def fget(self):
            return self._xyz_prec
        def fset(self, value):
            if value is None or value in 'fd':
                self._xyz_prec = value
            else:
                raise ValueError("Incorrect precision specified: %r" % value)
        doc = "Possible precision values are 'f' or 'd'"
        return locals()

    @property
    def XYZ_BLEN(self):
        return self._byte_dict[self._xyz_prec]

    @property
    def has_xyz(self):
        return not self._xyz_prec is None

    @rwprop
    def D4_PREC():
        def fget(self):
            return self._d4_prec
        def fset(self, value):
            if value is None or value in 'fd':
                self._d4_prec = value
            else:
                raise ValueError("Incorrect precision specified: %r" % value)
        doc = "Possible precision values are 'f' or 'd'"
        return locals()

    @property
    def D4_BLEN(self):
        return self._byte_dict[self._d4_prec]

    @property
    def has_d4(self):
        return not self._d4_prec is None

    @rwprop
    def Q_PREC():
        def fget(self):
            return self._q_prec
        def fset(self, value):
            if value is None or value in 'fd':
                self._q_prec = value
            else:
                raise ValueError("Incorrect precision specified: %r" % value)
        doc = "Possible precision values are 'f' or 'd'"
        return locals()

    @property
    def Q_BLEN(self):
        return self._byte_dict[self._q_prec]

    @property
    def has_q(self):
        return not self._q_prec is None
