"""
TODO:
    func:open_dcd() -- finish writing doc_string
    func:open_dcd() -- validate byte related variables
    func:guess_natoms() -- untested
    guess_xtl() -- untested
    read_header() -- rewrite as classs for mod via subclass
    pack_header() -- rewrite as class for mod via sublclass

    write:
        DCDFile.info() -- a method to pretty print header information
"""

from __future__ import division
############################################ This still breaks numpy
## from __future__ import unicode_literals # this is a py3k
############################################ compatibility problem.
__author__ = ("Frank C. Pickard <frank.pickard@nih.gov>")
__all__ = ["open_dcd", "read_header", "pack_header"]

from array import array
import os
import struct
import warnings

import numpy as np

from pychm.future.io.base import open_fort
from pychm.future.io.charmm.base import CharmmBin
from pychm.future.tools import rwprop


def open_dcd(fname, mode='r', buffering=None, header=True, has_rec=True,
            endian='@', rec_head_prec='i', c_array_prec='i', xtl_prec='d',
            xyz_prec='f', d4_prec='f', q_prec='f', **kwargs):
    """The public function responsible for mediating access to DCD file-
    like objects. Opens a file and returns a stream. If the file cannot be
    opened, an IOError is raised.

    Parameters
    ----------
    fname : a string representing the path to a CHARMM DCD file
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
    endian: an optional character specifying the endian-ness of the file.
        Possible values are '>', '<', '@', '=' or '!'. Consult :mod:`struct`
        for more details. Defaults to '@', native byte order.
    rec_head_prec: an optional character specifying the fortran record header's
        precision. Possible values are 'h', 'i', 'l' or 'q'. Consult
        :mod:`struct` for more details. Defaults to 'i', single precision
        integers.
    TODO Finish the doc string...
    """
    if not isinstance(fname, basestring):
        raise TypeError("Invalid fname: %r" % fname)
    if not isinstance(mode, basestring):
        raise TypeError("Invalid mode: %r" % mode)
    if buffering is not None and not isinstance(buffering, int):
        raise TypeError("Invalid buffering: %r" % buffering)
    if not isinstance(header, bool):
        raise TypeError("Invalid header: %r" % header)
    if not isinstance(has_rec, bool):
        raise TypeError("Invalid has_rec: %r" % has_rec)

    # #################################################################
    # TODO
    # We should probably validate all the byte related variables here.
    # While I am lazy, these variables do get validated eventually
    # by descriptors in the TRJFile class.
    # #################################################################

    # parse modes
    modes = set(mode)
    if modes - set("arw+x") or len(mode) > len(modes):
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
    # parse has_rec
    if has_rec:
        # kwargs can explicitly overwrite parsed metadata from header
        if header:
            tmp = read_header(fname, endian, rec_head_prec, c_array_prec)
            tmp.update(kwargs)
            kwargs = tmp
            del tmp
        else:
            # has_rec = True && header = False
            rec_head_prec = None
            c_array_prec = None
            kwargs['header_size'] = 0
            pass
    else:
        if header:
            raise ValueError("reader does not yet support files with a header and without record markers")
        else:
            # has_rec = False && header = False
            rec_head_prec = None
            c_array_prec = None
            kwargs['header_size'] = 0
    # any metadata (if it exists) should be in kwargs by now
    # parse kwargs
    if 'dcdtype' not in kwargs:
        kwargs['dcdtype'] = 'taco'
    # figure out natoms
    if 'natoms' not in kwargs and 'ndegfree' in kwargs:
        kwargs['natoms'] = kwargs['ndegfree'] / 3
    elif 'natoms' in kwargs and 'ndegfree' not in kwargs:
        kwargs['ndegfree'] = kwargs['natoms'] * 3
    elif 'natoms' not in kwargs and 'ndegfree' not in kwargs:
        if header:
            # if we have a header, and we dont know natoms, this is bad
            raise ValueError("no natom info found in header, aborting...")
        if not has_rec:
            # if we have no record markers, we cant possibly guess natoms
            raise ValueError("no natom information provided, aborting...")
        warnings.warn("`natoms` not specified, trying to autodetect", UserWarning)
        # has_rec = True && header = False
        try:
            kwargs['natoms'] = guess_natoms(fname, endian, rec_head_prec,
                                        xyz_prec)
        except (ValueError, IOError):
            raise ValueError("unable to guess natoms, aborting...")
        kwargs['ndegfree'] = kwargs['natoms'] * 3
        warnings.warn("detected natoms = %d" % kwargs['natoms'], UserWarning)
    # figure out has_xtl
    if 'has_xtl' not in kwargs:
        if header:
            # if we have a header, and we dont know has_xtl, this is bad
            raise ValueError("no has_xtl information found in header, aborting...")
        if not has_rec:
            # if we have no record markers, we cant possible guess has_xtl
            # lets assume false (most common use case), and throw a warning
            warnings.warn("has_xtl not specified, assuming False", UserWarning)
            kwargs['has_xtl'] = False
        else:
            # has_rec = True && header = False
            # try to guess has_xtl, if we fail assume false and throw a warning
            warnings.warn("has_xtl not specified, trying to autodetect", UserWarning)
            try:
                kwargs['has_xtl'] = guess_xtl(fname, endian, rec_head_prec,
                                            xtl_prec)
            except (ValueError, IOError):
                warnings.warn("unable to guess has_xtl, assuming False", UserWarning)
            else:
                warnings.warn("detected has_xtl = %s" % kwargs['has_xtl'], UserWarning)
    if 'has_d4' not in kwargs:
        warnings.warn("has_d4 not specified, assuming False", UserWarning)
        kwargs['has_d4'] = False
    if 'validated' not in kwargs:
        warnings.warn("validated not specified, assuming False", UserWarning)
        kwargs['validated'] = False
    # other kwarg processing
    kwargs['has_rec'] = has_rec
    kwargs['ENDIAN'] = endian
    kwargs['REC_HEAD_PREC'] = rec_head_prec
    kwargs['C_ARRAY_PREC'] = c_array_prec
    kwargs['XTL_PREC'] = xtl_prec
    kwargs['XYZ_PREC'] = xyz_prec
    kwargs['D4_PREC'] = d4_prec
    kwargs['Q_PREC'] = q_prec
    if 'has_q' not in kwargs:
        kwargs['has_q'] = False
    # instantiate!
    return DCDFile(fname,
                    (reading and "r" or "") +
                    (writing and "w" or "") +
                    (appending and "a" or "") +
                    (updating and "+" or ""),
                    buffering, **kwargs)


def read_header(fname, endian='@', rec_head_prec='i', c_array_prec='i'):
    """Attempts to parse out metadata from `fname`s header information. Returns
    a dictionary. This function should *only* be called on a charmm trajectory
    *with* record markers and *with* a trajectory header.

    Parameters
    ----------
    fname : a string representing the path to a fortran binary file
    endian: an optional character specifying the endian-ness of the file.
        Possible values are '>', '<', '@', '=' or '!'. Consult :mod:`struct`
        for more details. Defaults to '@', native byte order.
    rec_head_prec: an optional character specifying the fortran record header's
        precision. Possible values are 'h', 'i', 'l' or 'q'. Consult
        :mod:`struct` for more details. Defaults to 'i', single precision
        integers.
    c_array_prec: an optional character specifying the control array's
        precision. Possible values are 'h', 'i', 'l' or 'q'. Consult
        :mod:`struct` for more details. Defaults to 'i', single precision
        integers.
    """
    if not isinstance(fname, basestring):
        raise TypeError("Invalid mode: %r" % fname)
    if not isinstance(endian, basestring):
        raise TypeError("Invalid endian: %r" % endian)
    if not isinstance(rec_head_prec, basestring):
        raise TypeError("Invalid rec_head_prec: %r" % rec_head_prec)
    if not isinstance(c_array_prec, basestring):
        raise TypeError("Invalid c_array_prec: %r" % c_array_prec)
    #
    if not c_array_prec in 'hilq':
        raise ValueError("Invalid precision specified: %r" % c_array_prec)
    # defaults
    tmp = {
        'dcdtype':None,
        'nframes':None,
        'npriv':None,           # number of previous dynamics steps
        'nsav':None,
        'nsteps':None,
        'ndegfree':None,
        'has_xtl':None,
        'has_d4':None,
        'validated':None,       # NOT OFFICIAL -- Flag if we know header info matches data
        'title':'',
        'natoms':None,
        'header_size':0
        }
    # read file
    with open_fort(fname, mode='r', buffering=None, endian=endian,
                rec_head_prec=rec_head_prec) as fp:
        control_array = fp.read_record()
        title = fp.read_record()
        natoms = fp.read_record()
    # parse header
    ## rec0 - control Array
    tmp['dcdtype'] = ''.join(struct.unpack('cccc', control_array[:4])).lower()
    c_array = array(c_array_prec, control_array[4:])
    tmp['nframes'] = c_array[0]
    tmp['npriv'] = c_array[1]
    tmp['nsav'] = c_array[2]
    tmp['nsteps'] = c_array[3]
    tmp['ndegfree'] = c_array[7]
    tmp['has_xtl'] = c_array[10] == 1
    tmp['has_d4'] = c_array[11] == 1
    tmp['validated'] = c_array[18] == 1
    ## rec1 - title
    tmp['title'] = array('c', title)
    tmp['title'] = [ char for char in tmp['title'] if True ]
    tmp['title'] = ''.join(tmp['title'])
    ## rec2 - natoms
    if len(natoms) == 4:
        tmp['natoms'] = struct.unpack(c_array_prec, natoms)[0]
    else:
        raise IOError("Something went wrong when reading natoms")
    # check first four chars in rec0 for gibberish
    if tmp['dcdtype'] not in ['cord', 'taco']:
        raise ValueError("This function expected a dcd file with a header, something has gone wrong.")
    # calc header_size and return
    tmp['header_size'] = len(control_array) + len(title) + len(natoms) + 24
    return tmp


def pack_header(header_size=None, **kwargs):
    """DOCME
    """
    if header_size is None:
        warnings.warn("Not specifying a header_size may corrupt your traj file if in update mode")
    else:
        if not isinstance(header_size, int):
            raise TypeError("invalid header_size: %r" % header_size)
    #
    defaults = {
        'dcdtype':'taco',
        'nframes':0,
        'npriv':0,
        'nsav':0,
        'nsteps':0,
        'ndegfree':0,
        'has_xtl':False,
        'has_4d':False,
        'validated':False,
        'title':' ' * 78,
        'natoms':0
        }
    defaults.update(kwargs)
    kwargs = defaults
    del defaults
    # guess a reasonable value for the header_size if its not present
    if header_size is None:
        header_size = len(kwargs['title']) + 112 # (4+84+4) + (4+len(title)+4) + (4+4+4)
    # check if title is too long
    if len(kwargs['title']) + 112 > header_size:
        raise ValueError('title is too long for header: %d vs %d' % (len(kwargs['title']), header_size))
    # Pad title if too short
    while len(kwargs['title']) + 112 < header_size:
        kwargs['title'] += ' '
    ##
    # First record
    rec0 = []
    rec0.append(84)
    rec0.extend(list(kwargs['dcdtype'][:4].upper()))
    rec0.append(int(kwargs['nframes']))             # carray 1
    rec0.append(int(kwargs['npriv']))
    rec0.append(int(kwargs['nsav']))
    rec0.append(int(kwargs['nsteps']))
    rec0.append(0)                                  # carray 5
    rec0.append(0)
    rec0.append(0)
    rec0.append(int(kwargs['ndegfree']))
    rec0.append(0)
    rec0.append(0)                                  # carray 10
    rec0.append(int(kwargs['has_xtl']))
    rec0.append(int(kwargs['has_4d']))
    rec0.append(0)
    rec0.append(0)
    rec0.append(0)                                  # carray 15
    rec0.append(0)
    rec0.append(0)
    rec0.append(0)
    rec0.append(int(kwargs['validated']))
    rec0.append(0)                                  # carray 20
    rec0.append(84)
    rec0_string = 'i 4c 20i i'
    # Second record
    rec1 = []
    rec1.append(header_size - 112)
    rec1.extend(list(kwargs['title']))
    rec1.append(header_size - 112)
    rec1_string = 'i %dc i' % (header_size - 112)
    ####
    ####
    # Third record
    rec2 = []
    rec2.append(4)
    rec2.append(int(kwargs['natoms']))
    rec2.append(4)
    rec2_string = 'i i i'
    #
    final_string = rec0_string + ' ' + rec1_string + ' ' + rec2_string
    final_rec = rec0 + rec1 + rec2
    return struct.pack(final_string, *final_rec)


def guess_natoms(fname, endian='@', rec_head_prec='i', xyz_prec='f'):
    """Attempts to guess the number of atoms in the file, by using
    fortran record marker information. This function should *only* be
    called on a charmm trajectory *with* record markers and *without* a
    trajectory header.

    Parameters
    ----------
    fname : a string representing the path to a fortran binary file
    endian: an optional character specifying the endian-ness of the file.
        Possible values are '>', '<', '@', '=' or '!'. Consult :mod:`struct`
        for more details. Defaults to '@', native byte order.
    rec_head_prec: an optional character specifying the fortran record header's
        precision. Possible values are 'h', 'i', 'l' or 'q'. Consult
        :mod:`struct` for more details. Defaults to 'i', single precision
        integers.
    xyz_prec: an optional character specifying the fortran record data's
        precision. Possible values are 'f' and 'd'. Consult :mod:`struct`
        for more details. Defaults to 'f', single precision floats.
    """
    with open_fort(fname, mode='r', buffering=None, endian=endian,
                rec_head_prec=rec_head_prec) as fp:
        tmp = fp.read_record() # ignore the first record, as it could contain xtl data
        l = fp._read_chk()
        if l % fp._byte_dict[xyz_prec]:
            raise ValueError("Invalid record length: %r" % l)
        tmp = fp._read_exactly(l)
        chk = fp._read_chk()
        if chk != l:
            raise IOError("Chk failure.")
        return l / fp._byte_dict[xyz_prec]


def guess_xtl(fname, endian='@', rec_head_prec='i', xtl_prec='d'):
    """Attempts to guess if the charmm trajectory has xtl information by
    examining the length of the first two fortran records. This function
    should *only* be called on a charmm trajectory *with* record markers
    and *without* a trajectory header.

    Parameters
    ----------
    fname : a string representing the path to a fortran binary file
    endian: an optional character specifying the endian-ness of the file.
        Possible values are '>', '<', '@', '=' or '!'. Consult :mod:`struct`
        for more details. Defaults to '@', native byte order.
    rec_head_prec: an optional character specifying the fortran record header's
        precision. Possible values are 'h', 'i', 'l' or 'q'. Consult
        :mod:`struct` for more details. Defaults to 'i', single precision
        integers.
    xtl_prec: an optional character specifying the fortran record data's
        precision. Possible values are 'f' and 'd'. Consult :mod:`struct`
        for more details. Defaults to 'd', single precision floats.
    """
    with open_fort(fname, mode='r', buffering=None, endian=endian,
                rec_head_prec=rec_head_prec) as fp:
        # read first frame, and see if it could be 6 real numbers
        l = fp._read_chk()
        if l % fp._byte_dict[xtl_prec]:
            raise ValueError("Invalid record length: %r" % l)
        tmp = fp._read_exactly(l)
        chk = fp._read_chk()
        if chk != l:
            raise IOError("Chk failure.")
        if l != 6 * fp._byte_dict[xtl_prec]: # could this be 6 reals?
            return False                            # no, then not xtl
        # first frame was 6 * 8 (if 'd') = 56 bytes long
        # this could be xtl data or fortuitous natoms
        # read second frame, if length is different, first frame was xtl
        l = fp._read_chk()
        if l % fp._byte_dict[xtl_prec]:
            raise ValueError("Invalid record length: %r" % l)
        tmp = fp._read_exactly(l)
        chk = fp._read_chk()
        if chk != l:
            raise IOError("Chk failure.")
        if l != 6 * fp._byte_dict[xtl_prec]: # could this be 6 reals?
            return True
        # first and second frame were 6 * 8 = 56 bytes long
        # this means we have no idea about xtl data
        raise ValueError("unable to guess has_xtl")


class DCDFile(CharmmBin):
    def __init__(self, fname, mode='rb', buffering=None,
                endian='@', rec_head_prec='i', **kwargs):
        super(DCDFile, self).__init__(fname=fname, mode=mode,
            buffering=buffering, endian=endian, rec_head_prec=rec_head_prec)
        """The `func:open_dcd` should be used in lieu of this constructor.
        """
        # assign kwargs
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
        # build frame data structure
        self._frame_dt = self.compile_npdt()
        self._frame_size = self._frame_dt.itemsize
        # etc
        self._leftovers = []
        if 'r' in mode:
            self.seek(0, whence=0)

    # Highest level methods ###################################################
    def iter_frame(self, begin=0, end=None, stride=None):
        if end is not None and not isinstance(end, int):
            raise TypeError("invalid type for end: %r" % end)
        if stride is not None and not isinstance(stride, int):
            raise TypeError("invalid type for stride: %r" % end)
        if stride is not None and stride < 1:
            raise ValueError("invalid value for stride: %r" % end)
        #
        self.seek(begin, whence=0)
        if end is None:
            if stride is None:
                while 1:
                    yield self.next()
            else:
                while 1:
                    yield self.next()
                    self.seek(stride-1, whence=1)
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
                    self.seek(stride-1, whence=1)
                    begin += stride

    def iter_nparray(self, begin=0, end=None, stride=None):
        if end is not None and not isinstance(end, int):
            raise TypeError("invalid type for end: %r" % end)
        if stride is not None and not isinstance(stride, int):
            raise TypeError("invalid type for stride: %r" % end)
        if stride is not None and stride < 1:
            raise ValueError("invalid value for stride: %r" % end)
        #
        self.seek(begin, whence=0)
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
                        self.seek(stride-1, whence=1)
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
                        self.seek(stride-1, whence=1)
                        begin += stride
                except MemoryError:
                    raise StopIteration

    def get_massive_dump(self):
        self.seek(0, whence=0)
        tmp = np.fromfile(self.fp, dtype=self.frame_dt, count=-1)
        self.readline()
        return tmp

    def read_frame(self):
        return self.readline()

    def read_nparray(self):
        return np.fromfile(self.fp, dtype=self.frame_dt, count=1)

    def put_massive_dump(self, dump, order='C'):
        # untested
        for frame in dump:
            self.write_nparray(frame, order)

    def update_header(self):
        tmp = pack_header(header_size=self.header_size,
                        dcdtype=self.dcdtype,
                        nframes=self.nframes,
                        npriv=self.npriv,
                        nsav=self.nsav,
                        nsteps=self.nsteps,
                        ndegfree=self.ndegfree,
                        has_xtl=int(self.has_xtl),
                        has_d4=int(self.has_d4),
                        validated=int(self.validated),
                        title=self.title,
                        natoms=self.natoms)
        self.fp.seek(0, 0)
        self.fp.write(tmp)

    def write_nparray(self, nparray, order='C'):
        if not isinstance(nparray, np.ndarray):
            raise TypeError("invalid nparray")
        self.fp.write(nparray.tostring(order))

    @property
    def leftovers(self):
        return self._leftovers

    # Frame (meta)data ########################################################
    def compile_npdt(self):
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
        if self.has_rec:
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
        tmp = self.fp.read(self.frame_size)
        if len(tmp) == self.frame_size:
            return tmp
        elif len(tmp) == 0:
            return b''
        else:
            warnings.warn("Did not read enough data, wanted %r bytes, received %r bytes." % (self.frame_size, len(tmp)))
            warnings.warn("Dumping fractional frame to leftovers.")
            self._leftovers.append(tmp)

    def seek(self, offset, whence=0):
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
        >>> taco.seek(0, 0)

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

    def tell(self):
        """Returns the current position of the filepointer as a function
        of the frame. If an integer is returned the filepointer is between
        frames, if a float is returned the filepointer is in a frame.
        """
        tmp = (self.fp.tell() - self.header_size) / self.frame_size
        if tmp == int(tmp):
            return int(tmp)
        else:
            return tmp

    def truncate(self, offset):
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
