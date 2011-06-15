"""
:Author: fcp
:Date: 02/22/2011
"""


from commands import getstatusoutput
from os.path import abspath, dirname, expanduser, exists
from tempfile import NamedTemporaryFile
from charmming.tools import Property, expandPath, lowerKeys


class INPFile(object):
    """
    A base class for writing CHARMM .inp text files within charmminglib.


    This class provides functionality and book keeping for data common
    to nearly all CHARMM jobs, such as tracking .rtf, .prm, .psf and
    .crd file locations, writing .inp file headers and querying .dcd files.

    **Class Attributes:**
        | ``defaultBin``

    **kwargs:**
        | ``charmmbin`` :: Path to binary used to execute CHARMM,
            defaults to value in ``defaultBin``
    """

    defaultBin = 'cg_mscale'
    """
    If not explicitly specified as a kwarg, this is path used by the class
    to run CHARMM jobs through the shell.
    """

    def __init__(self, pdbFilename=None, **kwargs):
        super(INPFile, self).__init__()
        # kwargs
        kwargs = lowerKeys(kwargs)
        self.charmmBin = kwargs.get('charmmbin', self.__class__.defaultBin)
        # Main
        if pdbFilename is not None:
            self.pdbFilename = pdbFilename
        #
        self._rtfFilename = None
        self._prmFilename = None
        self._psfFilename = None
        self._crdFilename = None

# charmminglib file locations
    @Property
    def inpFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .inp file to be written.
        """
        def fget(self):
            return self._inpFilename
        def fset(self, value):
            self._inpFilename = self.expandPath(value)
        return locals()

    @Property
    def outFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .out file to be written.
        """
        def fget(self):
            return self._outFilename
        def fset(self, value):
            self._outFilename = self.expandPath(value)
        return locals()

# CHARMM File/Path locations
    @Property
    def crdFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .crd file to be read from.
        """
        def fget(self):
            return self._crdFilename
        def fset(self, value):
            value = self.expandPath(value)
            if not exists(value):
                raise IOError("No such file or directory: '%s'" % value)
            self._crdFilename = value
        return locals()

    @Property
    def dcdFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .dcd file to be read from.
        """
        def fget(self):
            return self._dcdFilename
        def fset(self, value):
            value = self.expandPath(value)
            if not exists(value):
                raise IOError("No such file or directory: '%s'" % value)
            self._dcdFilename = value
            self._maxatom, self._nstep, self._timestep = self.query_dcd(value)
        return locals()

    @Property
    def pdbFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .pdb file to be read from.
        """
        def fget(self):
            return self._pdbFilename
        def fset(self, value):
            value = expandPath(value)
            if not exists(value):
                raise IOError("No such file or directory: '%s'" % value)
            self._pdbFilename = value
        return locals()

    @Property
    def prmFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .prm file to be read from.
        """
        def fget(self):
            return self._prmFilename
        def fset(self, value):
            value = self.expandPath(value)
            if not exists(value):
                raise IOError("No such file or directory: '%s'" % value)
            self._prmFilename = value
        return locals()

    @Property
    def psfFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .psf file to be read from.
        """
        def fget(self):
            return self._psfFilename
        def fset(self, value):
            value = self.expandPath(value)
            if not exists(value):
                raise IOError("No such file or directory: '%s'" % value)
            self._psfFilename = value
        return locals()

    @Property
    def rootPath():
        doc =\
        """
        A read-only ``property`` for the directory where ``self.pdbFilename``
        resides.  This allows for easy hierarchical organization in
        conjunction with the ``self.expandPath`` method.
        """
        def fget(self):
            return dirname(self.pdbFilename)
        return locals()

    @Property
    def rtfFilename():
        doc =\
        """
        A ``property`` for the name of the CHARMM .rtf file to be read from.
        """
        def fget(self):
            return self._rtfFilename
        def fset(self, value):
            value = self.expandPath(value)
            if not exists(value):
                raise IOError("No such file or directory: '%s'" % value)
            self._rtfFilename = value
        return locals()

# DCD Properties
    @Property
    def maxatom():
        doc =\
        """
        If a valid .dcd file and valid charmmBin are specified, this
        ``property`` will return the number of atoms present in the
        trajectory.
        """
        def fget(self):
            return self._maxatom
        return locals()

    @Property
    def nstep():
        doc =\
        """
        If a valid .dcd file and valid charmmBin are specified, this
        ``property`` will return the number of frames present in the
        trajectory.
        """
        def fget(self):
            return self._nstep
        return locals()

    @Property
    def timestep():
        doc =\
        """
        If a valid .dcd file and valid charmmBin are specified, this
        ``property`` will return the timestep used in the trajectory.
        """
        def fget(self):
            return self._timestep
        return locals()

##################
# Public Methods #
##################

    def expandPath(self, arg):
        """
        A better version of the default path expander.  In addition to
        expanding ``~`` to the current user's home directory, it will
        also expand ``&`` to the current ``rootPath``.
        """
        arg = arg.strip()
        if arg.startswith('&'):
            try:
                return abspath(arg.replace('&', self.rootPath))
            except AttributeError:
                pass
        elif '~' in arg:
            return expanduser(arg)
        else:
            return abspath(arg)

    def get_inputHeader(self, header=[]):
        """
        Returns a list of strings which correspond to the typical
        "open unit, read unit" lines for .rtf, .prm, .psf and .crd
        files.

        ``header`` is an optional list which is prepended to the list
        returned by this method.
        """
        String = []
        for line in header:
            String.append('* %s' % line)
        String.append('*')
        String.append('bomlev -1')
        String.append('wrnlev 5')
        String.append('')
        if self.rtfFilename:
            String.append('! rtf')
            String.append('  open unit 1 read card name %s' % self.rtfFilename)
            String.append('  read rtf card unit 1')
            String.append('  close unit 1')
            String.append('')
        if self.prmFilename:
            String.append('! prm')
            String.append('  open unit 1 read card name %s' % self.prmFilename)
            String.append('  read para card unit 1')
            String.append('  close unit 1')
            String.append('')
        if self.psfFilename:
            String.append('! psf')
            String.append('  open unit 1 read card name %s' % self.psfFilename)
            String.append('  read psf card unit 1')
            String.append('  close unit 1')
            String.append('')
        if self.crdFilename:
            String.append('! crd')
            String.append('  open unit 1 read card name %s' % self.crdFilename)
            String.append('  read coor card unit 1')
            String.append('  close unit 1')
            String.append('')
        return String

    def query_dcd(self, dcdFilename=None):
        """
        Queries a .dcd file using the specified ``charmmBin`` and
        returns a 3-tuple of integers corresponding to the number of
        atoms, the number of steps, and the timestep in fs of the
        specified .dcd file.  Defaults to ``self.dcdFilename``.
        """
        if dcdFilename is None:
            dcdFilename = self.dcdFilename
        else:
            dcdFilename = self.expandPath(dcdFilename)
        # tmp file
        tmp = NamedTemporaryFile()
        tmpOut = []
        tmpOut.append('* trj query')
        tmpOut.append('*')
        tmpOut.append('open unit 10 read unform name %s' % dcdFilename)
        tmpOut.append('traj query unit 10')
        tmpOut.append('stop')
        tmpOut = '\n'.join(tmpOut)
        tmp.file.write(tmpOut)
        tmp.file.flush()
        # execute charmm
        charmmOut = getstatusoutput('%s < %s' % (self.charmmBin, tmp.name))[1].split('\n')
        # parse output
        nCrd = None
        freq = None
        nDeg = None
        nFix = None
        ts = None
        iterator = ( line.strip().lower() for line in charmmOut )
        for line in iterator:
            if line.startswith('number of coordinate'):
                nCrd = int(line.split()[-1])
            elif line.startswith('frequency'):
                freq = int(line.split()[-1])
            elif line.startswith('number of degrees of freedom'):
                nDeg = int(line.split()[-1])
            elif line.startswith('number of fixed atoms'):
                nFix = int(line.split()[-1])
            elif line.startswith('the integration time'):
                ts = float(line.split()[-1])
        maxatom = int(nDeg / 3 + nFix) + 1
        nstep = int(nCrd * freq)
        timestep = int(ts * 1000)
        return (maxatom, nstep, timestep)
