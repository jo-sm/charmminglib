"""
DOCME
"""
# fcp
# 02/22/2011


from commands import getstatusoutput
from os.path import abspath, dirname, expanduser
from tempfile import NamedTemporaryFile
from charmming.tools import Property, expandPath


class BaseCHARMMFile(object):
    """
    DOCME
    """

    charmmBin = 'cg_mscale'

    def __init__(self, arg=None):
        super(BaseCHARMMFile, self).__init__()
        if arg is not None:
            self.pdbFilename = arg
        #
        self._rtfFilename = None
        self._prmFilename = None
        self._psfFilename = None
        self._crdFilename = None

# File/Path locations
    @Property
    def crdFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._crdFilename
        def fset(self, value):
            self._crdFilename = self.expandPath(value)
        return locals()

    @Property
    def dcdFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._dcdFilename
        def fset(self, value):
            value = self.expandPath(value)
            self._dcdFilename = value
            self._maxatom, self._nstep, self._timestep = self.query_dcd(value)
        return locals()

    @Property
    def pdbFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._pdbFilename
        def fset(self, value):
            self._pdbFilename = expandPath(value)
        return locals()

    @Property
    def prmFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._prmFilename
        def fset(self, value):
            self._prmFilename = self.expandPath(value)
        return locals()

    @Property
    def psfFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._psfFilename
        def fset(self, value):
            self._psfFilename = self.expandPath(value)
        return locals()

    @Property
    def rootPath():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return dirname(self.pdbFilename)
        return locals()

    @Property
    def rtfFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._rtfFilename
        def fset(self, value):
            self._rtfFilename = self.expandPath(value)
        return locals()

# DCD Properties
    @Property
    def maxatom():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._maxatom
        return locals()

    @Property
    def nstep():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._nstep
        return locals()

    @Property
    def timestep():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._timestep
        return locals()

##################
# Public Methods #
##################

    def expandPath(self, arg):
        """
        A better version of the default path expander.
        """
        if '&' in arg:
            try:
                return arg.replace('&', self.rootPath)
            except AttributeError:
                pass
        elif '~' in arg:
            return expanduser(arg)
        else:
            return abspath(arg)

    def get_inputHeader(self, header=[]):
        """
        DOCME
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
            String.append(' read rtf card unit 1')
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

    def query_dcd(self, dcdFilename):
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
        maxatom = int(nDeg / 3 + nFix)
        nstep = int(nCrd * freq)
        timestep = int(ts * 1000)
        return (maxatom, nstep, timestep)
