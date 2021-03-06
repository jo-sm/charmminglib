"""
DOCME
"""
# fcp
# 02/22/2011


from os import remove
from os.path import exists
import numpy as np
from pychm.tools import Property, walk, lowerKeys, mkdir
from pychm.io.inp import INPFile


def load_correlOutput(filename):
    """
    Reads a charmm formatted correl output reporting a bondlength as a function
    of time.  Data is returned as a 1darray.
    """
    def gen():
        start_parse = False
        for line in open(filename):
            if start_parse:
                value = line.split()[1]
                if value.lower() == 'nan':
                    raise AssertionError("parsed NaN in %s" % filename)
                yield float(value)
            else:
                try:
                    dummy = int(line.split()[0])
                    value = line.split()[1]
                    if value.lower() == 'nan':
                        raise AssertionError('load_correlOutput: Parsed NaN in %s' % filename)
                    start_parse = True
                    yield float(value)
                except ValueError:
                    pass
    return np.fromiter(gen(), np.float)


class BaseAnalysis(INPFile):
    """
    DOCME
    """
    def __init__(self, pdbFilename=None, **kwargs):
        super(BaseAnalysis, self).__init__(pdbFilename, **kwargs)
        # kwargs
        kwargs = lowerKeys(kwargs)
        self.correlAtomSelection = kwargs.get('selection', 'all')
        #
        self._correlStart = 0
        self._correlStop = -1

# File/Path locations
    @Property
    def anlFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._anlFilename
        def fset(self, value):
            self._anlFilename = self.expandPath(value)
        return locals()

    @Property
    def anlPathname():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._anlPathname
        def fset(self, value):
            self._anlPathname = self.expandPath(value)
            mkdir(self._anlPathname)
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
            self._maxatom, self.nstep, self._timestep = self.query_dcd(value)
        return locals()


# Correl Properties
    @Property
    def correlArrayLength():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._correlArrayLength
        def fset(self, value):
            self._correlArrayLength = value
            try:
                self._correlSkip = self._nstep / value
            except AttributeError:
                pass
        return locals()

    @Property
    def correlSkip():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._correlSkip
        def fset(self, value):
            self._correlSkip = value
            try:
                self._correlArrayLength = self._nstep / value
            except AttributeError:
                pass
        return locals()

    @Property
    def correlStart():
        doc =\
        """
        dcd frame number to start correl processing, in quanta of `correlSkip`
        """
        def fget(self):
            if self._correlStart >= 0:
                return self.correlSkip * (self._correlStart + 1)
            else:
                return self.correlSkip * (self.correlArrayLength + self._correlStart + 1)
            return self._correlStart
        def fset(self, value):
            self._correlStart = int(value)
        return locals()

    @Property
    def correlStop():
        doc =\
        """
        dcd frame number to stop correl processing, in quanta of `correlSkip`
        """
        def fget(self):
            if self._correlStop >= 0:
                return self.correlSkip * (self._correlStop + 1)
            else:
                return self.correlSkip * (self.correlArrayLength + self._correlStop + 1)
            return self._correlStop
        def fset(self, value):
            self._correlStop = int(value)
        return locals()

    @Property
    def nstep():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._nstep
        def fset(self, value):
            self._nstep = value
            try:
                self._correlSkip = value / self._correlArrayLength
            except AttributeError:
                try:
                    self._correlArrayLength = value / self._correlSkip
                except AttributeError:
                    pass
        return locals()

##################
# Public Methods #
##################

    def get_correlInputHeader(self, header=[]):
        """
        DOCME
        """
        String = self.get_inputHeader(header)
        String.append('bomlev -2')
        String.append('')
        String.append('! dcd :: read')
        String.append('open unit 10 read unform name %s' % self.dcdFilename)
        String.append('')
        return String

    def rm_pickles(self, directory=None):
        """
        Recursively removes .pickle files starting at `directory`.  Defaults
        to `rootPath`.
        """
        if directory is None:
            directory = self.rootPath
        #
        iterator = ( path for path in walk(directory) if path.endswith('pickle') )
        for pickleFile in iterator:
            print 'Removing %s' % pickleFile
            remove(pickleFile)
