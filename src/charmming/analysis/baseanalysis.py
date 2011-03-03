"""
DOCME
"""
# fcp
# 02/22/2011


import numpy as np
from charmming.tools import Property
from charmming.io.basecharmm import BaseCHARMMFile


def load_correlOutput(filename):
    """
    Reads a charmm formatted correl output reporting a bondlength as a function
    of time.  Data is returned as a 1darray.
    """
    def gen():
        for line in open(filename):
            try:
                dummy = int(line.split()[0])
                value = line.split()[1]
                if value.lower() == 'nan':
                    raise ValueError('load_correlOutput: Parsed NaN in %s' % filename)
                yield float(value)
            except ValueError:
                pass
    return np.fromiter(gen(), np.float)


class BaseAnalysis(BaseCHARMMFile):
    """
    DOCME
    """
    def __init__(self, arg=None):
        super(BaseAnalysis, self).__init__(arg)
        #
        self._correlStart = 0
        self._correlStop = -1

# File/Path locations
    @Property
    def anlPath():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._anlPath
        def fset(self, value):
            self._anlPath = self.expandPath(value)
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
            self._correlSkip = self._nstep / value
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
            self._correlArrayLength = self._nstep / value
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

##################
# Public Methods #
##################

    def get_correlInputHeader(self, header=[]):
        """
        DOCME
        """
        String = super(BaseCHARMMFile, self).get_inputHeader(header)
        String.append('bomlev -2')
        String.append('')
        String.append('! dcd :: read')
        String.append('open unit 10 read unform name %s' % self.dcdFilename)
        String.append('')
        return String
