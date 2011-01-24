

import os
import numpy as np
from charmming.tools import Property,get_inpProp,walk,mkdir,out2inp,logicalLines


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
            except ValueError: pass
    return np.fromiter(gen(),np.float)


class BaseAnalysis(object):
    """docstring for NativeQ"""
    def __init__(self,aaCrdFile):
        """
        """
        super(BaseAnalysis,self).__init__()
        #
        self.charmmBin  = 'cg_charmm'
        # File/Path Defaults
        self.aaCrdFile  = aaCrdFile

        # TODO
        #self.atomSelect         = 'all' # Select by chainid, eg. 'ABD' or 'all'

        # Correl Defaults
        self.correlAtomSelect   = 'all' # Select by chainid, eg. 'ABD' or 'all'
        self.correlMaxAtom      = 300
        self.correlArrayLength  = 10000
        self.correlSkip         = 'auto'
        self.correlStart        = 0     # Pythonic indexing of "_correlArray"
        self.correlStop         = -1    # ie, [0:-1] will get you everything
                                        # in quanta of _correlSkip
        self.nstep              = 'auto'
        self.timestep           = 'auto'

# File/Path locations
    @Property
    def aaCrdFile():
        doc = "The aaCrdFile property."
        def fget(self):
            return self._aaCrdFile
        def fset(self, value):
            self._aaCrdFile = self.expandPath(value)
        return locals()

    @Property
    def aaCrdPath():
        doc = "The aaCrdPath property."
        def fget(self):
            return os.path.dirname(self.aaCrdFile)
        return locals()

    @Property
    def anlPath():
        doc = "The anlPath property."
        def fget(self):
            return self._anlPath
        def fset(self, value):
            value = self.expandPath(value)
            mkdir(value)
            if os.path.isdir:
                self._anlPath = value
            else:
                print 'anlPath: failed to create %s' % value
        return locals()

    @Property
    def rtfFile():
        doc = "A CHARMM .rtf file."
        def fget(self):
            return self._rtfFile
        def fset(self, value):
            self._rtfFile = self.expandPath(value)
        return locals()

    @Property
    def prmFile():
        doc = "A CHARMM .prm file."
        def fget(self):
            return self._prmFile
        def fset(self, value):
            self._prmFile = self.expandPath(value)
        return locals()

    @Property
    def psfFile():
        doc = "A CHARMM .psf file."
        def fget(self):
            return self._psfFile
        def fset(self, value):
            self._psfFile = self.expandPath(value)
        return locals()

    @Property
    def crdFile():
        doc = "A CHARMM .crd  file."
        def fget(self):
            return self._crdFile
        def fset(self, value):
            self._crdFile = self.expandPath(value)
        return locals()

    @Property
    def dcdFile():
        doc = "A CHARMM .dcd file."
        def fget(self):
            return self._dcdFile
        def fset(self, value):
            self._dcdFile = self.expandPath(value)
        return locals()

    @Property
    def dcdPath():
        doc = "A path containing all of the CHARMM .dcd files."
        def fget(self):
            return self._dcdPath
        def fset(self, value):
            self._dcdPath = self.expandPath(value)
        return locals()

    @Property
    def inpFile():
        doc = "A CHARMM .inp file."
        def fget(self):
            return self._inpFile
        def fset(self, value):
            self._inpFile = self.expandPath(value)
        return locals()

    @Property
    def outFile():
        doc = "A CHARMM .out file."
        def fget(self):
            return self._outFile
        def fset(self, value):
            self._outFile = self.expandPath(value)
        return locals()

    @Property
    def outPath():
        doc = "The path containing all of the CHARMM .out files."
        def fget(self):
            return self._outPath
        def fset(self, value):
            self._outPath = self.expandPath(value)
        return locals()

# DCD Properties
    @Property
    def nstep():
        doc = "Number of frames in dcd file, defaults to 'auto'"
        def fget(self):
            if self._nstep == 'auto':
                inp = open(self.inpFile)
                inp = logicalLines(inp)
                return int(get_inpProp('nste',inp))
            else:
                return int(self._nstep)
        def fset(self, value):
                self._nstep = value
        return locals()

    @Property
    def timestep():
        doc = "Trajectory timestep in picoseconds, defaults to 'auto'"
        def fget(self):
            if self._timestep == 'auto':
                inp = open(self.inpFile)
                inp = logicalLines(inp)
                return float(get_inpProp('time',inp))
            else:
                return float(self._nstep)
            return self._timestep
        def fset(self, value):
            self._timestep = value
        return locals()

# Correl Properties
    @Property
    def correlArrayLength():
        doc = "Number of correl writes per dcd file, defaults to 10,000"
        def fget(self):
            return self._correlArrayLength
        def fset(self, value):
            self._correlArrayLength = int(value)
        return locals()

    @Property
    def correlMaxAtom():
        doc = "Maximum number of atoms present in the dcd file."
        def fget(self):
            return self._correlMaxAtom
        def fset(self, value):
            self._correlMaxAtom = int(value)
        return locals()

    @Property
    def correlSkip():
        doc = "Number of dcd frames to skip between correl writes, defaults to 'auto'"
        def fget(self):
            if self._correlSkip == 'auto':
                return int(self.nstep / self.correlArrayLength)
            else:
                return self._correlSkip
        def fset(self, value):
            self._correlSkip = value
        return locals()

    @Property
    def correlStart():
        doc = "dcd frame number to start correl processing, in quanta of correlSkip"
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
        doc = "dcd frame number to stop correl processing, in quanta of correlSkip"
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
    def correlAtomSelect():
        doc = "The correlAtomSelect property."
        def fget(self):
            return self._correlAtomSelect
        def fset(self, value):
            self._correlAtomSelect = value
        return locals()

# Analysis Properties
    #@Property
    #def atomSelect():
    #    doc = "The atomSelect property."
    #    def fget(self):
    #        return self._atomSelect
    #    def fset(self, value):
    #        self._atomSelect = 'all'
    #    return locals()

# Public Methods
    def rm_pickles(self,path):
        """
        Recursively searches and removes *.pickle from the specified path.
        """
        # Purgation
        pickleFiles = ( fileName for fileName in walk(path) if fileName.endswith('pickle') )
        for pickleFile in pickleFiles:
            print 'removing %s' % pickleFile
            os.remove(pickleFile)

    def get_correlInputHeader(self, header=[]):
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
        String.append('! toppar')
        String.append('read rtf  card name %s' % self.rtfFile)
        String.append('read para card name %s' % self.prmFile)
        String.append('')
        String.append('! psfcor')
        String.append('read psf  card name %s' % self.psfFile)
        String.append('read coor card name %s' % self.crdFile)
        String.append('')
        String.append('bomlev -2')
        String.append('')
        String.append('! open trajectory for reading')
        String.append('open unit 10 read unform name %s' % self.dcdFile)
        String.append('')
        return '\n'.join(String)

    def expandPath(self,string):
        """
        A better version of the default path expander.
        """
        if '&' in string:
            try:
                return string.replace('&',self.aaCrdPath)
            except AttributeError: pass
        elif '~' in string:
            return os.path.expanduser(string)
        else:
            return os.path.abspath(string)
