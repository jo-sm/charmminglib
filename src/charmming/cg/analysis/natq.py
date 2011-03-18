

import os
import numpy as np
from copy import deepcopy
from charmming.tools import Property, mkdir, lowerKeys
from charmming.analysis.baseanalysis import BaseAnalysis,load_correlOutput
from charmming.io.pdb import PDBFile
from charmming.cg.ktgo import KTGo


class NatQ(BaseAnalysis):
    """
    DOCME
    """
    def __init__(self, arg=None, **kwargs):
        super(NatQ, self).__init__(arg, **kwargs)
        self.inpFilename = '&/natq.inp'
        self.outFilename = '&/natq.out'
        self._anlPathname = None
        # kwargs
        self._nativeRad = kwargs.get('nativerad', 8.)
        # Main
        if arg is not None:
            self.aa = PDBFile(self.pdbFilename)[0]
            self.aa.parse()
            self.cg = KTGo(self.aa)
            self.natq = self.cg.get_nativeSCSC()

    @Property
    def nativeRad():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._nativeRad
        def fset(self, value):
            if value <= 1:
                raise ValueError("The native contact radius must be at least 1.")
            self._nativeRad = float(value)
        return locals()

    def calc_timeseries(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        filterFunction = kwargs.get('filterfunction', lambda x: True)
        name = kwargs.get('name', 'all')
        #
        contacts = []
        for i, contact in enumerate(self.cg.get_nativeSCSC()):
            if filterFunction(contact):
                contacts.append(i)
        iterator = ( '%s/nat%03d.anl' % (self.anlPathname, i) for i in contacts )
        matrix = np.array([ load_correlOutput(file) for file in iterator ])
        matrix[matrix <= self.nativeRad] = 1
        matrix[matrix > self.nativeRad] = 0
        return matrix.mean(axis=0)

    def do_correl(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        #
        self.write_correlInput(**kwargs)
        self.run_correlInput(**kwargs)

    def write_correlInput(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        #
        if self.anlPathname is None:
            raise IOError("No directory specified.")
        mkdir(self.anlPathname)
        #
        contacts = map(lambda (x, y): (x.atomNum, y.atomNum), self.natq)
        #
        String = self.get_correlInputHeader(header)
        String.append('! anl :: write')
        for i, contact in enumerate(contacts):
            String.append('open unit %03d write card name %s/nat%03d.anl' % (i+100, self.anlPathname, i))
        String.append('')
        String.append('correl maxtimesteps %d maxatom %d maxseries %d' % (self.correlArrayLength, 2 * self.maxatom, len(self.natq)))
        for i, contact in enumerate(contacts):
            String.append('enter n%03d bond bynum %d bynum %d geometry' % (i, contact[0], contact[1]))
        String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select all end' % (self.correlStart, self.correlStop, self.correlSkip))
        String.append('')
        for i, contact in enumerate(contacts):
            String.append('write n%03d card unit %03d' % (i, i+100))
            String.append('* Native contact %d: between cgAtoms %d and %d' % (i, contact[0], contact[1]))
            String.append('*')
            String.append('')
        String.append('stop')
        String = '\n'.join(String)
        # Write file
        write_to = open('%s' % self.inpFilename, 'w')
        write_to.write(String)
        write_to.close()

    def run_correlInput(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        #
        try:
            os.remove(self.outFilename)
        except OSError:
            pass
        if not os.path.exists(self.inpFilename):
            raise IOError("No such file or directory: '%s'" % self.inpFilename)
        os.system('%s < %s > %s' % (self.charmmBin, self.inpFilename, self.outFilename))
