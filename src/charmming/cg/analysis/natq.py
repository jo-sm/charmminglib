

import os
import numpy as np
from copy import deepcopy
from charmming.tools import Property, mkdir, lowerKeys, grouper
from charmming.lib.bond import Bond
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
            self.contacts = self.get_contacts()

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

    def calc_avg_contactMatrix(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        t0 = kwargs.get('starttime', None)
        tf = kwargs.get('stoptime', None)
        # Build Null Matrix
        poly = np.lib.polynomial.poly1d([0.5, -0.5, -1 * len(self.contacts)])
        dim = -1
        for i in poly.r:
            if i > 0:
                dim = i
        matrix = np.zeros((dim, dim))
        # Build atomNum -> index Map
        indexMap = dict(( (atom.atomNum, i) for i, atom in enumerate(self.cg.find(atomtype='s')) ))
        #
        for n, contact in enumerate(self.contacts):
            filename = '%s/nat%04d.anl' % (self.anlPathname, n)
            tmp = load_correlOutput(filename)
            tmp = tmp[t0:tf]
            tmp[tmp <= self.nativeRad] = 1
            tmp[tmp > self.nativeRad] = 0
            i = indexMap[contact.i.atomNum]
            j = indexMap[contact.j.atomNum]
            matrix[i][j] = matrix[j][i] = tmp.mean()
        return matrix

    def calc_natq_of_t(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        filterFunction = kwargs.get('filterfunction', lambda x: True)
        name = kwargs.get('name', 'all')
        #
        index = []
        iterNative = ( (i, contact) for i, contact in enumerate(self.contacts) if contact.native )
        for i, contact in iterNative:
            if filterFunction(contact):
                index.append(i)
        iterator = ( '%s/nat%04d.anl' % (self.anlPathname, i) for i in index )
        matrix = np.array([ load_correlOutput(filename) for filename in iterator ])
        matrix[matrix <= self.nativeRad] = 1
        matrix[matrix > self.nativeRad] = 0
        return matrix.mean(axis=0)

    def do_correl(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        #
        self.write_correlInput(**kwargs)
        self.run_correlInput(**kwargs)

    def get_contacts(self):
        iterator = [ atom for atom in self.cg if atom.atomType == 's' ]
        tmp = [ Bond(atom_i, atom_j) for i, atom_i in enumerate(iterator) for
                    j, atom_j in enumerate(iterator) if i < j ]
        nativeContacts = set(self.cg.get_nativeSCSC())
        for contact in tmp:
            if contact in nativeContacts:
                contact.native = True
            else:
                contact.native = False
        return tmp

    def write_correlInput(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        #
        if self.anlPathname is None:
            raise IOError("No directory specified.")
        mkdir(self.anlPathname)
        #
        String = self.get_correlInputHeader(header)
        for i, group in enumerate(grouper(self.contacts, 100)):
            String.append('! anl :: write')
            for j, contact in enumerate(group):
                String.append('open unit %03d write card name %s/nat%02d%02d.anl' % (j+100, self.anlPathname, i, j))
            String.append('')
            String.append('correl maxtimesteps %d maxatom %d maxseries %d' % (self.correlArrayLength, 1000, len(group)))
            for j, contact in enumerate(group):
                String.append('enter n%03d bond bynum %d bynum %d geometry' % (j, contact.i.atomNum, contact.j.atomNum))
            String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select all end' % (self.correlStart, self.correlStop, self.correlSkip))
            String.append('')
            for j, contact in enumerate(group):
                String.append('write n%03d card unit %03d' % (j, j+100))
                String.append('* Contact %02d%02d: between cgAtoms %s and %s' % (i, j, contact.i.addr, contact.j.addr))
                if contact.native:
                    String.append('* Native Contact - Interaction between %s and %s' % (contact.i.prmString, contact.j.prmString))
                String.append('*')
                String.append('')
            String.append('end')
            String.append('')
            String.append('rewind unit 10')
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

    def bak_write_correlInput(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        #
        if self.anlPathname is None:
            raise IOError("No directory specified.")
        mkdir(self.anlPathname)
        #
        native = self.cg.get_nativeSCSC()
        #
        String = self.get_correlInputHeader(header)
        String.append('! anl :: write')
        for i, contact in enumerate(native):
            String.append('open unit %03d write card name %s/nat%d.anl' % (i+100, self.anlPathname, i))
        String.append('')
        String.append('correl maxtimesteps %d maxatom %d maxseries %d' % (self.correlArrayLength, 40 * self.maxatom, len(native)))
        for i, contact in enumerate(native):
            String.append('enter n%03d bond bynum %d bynum %d geometry' % (i, native.i.atomNum, native.j.atomNum))
        String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select all end' % (self.correlStart, self.correlStop, self.correlSkip))
        String.append('')
        for i, contact in enumerate(native):
            String.append('write n%03d card unit %03d' % (i, i+100))
            String.append('* Native contact %d: between cgAtoms %s and %s' % (i, native.i.addr, native.j.addr))
            String.append('*')
            String.append('')
        String.append('stop')
        String = '\n'.join(String)
        # Write file
        write_to = open('%s' % self.inpFilename, 'w')
        write_to.write(String)
        write_to.close()
