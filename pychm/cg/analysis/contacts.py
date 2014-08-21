

import os
import numpy as np
import cPickle as pickle
from copy import deepcopy
from pychm.tools import Property, mkdir, lowerKeys, grouper
from pychm.lib.bond import Bond
from pychm.analysis.baseanalysis import BaseAnalysis, load_correlOutput
from pychm.io.pdb import PDBFile
from pychm.cg.ktgo import KTGo


class Contacts(BaseAnalysis):
    """
    DOCME
    """
    def __init__(self, pdbFilename=None, **kwargs):
        super(Contacts, self).__init__(pdbFilename, **kwargs)
        # kwargs
        self.nativeRad = kwargs.get('nativerad', 8.)
        # Main
        if pdbFilename is not None:
            self.aa = PDBFile(self.pdbFilename)[0]
            self.aa.parse()
            self.cg = KTGo(self.aa)
            self.contacts = self.get_contacts()

    @Property
    def data():
        doc =\
        """
        """
        def fget(self):
            filename = '%s/contacts.pickle' % self.anlPathname
            try:
                self._data = pickle.load(open(filename))
                print 'found pickled data in: %s' % filename
                return self._data
            except IOError:
                print 'processing data in: %s' % self.anlPathname
                self._data = self.build_contactMatrix()
                pickle.dump(self._data, open(filename, 'w'))
                return self._data
        def fset(self, value):
            assert isinstance(value, np.ndarray)
            self._data = value
        def fdel(self):
            filename = '%s/contacts.pickle' % self.anlPathname
            del self._data
            try:
                os.remove(filename)
            except IOError:
                pass
        return locals()

    @Property
    def inpPathname():
        doc =\
        """
        """
        def fget(self):
            return self._inpPathname
        def fset(self, value):
            self._inpPathname = self.expandPath(value)
        return locals()

    @Property
    def nativeRad():
        doc =\
        """
        The cutoff distance between two CG sidechain atom objects which is used
        to determine if their parent residues are in 'contact'.
        """
        def fget(self):
            return self._nativeRad
        def fset(self, value):
            if value <= 1:
                raise ValueError("The native contact radius must be at least 1.")
            self._nativeRad = float(value)
        return locals()

    def get_contacts(self):
        """
        Returns a complete :class:`list` of :class:`Bond` objects, each bond
        object represents a unique sidechain-sidechain pair. Additionally,
        and most importantly, each :class:`Bond` object has been assigned the
        additional ``native`` attribute, to ``True`` or ``False``.  This value
        is true, if the corresponding all-atom representation of the residue
        sidechains have any heavy atoms within 4.5 Angstrom, as determined by
        :meth:`KTGo.get_nativeSCSC`.
        """
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

    def build_contactMatrix(self):
        """
        Builds a 2D numpy array from available correl output files. The first
        index references the contact number, the second index references the
        correl time array value (which in turn corresponds to nstep of the
        original dynamics in quanta of ``correlSkip``. This method can take
        awhile because of disk IO.
        """
        print "Building full contact matrix."
        rowLength = (self.correlStop - self.correlStart) / self.correlSkip + 1
        tmp = []
        for i in xrange(len(self.contacts)):
            try:
                tmp.append(load_correlOutput('%s/contact%04d.anl' % (self.anlPathname, i)))
            except IOError:
                tmp.append(np.zeros(rowLength))
                print "Can't find correl output for contact number: %04d." % i
        return np.array(tmp, dtype=np.float64)

    def get_nativeContactMatrix(self):
        """
        Projects out elements of the contactMatrix which are not derived from
        'native' contacts.
        """
        nativeIndex = [ i for i, taco in enumerate(self.contacts) if taco.native ]
        return self.data[nativeIndex]

    def get_natQofT(self):
        """
        """
        nativeContacts = self.get_nativeContactMatrix()
        tmp = np.zeros(nativeContacts.shape)
        tmp[nativeContacts <= self.nativeRad] = 1
        tmp[nativeContacts > self.nativeRad] = 0
        return tmp.mean(axis=0)

# Charmm input creation
    def do_correl(self, **kwargs):
        self.write_correlInput(**kwargs)
        self.run_correlInput()

    def write_correlInput(self, **kwargs):
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        #
        if self.anlPathname is None:
            raise IOError("No directory specified.")
        mkdir(self.anlPathname)
        #
        for i, group in enumerate(grouper(self.contacts, 100)):
            String = self.get_correlInputHeader(header)
            String.append('! anl :: write')
            for j, contact in enumerate(group):
                String.append('open unit %03d write card name %s/contact%02d%02d.anl' % (j+100, self.anlPathname, i, j))
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
            self.inpFilename = '%s/contact%02d.inp' % (self.inpPathname, i)
            write_to = open(self.inpFilename, 'w')
            write_to.write(String)
            write_to.close()

    def run_correlInput(self):
        #
        for i, group in enumerate(grouper(self.contacts, 100)):
            self.inpFilename = '%s/contact%02d.inp' % (self.inpPathname, i)
            self.outFilename = '%s/contact%02d.out' % (self.outPathname, i)
            #
            try:
                os.remove(self.outFilename)
            except OSError:
                pass
            if not os.path.exists(self.inpFilename):
                raise IOError("No such file or directory: '%s'" % self.inpFilename)
            os.system('%s < %s > %s' % (self.charmmBin, self.inpFilename, self.outFilename))
