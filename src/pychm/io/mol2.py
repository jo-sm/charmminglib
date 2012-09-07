"""
Author: btm
:Date: 08/17/2012
"""

from pychm.tools import Property, cleanStrings
from pychm.lib.bond import Bond
from pychm.lib.mol import Atom, Mol

class MOL2Bond(Bond):

    def writeOut(self):
        return '%6s%6s%6s%5s' % (self.seq,self.i.atomNum,self.j.atomNum,self.order)

class MOL2File(object):

    def _partition(self):
        inHeader = False
        inCrd = False
        inBonds = False

        iterator = ( line for line in cleanStrings(open(self.filename) ))

        self._header = []
        self._crd = []
        self._bondstext = []

        for line in iterator:

            if not line:
                continue

            if line.startswith('@<tripos>'):
                sectionName = line[9:].strip()
                if sectionName == 'molecule':
                    inHeader = True
                    inCrd = False
                    inBonds = False

                    if self._header:
                        raise('Multiple header sections')
                elif sectionName == 'atom':
                    inHeader = False
                    inCrd = True
                    inBonds = False
                elif sectionName == 'bond':
                    if not self._crd:
                        raise('')

                    inHeader = False
                    inCrd = False
                    inBonds = True
                else:
                    # Throw away other sections
                    inHeader = False
                    inCrd = False
                    inBonds = False
            else:
                if inHeader:
                    self._header.append(line)
                elif inCrd:
                    self._crd.append(line)
                elif inBonds:
                    self._bondstext.append(line)

    def _buildbonds(self):
        self._bonds = []

        for bondline in self._bondstext:
            seq, atmi, atmj, order = bondline.split()
            tmp_bond = MOL2Bond(self._mymol[int(atmi)-1], self._mymol[int(atmj)-1])
            tmp_bond.seq = seq
            tmp_bond.order = order
            self._bonds.append(tmp_bond)

    def _buildmodel(self):
        iterator = ( Atom(text=line, informat='mol2', index=i) for i, line in enumerate(self._crd) )
        self._mymol = Mol(iterable=iterator, name='model0', autofix=True)
        self._buildbonds()

    def write_out(self, filename, **kwargs):
        pass

    @Property
    def mol():
        def fget(self):
            return self._mymol
        return locals()

    @Property
    def header():
        def fget(self):
            return self._header
        return locals()

    @Property
    def bonds():
        def fget(self):
            return self._bonds
        return locals()

    def __init__(self, filename=None, **kwargs):
        super(MOL2File, self).__init__()

        if filename is not None:
            self.filename = filename
            self._partition()
            self._buildmodel()
        else:            
            self.filename = 'null'
            self._header = 'null'
            self._crd = 'null'
            self._bonds = 'null'
