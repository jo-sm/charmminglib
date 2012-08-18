"""
:Author: btm
:Date: 08/17/2012
"""

class MOL2File(object):

    def _partition(self):
        inHeader = False
        inCrd = False
        inBonds = False

        iterator = ( line for line in cleanStrings(open(self.filename) ))

        self._header = []
        self._crd = []
        self._bonds = []

        for line in iterator:

            line = line.strip()
            if not line:
                continue
            if line.startswith('@<TRIPOS>'):
                sectionName = line[9:].strip()
                if sectionName == 'MOLECULE'
                    inHeader = True
                    inCrd = False
                    inBonds = False

                    if self._header:
                        raise('Multiple header sections')
                elif sectionName == 'ATOM':
                    inHeader = False
                    inCrd = True
                    inBonds = False
                elif sectionName == 'BOND':
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
                    self._bonds.append(line)

    def _buildmodel(self):
        iterator = ( Atom(text=line, informat='mol2', index=i) ) \
                   for i, line in enumerate(self._crd)
        self._mymol = Mol(iterable=iterator, name='model0', autofix=True)

    def write_out(self, filename, **kwargs):
        pass

    @Property
    def mol(self):
        return self._mymol

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
