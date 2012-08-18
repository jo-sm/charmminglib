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

    def __init__(self, filename=None, **kwargs):
        super(MOL2File, self).__init__()

        if filename is not None:
            self.filename = filename
            self._partition()
        else:            
            self.filename = 'null'
            self._header = 'null'
            self._crd = 'null'
            self._bonds = 'null'
