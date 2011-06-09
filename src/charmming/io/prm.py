"""
:Author: fcp
:Date: 02/22/2011
"""


from charmming.tools import paragraphs
from charmming.lib.toppar import AnglePRM, BondPRM, DihedralPRM, ImproperPRM, \
        MassPRM, NonBondPRM, NBFixPRM
from charmming.io.basecharmm import BaseCHARMMFile


class PRMFile(BaseCHARMMFile):
    """
    **STUB**
    """

    _sections = ['atom','bond','angl','thet','dihe','phi','impr','imph','cmap',
                'nbon','nonb','nbfi','hbon','end']

    _prmClass = {
        'atom': MassPRM,
        'angl': AnglePRM,
        'bond': BondPRM,
        'dihe': DihedralPRM,
        'impr': ImproperPRM,
        'nbon': NonBondPRM,
        'nbfi': NBFixPRM
        }

    def __init__(self, arg=None, **kwargs):
        super(PRMFile, self).__init__(arg, **kwargs)

    def _parse(self):
        """
        Reads the .prm text file, and parses the data into their relevant sections.
        """
        self.sectionCmd = {}
        self.sectionPrm = {}
        #
        iterator = ( line for line in self.body if not line.startswith('!') )
        for taco in paragraphs(iterator, self._sections):
            try:
                self.sectionCmd[taco[0][:4]] = taco[0]
                self.sectionPrm[taco[0][:4]] = taco[1:]
            except IndexError:
                pass
        # discard empties
        tmp = []
        for key, value in self.sectionPrm.iteritems():
            if not value:
                tmp.append(key)
        for key in tmp:
            del self.sectionCmd[key]
            del self.sectionPrm[key]
        # rename sections
        def rename(old, new):
            if old in self.sectionCmd.keys():
                self.sectionCmd[new] = self.sectionCmd[old]
                del self.sectionCmd[old]
                self.sectionPrm[new] = self.sectionPrm[old]
                del self.sectionPrm[old]
        rename('thet', 'angl')
        rename('phi', 'dihe')
        rename('imph', 'impr')
        rename('nonb', 'nbon')

    def objectify(self):
        """
        convert text to PRM objects
        """
        self.sectionObj = {}
        for key, Value in self._prmClass.iteritems():
            try:
                iterator = ( Value(line) for line in self.sectionPrm[key] )
                self.sectionObj[key] = sorted(list(set(iterator)))
            except KeyError:
                pass
