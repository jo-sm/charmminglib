"""
:Author: fcp
:Date: 02/22/2011
"""


from charmming.tools import cleanStrings, logicalLines, paragraphs
from charmming.lib.toppar import AnglePRM, BondPRM, DihedralPRM, ImproperPRM, \
        MassPRM, NonBondPRM, NBFixPRM
from charmming.io.basecharmm import BaseCHARMMFile


class PRMFile(BaseCHARMMFile):
    """
    **STUB**
    """

    _sections = ['atom','bond','angl','thet','dihe','phi','impr','imph','cmap',
                'nbon','nonb','nbfi','hbon','end']

    def __init__(self, filename, **kwargs):
        super(PRMFile, self).__init__(filename, **kwargs)
        self.parse()

    def parse(self):
        super(PRMFile, self).parse()
        self.sectionCmd = {}
        self.sectionPrm = {}
        self.comments = []
        self._parse_mass()

        for taco in paragraphs(logicalLines(self.body), self._sections):
            try:
                if taco[0][:4] in self._sections:
                    self.sectionCmd[taco[0][:4]] = taco[0]
                    self.sectionPrm[taco[0][:4]] = taco[1:]
                else:
                    self.comments.append(taco)
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

        # parse sections
        self._parse_bond()
        self._parse_angl()
        self._parse_dihe()
        self._parse_impr()
        self._parse_nbon()
        self._parse_nbfi()
        self._parse_hbon()
        self._parse_cmap()

###################
# Private Methods #
###################

    def _parse_mass(self):
        tmp = []
        for line in self.body:
            if line.startswith('mass'):
                tmp.append(MassPRM(line))
        if tmp:
            self.sectionCmd['mass'] = ''
            self.sectionPrm['mass'] = tmp
        self.mass = tmp

    def _parse_bond(self):
        try:
            self.bond = [ BondPRM(line) for line in cleanStrings(
                self.sectionPrm['bond'], CC='!') ]
        except KeyError:
            self.bond = []

    def _parse_angl(self):
        try:
            self.angl = [ AnglePRM(line) for line in cleanStrings(
                self.sectionPrm['angl'], CC='!') ]
        except KeyError:
            self.angl = []

    def _parse_dihe(self):
        try:
            self.dihe = [ DihedralPRM(line) for line in cleanStrings(
                self.sectionPrm['dihe'], CC='!') ]
        except KeyError:
            self.dihe = []

    def _parse_impr(self):
        try:
            self.impr = [ ImproperPRM(line) for line in cleanStrings(
                self.sectionPrm['impr'], CC='!') ]
        except KeyError:
            self.impr = []

    def _parse_nbon(self):
        try:
            self.nbon = [ NonBondPRM(line) for line in cleanStrings(
                self.sectionPrm['nbon'], CC='!') ]
        except KeyError:
            self.nbon = []

    def _parse_nbfi(self):
        try:
            self.nbfi = [ NBFixPRM(line) for line in cleanStrings(
                self.sectionPrm['nbfi'], CC='!') ]
        except KeyError:
            self.nbfi = []

    def _parse_hbon(self):
        try:
            self.hbon = [ line for line in cleanStrings(
                self.sectionPrm['hbon'], CC='!') ]
        except KeyError:
            self.hbon = []

    def _parse_cmap(self):
        try:
            self.cmap = self.sectionCmd['cmap']
        except KeyError:
            self.cmap = []

