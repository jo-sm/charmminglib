"""
:Author: fcp
:Date: 02/22/2011
"""


from copy import deepcopy
from pychm.tools import cleanStrings, logicalLines, paragraphs
from pychm.lib.toppar import AnglePRM, BondPRM, DihedralPRM, ImproperPRM, \
        MassPRM, NonBondPRM, NBFixPRM
from pychm.io.basecharmm import BaseCHARMMFile


class PRMFile(BaseCHARMMFile):
    """
    **STUB**
    """

    _sections = ['atom','bond','angl','thet','dihe','phi','impr','imph','cmap',
                'nbon','nonb','nbfi','hbon','end']

    def __init__(self, filename=None, **kwargs):
        super(PRMFile, self).__init__(filename, **kwargs)
        self.sectionCmd = {}
        self.sectionPrm = {}
        self.comments = []
        if filename is not None:
            self.parse()

    def parse(self):
        super(PRMFile, self).parse()
        self._parse_mass()
        iterator = ( line for line in self.body if not line.lstrip().startswith('!') )
        iterator = paragraphs(logicalLines(iterator), self._sections)
        for taco in iterator:
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
        # self._parse_hbon() TODO -- implement HBondPRM objects in lib.toppar
        self._parse_cmap()

    def write(self, filename=None):
        """
        """
        def get_section(name):
            tmp = []
            try:
                tmp.append(self.sectionCmd[name])
                iterator = ( prm.Print() for prm in getattr(self, name) )
                tmp.extend(iterator)
                tmp.append('')
            except KeyError:
                pass
            return tmp
        #
        tmp = []
        # header
        iterator = ( '* %s' % line for line in self.header if line )
        tmp.extend(iterator)
        tmp.append('*')
        tmp.append('')
        #
        # atom
        i = 1
        for prm in self.atom:
            prm.index = i
            i += 1
        tmp.extend(get_section('atom'))
        tmp.extend(get_section('bond'))
        tmp.extend(get_section('angl'))
        tmp.extend(get_section('dihe'))
        tmp.extend(get_section('impr'))
        tmp.extend(get_section('nbon'))
        tmp.extend(get_section('nbfi'))

        # cmap -- OO cmap not implemented
        # TODO hbon -- hbon not yet implemented

        #
        if filename is None:
            for line in tmp:
                print line.upper()
        else:
            writeTo = open(filename, 'w')
            writeTo.write('\n'.join(tmp).upper())
            writeTo.close()

###################
# Private Methods #
###################

    def _parse_mass(self):
        tmp = [ line for line in self.body if line.startswith('mass') ]
        if tmp:
            self.sectionPrm['atom'] = tmp
            self.atom = [ MassPRM(line) for line in tmp ]

    def _parse_bond(self):
        try:
            self.bond = [ BondPRM(line) for line in self.sectionPrm['bond'] ]
        except KeyError:
            self.bond = []

    def _parse_angl(self):
        try:
            self.angl = [ AnglePRM(line) for line in self.sectionPrm['angl'] ]
        except KeyError:
            self.angl = []

    def _parse_dihe(self):
        try:
            self.dihe = [ DihedralPRM(line) for line in self.sectionPrm['dihe'] ]
        except KeyError:
            self.dihe = []

    def _parse_impr(self):
        try:
            self.impr = [ ImproperPRM(line) for line in self.sectionPrm['impr'] ]
        except KeyError:
            self.impr = []

    def _parse_nbon(self):
        try:
            self.nbon = [ NonBondPRM(line) for line in self.sectionPrm['nbon'] ]
        except KeyError:
            self.nbon = []

    def _parse_nbfi(self):
        try:
            self.nbfi = [ NBFixPRM(line) for line in self.sectionPrm['nbfi'] ]
        except KeyError:
            self.nbfi = []

    def _parse_hbon(self):
        try:
            self.hbon = [ HBondPRM(line) for line in self.sectionPrm['hbon'] ]
        except KeyError:
            self.hbon = []

    def _parse_cmap(self):
        try:
            self.cmap = self.sectionCmd['cmap']
        except KeyError:
            self.cmap = []

###################
# Special Methods #
###################

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.sectionPrm.keys())

    def __getitem__(self, key):
        return getattr(self, key)

    def __add__(self, other):
        tmp = deepcopy(self)
