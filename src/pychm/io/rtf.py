"""
:Author: fcp

**STUB**
"""


from pychm.lib.toppar import MassPRM, TopRes
from pychm.tools import logicalLines, paragraphs
from pychm.io.basecharmm import BaseCHARMMFile


class RTFFile(BaseCHARMMFile):
    """
    DOCME
    """

    _sections = ['resi', 'pres', 'end']

    def __init__(self, filename=None, **kwargs):
        """
        DOCME
        """
        super(RTFFile, self).__init__(filename, **kwargs)
        self.parse()

    def parse(self):
        """
        DOCME
        """
        super(RTFFile, self).parse()
        self._parse_mass()
        self._parse_decl()
        self._parse_defa()
        self._parse_auto()
        iterator = ( line for line in self.body if not line.startswith('!') )
        iterator = paragraphs(logicalLines(iterator), self._sections)
        self.resi = {}
        self.pres = {}
        self._wtf = []
        try:
            for taco in iterator:
                if taco[0][:4] == 'resi':
                    tmp = taco[0].split()[1]
                    self.resi[tmp] = TopRes(taco)
                elif taco[0][:4] == 'pres':
                    tmp = taco[0].split()[1]
                    self.pres[tmp] = TopRes(taco)
                else:
                    self._wtf.append(taco)
        except IndexError:
            print taco
            raise IndexError
        self.chemResiDict = dict(( (k, v.chemDict) for k, v in self.resi.iteritems() ))
        self.chemPresDict = dict(( (k, v.chemDict) for k, v in self.pres.iteritems() ))
        self.chargeResiDict = dict(( (k, v.chargeDict) for k, v in self.resi.iteritems() ))
        self.chargePresDict = dict(( (k, v.chargeDict) for k, v in self.pres.iteritems() ))

    def _parse_mass(self):
        self.atom = []
        tmp = []
        for line in self.body:
            if line.startswith('mass'):
                self.atom.append(line)
            else:
                tmp.append(line)
        self.atom = map(MassPRM, self.atom)
        self.body = tmp

    def _parse_decl(self):
        self.decl = []
        tmp = []
        for line in self.body:
            if line.startswith('decl'):
                self.decl.append(line)
            else:
                tmp.append(line)
        self.body = tmp

    def _parse_defa(self):
        self.defa = []
        tmp = []
        for line in self.body:
            if line.startswith('defa'):
                self.defa.append(line)
            else:
                tmp.append(line)
        self.body = tmp

    def _parse_auto(self):
        self.auto = []
        tmp = []
        for line in self.body:
            if line.startswith('auto'):
                self.auto.append(line)
            else:
                tmp.append(line)
        self.body = tmp

