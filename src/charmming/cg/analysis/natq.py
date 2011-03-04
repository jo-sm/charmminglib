

from charmming.tools import Property
from charmming.analysis.baseanalysis import BaseAnalysis,load_correlOutput
from charmming.io.pdb import PDBFile
from charmming.cg.ktgo import KTGo


class NatQ(BaseAnalysis):
    """
    docstring for NatQ
    """
    def __init__(self, arg=None):
        super(NatQ, self).__init__(arg)
        #
        if arg is not None:
            self.aa = PDBFile(self.pdbFilename)[0]
            self.aa.parse()
            self.cg = KTGo(self.aa)
            self.natq = self.cg.get_nativeSCSC()

