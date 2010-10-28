

from charmming.analysis.cg.basecganalysis import BaseAnalysis,load_correlOutput
from charmming.tools import Property


class NatQ(BaseAnalysis):
    """
    docstring for NatQ
    """
    def __init__(self, aaCrdFileName):
        super(NatQ, self).__init__(aaCrdFileName)

    @Property
    def all():
        doc = "The all property."
        def fget(self):
            return self._all
        def fset(self, value):
            self._all = value
        return locals()


# GoAtom class too
# WRITE THE CG CLASS ALREADY

# CG .pdb representation
# - Write cg.pdb file
# - Parse into cgMol object

# Determine all native contacts
# - List all contacts
# --- cgAtom.atomType == '   S'
# --- i < j pairs of cgAtoms
# - Wean then to all native contacts
# --- pairs which are <= 10A apart
# --- If any atoms between 2 residues are < 8A -> native
