"""
DOCME
"""
# fcp
# 10/26/2010


from copy import deepcopy
from numpy import array
from charmming.lib.res import Res
from charmming.tools import Property


class NoAlphaCarbonError(Exception):
    """
    Exception to raise when an alpha carbon is expected but not found.
    """
    def __init__(self):
        super(NoAlphaCarbonError, self).__init__()


class ProError(Exception):
    """
    Exception to raise when errors occur involving the Pro class.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Pro(Res):
    """
    Properties
        `scCom`
    Public Methods
        `get_alphaCarbon`
        `get_goBB`
        `get_goSC`
        `sanity`            TODO
    """
    def __init__(self, iterable=None, **kwargs):
        super(Pro, self).__init__(iterable, **kwargs)

##############
# Properties #
##############

    @Property
    def scCom():
        doc =\
        """
        The center of mass of the resdiue's "side chain" atoms.
        """
        def fget(self):
            if self.resName == 'gly':
                raise NoAlphaCarbonError
            scList = [ atom for atom in self
                    if atom.element != 'h' and not atom.is_backbone() ]
            result = array([ atom.mass * atom.cart for atom in scList ])
            result = result.sum(axis=0)
            mass = sum( ( atom.mass for atom in scList ) )
            return result / mass
        return locals()

##################
# Public Methods #
##################

    def get_alphaCarbon(self):
        """
        Return a residue's alpha carbon as a new Atom object, alternatively
        raises a NoAlphaCarbonError.
        """
        for atom in self:
            if atom.atomType == ' ca ':
                return copy.deepcopy(atom)
        raise NoAlphaCarbonError

    def iter_bbAtoms(self):
        """
        """
        for atom in self:
            if atom.is_backbone():
                yield atom

    def iter_scAtoms(self):
        """
        """
        for atom in self:
            if not atom.is_backbone():
                yield atom

    def sanity(self):
        raise NotImplementedError
