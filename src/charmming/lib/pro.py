"""
DOCME
"""
# fcp
# 10/26/2010


from copy import deepcopy
from numpy import array
from charmming.const.bio import aaMass
from charmming.lib.res import Res
from charmming.tools import Property


class NoAlphaCarbonError(Exception):
    """
    Exception to raise when an alpha carbon is expected but not found.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


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
        `sanity`            STUB
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
                raise NoAlphaCarbonError()
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
        raise NoAlphaCarbonError()

    def get_goBB(self):
        """
        DOCME
        """
        result = self.get_alphaCarbon()
        result.tag = 'atom'
        result.atomNum = result.resid
        result.atomType = '   b'
        result.weight = 0
        result.bFactor = 0
        result.mass = aaMass['gly']
        result.segType = 'ktgo'
        result.addr0 = result.addr
        # TODO Other things need checking
        return result

    def get_goSC(self):
        """
        DOCME
        """
        if self.resName == 'gly':
            raise NoAlphaCarbonError()
        result = self.get_alphaCarbon()
        result.tag = 'atom'
        result.atomNum = self.resid
        result.atomType = '   s'
        result.cart = self.scCom
        result.weight = 0
        result.bFactor = 0
        result.mass = aaMass[result.resName] - aaMass['gly']
        result.segType = 'ktgo'
        result.addr0 = result.addr
        del result.element
        # TODO Other things need checking
        return result

    def sanity(self):
        raise NotImplementedError

