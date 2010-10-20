from copy import deepcopy
from numpy import array
from charmming.const.bio import aaMass
from charmming.lib.res import Res
from charmming.tools import Property


class NoAlphaCarbonError(Exception):
    """
    Exception to raise when an alpha carbon is expected but not found.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class ProError(Exception):
    """
    Exception to raise when errors occur involving the Pro class.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Pro(Res):
    """
    Properties
        scCom
    Public Methods
        sanity              STUB
        get_alphaCarbon
        get_goBB
        get_goSC
    """
    def __init__(self,iterable=None,**kwargs):
        super(Pro,self).__init__(iterable,**kwargs)

# Properties
    @Property
    def scCom():
        doc = "The scCom property."
        def fget(self):
            if self.resName == 'gly':
                raise NoAlphaCarbonError()
            scList = [ atom for atom in self if atom.element != 'h' and not atom.isBackBone() ]
            tmp = array([ atom.mass * atom.cart for atom in scList ])
            tmp = tmp.sum(axis=0)
            mass = sum( ( atom.mass for atom in scList ) )
            return tmp/mass
        return locals()

# Public Methods
    def sanity(self):
        NotImplementedError

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
        tmp = self.get_alphaCarbon()
        tmp.tag        = 'atom'
        tmp.atomNum    = tmp.resid
        tmp.atomType   = '   b'
        tmp.weight     = 0
        tmp.bFactor    = 0
        tmp.mass       = Aux.aaMass['gly']
        tmp.segType    = 'ktgo'
        tmp.addr0      = tmp.addr
        del tmp.element
        # TODO Other things need checking
        return tmp

    def get_goSC(self):
        if self.resName == 'gly':
            raise ProError('get_alphaCarbon: No Alpha Carbon')
        tmp = self.get_alphaCarbon()
        tmp.tag        = 'atom'
        tmp.atomNum    = self.resid
        tmp.atomType   = '   s'
        tmp.cart       = self.scCom
        tmp.weight     = 0
        tmp.bFactor    = 0
        tmp.mass       = Aux.aaMass[tmp.resName] - Aux.aaMass['gly']
        tmp.segType    = 'ktgo'
        tmp.addr0      = tmp.addr
        del tmp.element
        # TODO Other things need checking
        return tmp

