"""
DOCME
"""
# fcp
# 01/25/2011


from pychm.tools import lowerKeys
from pychm.const.bio import aaMass
from pychm.lib.pro import NoAlphaCarbonError, ProError, Pro
from pychm.cg.cgatom import CGAtom


def isBBAA(arg):
    """
    Returns a bool, `True` if the atom in an all-atom representaiton would map
    into the backbone particle
    """
    return arg in set([' n  ', ' hn ', ' ca ', ' ha ', ' ha1', ' ha2', ' c  ',
                    ' o  '])


class CGPro(Pro):
    """
    DOCME
    """
    def __init__(self, iterable=None, **kwargs):
        """
        DOCME
        """
        super(CGPro, self).__init__(iterable, **kwargs)

##################
# Public Methods #
##################

    def get_alphaCarbon(self):
        """
        """
        for atom in self:
            if atom.atomType == ' ca ':
                return CGAtom(atom)
        raise NoAlphaCarbonError

    def get_goBB(self, **kwargs):
        """
        DOCME
        """
        kwargs = lowerKeys(kwargs)
        calcmass = kwargs.get('calcmass', None)
        #
        tmp = self.get_alphaCarbon()
        if calcmass is None:
            tmp.mass = aaMass['gly']
        else:
            if self.resName == 'gly':
                iterator = ( atom.mass for atom in self )
            else:
                iterator = ( atom.mass for atom in self if isBBAA(atom.atomType) )
            tmp.mass = sum(iterator)
        tmp.atomNum = 0
        tmp.atomType = 'b'
        tmp.derivedResName = tmp.resName
        tmp.resName = '%s%d' % (self.chainid, self.resid)
        tmp.segType = 'ktgo'
        return tmp

    def get_goSC(self, **kwargs):
        """
        DOCME
        """
        kwargs = lowerKeys(kwargs)
        calcmass = kwargs.get('calcmass', None)
        #
        tmp = self.get_alphaCarbon()
        tmp.cart = self.scCom
        if calcmass is None:
            tmp.mass = aaMass[self.resName] - aaMass['gly']
        else:
            iterator = ( atom.mass for atom in self if not isBBAA(atom.atomType) )
            tmp.mass = sum(iterator)
        tmp.atomNum = 1
        tmp.atomType = 's'
        tmp.derivedResName = tmp.resName
        tmp.resName = '%s%d' % (self.chainid, self.resid)
        tmp.segType = 'ktgo'
        return tmp

    def get_blnBB(self):
        pass

    def get_blnSC(self):
        pass

