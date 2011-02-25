"""
DOCME
"""
# fcp
# 01/25/2011


from charmming.cg.cgatom import CGAtom
from charmming.const.bio import aaMass
from charmming.lib.pro import NoAlphaCarbonError, ProError, Pro


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

    def get_goBB(self):
        """
        DOCME
        """
        tmp = self.get_alphaCarbon()
        tmp.mass = aaMass['gly']
        tmp.atomType = 'b'
        tmp.derivedResName = tmp.resName
        tmp.resName = '%s%d' % (self.chainid, self.resid)
        tmp.segType = 'ktgo'
        return tmp

    def get_goSC(self):
        """
        DOCME
        """
        tmp = self.get_alphaCarbon()
        tmp.cart = self.scCom
        tmp.mass = aaMass[self.resName]
        tmp.atomType = 's'
        tmp.derivedResName = tmp.resName
        tmp.resName = '%s%d' % (self.chainid, self.resid)
        tmp.segType = 'ktgo'
        return tmp

    def get_blnBB(self):
        pass

    def get_blnSC(self):
        pass

