"""
DOCME
"""
# fcp
# 10/26/2010


from numpy import array
from pychm.tools import Property
from pychm.lib.basestruct import BaseStruct


class Res(BaseStruct):
    """
    Properties
        `addr`
        `chainid`
        `heavyCom`
        `resid`
        `resIndex`
        `resName`
        `segid`
        `segType`
    """
    def __init__(self, iterable=None, **kwargs):
        super(Res, self).__init__(iterable, **kwargs)

##################
# Public methods #
##################

    def iter_atom(self):
        doc =\
        """
        A generator that returns each atom in the residue
        """
        for atom in self:
            yield atom

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        The `addr` property provides a human readable unique string
        representation for each `Res` instance.
        """
        def fget(self):
            return '%s.%s.%d' % (self.chainid, self.segType, self.resid)
        return locals()

    @Property
    def chainid():
        doc =\
        """
        DOCME
        """
        def fget(self):
            for atom in self:
                return atom.chainid
        return locals()

    @Property
    def heavyCom():
        doc =\
        """
        The center of mass of residue, computed using only "heavy"
        (non-hydrogen) atoms.

        Care should be taken with this method, as it filters with
        `BaseAtom.element` which won't necesarily be defined for all
        atoms.
        """
        def fget(self):
            result = array([ atom.mass * atom.cart for atom in self
                        if atom.element == 'h' ])
            result = result.sum(axis=0)
            mass = sum( ( atom.mass for atom in self if atom.element == 'h' ) )
            return result / mass
        return locals()

    @Property
    def resid():
        doc =\
        """
        DOCME
        """
        def fget(self):
            for atom in self:
                return atom.resid
        def fset(self, value):
            for atom in self:
                atom.resid = value
        return locals()

    @Property
    def resIndex():
        doc =\
        """
        DOCME
        """
        def fget(self):
            for atom in self:
                return atom.resIndex
        def fset(self, value):
            for atom in self:
                atom.resIndex = value
        return locals()

    @Property
    def resName():
        doc =\
        """
        DOCME
        """
        def fget(self):
            for atom in self:
                return atom.resName
        def fset(self, value):
            for atom in self:
                atom.resName = value
        return locals()

    @Property
    def segid():
        doc =\
        """
        DOCME
        """
        def fget(self):
            for atom in self:
                return atom.segid
        return locals()

    @Property
    def segType():
        doc =\
        """
        DOCME
        """
        def fget(self):
            for atom in self:
                return atom.segType
        return locals()

    @property
    def resid0():
        return self.iter_atom().next().resid0

###################
# private methods #
###################

    def _dogmans_rename(self):
        elements = set(( atom.element for atom in self ))
        for this_element in elements:
            tmp = [ atom for atom in self if atom.element == this_element ]
            if len(tmp) == 1:
                tmp[0].atomType = this_element
            else:
                for i, atom in enumerate(tmp):
                    atom.atomType = "%s%d" % (this_element, i+1)
