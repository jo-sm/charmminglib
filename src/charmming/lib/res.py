from numpy import array
from charmming.lib.basestruct import BaseStruct
from charmming.tools import Property


class Res(BaseStruct):
    """
    Properties
        chainid
        segType
        segid
        resid
        resName
        mass
        com
        heavyCom
    """
    def __init__(self,iterable=None,**kwargs):
        super(Res,self).__init__(iterable,**kwargs)

    @Property
    def chainid():
        doc = "The chainid property."
        def fget(self):
            for atom in self:
                return atom.chainid
        return locals()

    @Property
    def segType():
        doc = "The segType property."
        def fget(self):
            for atom in self:
                return atom.segType
        return locals()

    @Property
    def segid():
        doc = "The segid property."
        def fget(self):
            for atom in self:
                return atom.segid
        return locals()

    @Property
    def resid():
        doc = "The resid property."
        def fget(self):
            for atom in self:
                return atom.resid
        def fset(self, value):
            for atom in self:
                atom.resid = value
        return locals()

    @Property
    def resIndex():
        doc = "The resIndex property."
        def fget(self):
            for atom in self:
                return atom.resIndex
        def fset(self, value):
            for atom in self:
                atom.resIndex = value
        return locals()

    @Property
    def resName():
        doc = "The resName property."
        def fget(self):
            for atom in self:
                return atom.resName
        def fset(self, value):
            for atom in self:
                atom.resName = value
        return locals()

    @Property
    def addr():
        doc = "The addr property."
        def fget(self):
            return '%s.%4s.%04d' % (self.chainid,self.segType,self.resid)
        return locals()

    @Property
    def mass():
        doc = "The mass property."
        def fget(self):
            return sum( ( atom.mass for atom in self ) )
        return locals()

    @Property
    def com():
        doc = "The com property."
        def fget(self):
            tmp = array([ atom.mass * atom.cart for atom in self ])
            tmp = tmp.sum(axis=0)
            return tmp/self.mass
        return locals()

    @Property
    def heavyCom():
        doc = "The heavyCom property."
        def fget(self):
            tmp = array([ atom.mass * atom.cart for atom in self if atom.element == 'h' ])
            tmp = tmp.sum(axis=0)
            mass = sum( ( atom.mass for atom in self if atom.element == 'h' ) )
            return tmp/mass
        return locals()
