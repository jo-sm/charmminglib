from charmming.lib.basestruct import BaseStruct
from charmming.lib.pro import Pro
from charmming.lib.res import Res
from charmming.tools import Property

class Seg(BaseStruct):
    """
    Properties
        chainid
        segType
        segid
    """
    def __init__(self,iterable=None,**kwargs):
        """"""
        super(Seg,self).__init__(iterable,**kwargs)

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
        def fset(self, value):
            for atom in self:
                atom.segType = value
        return locals()

    @Property
    def segid():
        doc = "The segid property."
        def fget(self):
            for atom in self:
                return atom.segid
        return locals()

    @Property
    def addr():
        doc = "The addr property."
        def fget(self):
            return '%s.%4s' % (self.chainid,self.segType)
        return locals()

    def iter_res(self):
        """
        For looping over residues.
        >>> thisSeg.iter_res()
        """
        if self.segType == 'pro':
            newObj = Pro
        else:
            newObj = Res
        result = newObj(iterable=None,code=self.code,autoFix=False)
        for atom in self:
            if len(result) == 0:
                result.append(atom)
                lastresid = atom.resid
            elif atom.resid == lastresid:
                result.append(atom)
            else:
                if result:
                    yield result
                result = newObj(iterable=[atom],code=self.code,autoFix=False)
                lastresid = atom.resid
        if result:
            yield result

    def reindex_atomNum(self,start=1):
        """
        Reindex the atom numbers of the segment.  Indexing starts from the
        optional 'start' value, which defaults to 1.
        """
        i = start
        for atom in self:
            atom.atomNum = i
            i += 1

    def reindex_resid(self,start=1):
        """
        Reindex the residue numbers of the segment.  Indexing starts from the
        optional 'start' value, which defaults to 1.
        """
        i = start
        for res in self.iter_res():
            res.resid = i
            i += 1
