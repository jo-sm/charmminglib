"""
DOCME
"""
# fcp
# 10/26/2010


from charmming.tools import lowerKeys, Property
from charmming.lib.basestruct import BaseStruct
from charmming.lib.res import Res
from charmming.lib.pro import Pro


class Seg(BaseStruct):
    """
    Properties
        `addr`
        `chainid`
        `segid`
        `segType`
    Public Methods
        `iter_res`
        `reindex_atomNum`
        `reindex_resid`
    """
    def __init__(self, iterable=None, **kwargs):
        super(Seg, self).__init__(iterable, **kwargs)

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        The `addr` property provides a human readable unique string
        representation for each `Seg` instance.
        """
        def fget(self):
            return '%s.%4s' % (self.chainid, self.segType)
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
        def fset(self, value):
            for atom in self:
                atom.segType = value
        return locals()

##################
# Public Methods #
##################

    def iter_res(self,**kwargs):
        """
        A generator that returns one `Res` per iteration.

        kwargs:
                `restype`

        The kwarg `restype` allows you to specify which type of
        residue is iterated over by passing the class through
        Examples include: `Pro`, `Res` and `CGPro`.  The default
        behavior is to detect the appropriate one.

        >>> thisSeg.iter_res(restype=Pro)
        """
        # kwargs
        kwargs = lowerKeys(kwargs)
        resType = kwargs.pop('restype', None)
        if resType is None:
            if self.segType == 'pro':
                newObj = Pro
            else:
                newObj = Res
        else:
            newObj = resType
        result = newObj(iterable=None ,code=self.code, autofix=False)
        for atom in self:
            if len(result) == 0:
                result.append(atom)
                lastresid = atom.resid
            elif atom.resid == lastresid:
                result.append(atom)
            else:
                if result:
                    yield result
                result = newObj(iterable=[atom], code=self.code, autofix=False)
                lastresid = atom.resid
        if result:
            yield result

    def reindex_atomNum(self, start=1):
        """
        Reindex the atom numbers of the segment.  Indexing starts from
        the optional `start` value, which defaults to 1.
        """
        i = start
        for atom in self:
            atom.atomNum = i
            i += 1

    def reindex_resid(self, start=1):
        """
        Reindex the residue numbers of the segment.  Indexing starts from
        the optional `start` value, which defaults to 1.
        """
        i = start
        for res in self.iter_res():
            res.resid = i
            i += 1
