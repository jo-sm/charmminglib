"""
DOCME
"""
# fcp
# 10/26/2010


from charmming.lib.basestruct import BaseStruct
from charmming.lib.seg import Seg
from charmming.tools import Property, lowerKeys


class Chain(BaseStruct):
    """
    Properties
        `addr`
        `chainid`
    Public Methods
        `iter_res`
        `iter_seg`
    """
    def __init__(self, iterable=None, **kwargs):
        super(Chain, self).__init__(iterable, **kwargs)

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        The `addr` property provides a human readable unique string
        representation for each `Chain` instance.
        """
        def fget(self):
            return  '%s' % self.chainid
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
        def fset(self, value):
            for atom in self:
                atom.chainid = value
        return locals()

##################
# Public Methods #
##################

    def iter_res(self,**kwargs):
        """
        A generator that returns one `Res` per iteration.

        kwargs:
            `segtypes`

        The kwarg `segtypes` can be used to specify which residues are
        iterated over, and in which order.  By default, all residues
        are iterated over in alpha order of `segtype`.

        >>> thisChain.iter_res(segtypes=['pro','good'])
        """
        for seg in self.iter_seg(**kwargs):
            for res in seg.iter_res():
                yield res

    def iter_seg(self, **kwargs):
        """
        A generator that returns one `Seg` per iteration.

        kwargs:
            `segtypes`

        The kwarg `segtypes` can be used to specify which segments are
        iterated over, and in which order.  By default, all segments
        are iterated over in alpha order of `segtype`.

        >>> thisChain.iter_seg(segtypes=['pro','dna','rna','nuc','good','bad'])
        """
        # kwargs
        kwargs = lowerKeys(kwargs)
        segTypes = kwargs.pop('segtypes', None)
        # Default to 'all' segtypes
        if segTypes is None:
            segTypes = list(set(( atom.segType for atom in self )))
            segTypes.sort()
        # Do Work
        for segType in segTypes:
            iterator = ( atom for atom in self if atom.segType == segType )
            result = Seg(iterable=iterator,code=self.code,autoFix=False)
            if result:
                yield result
