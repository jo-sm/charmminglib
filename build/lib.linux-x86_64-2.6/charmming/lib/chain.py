from charmming.lib.basestruct import BaseStruct
from charmming.lib.seg import Seg
from charmming.tools import Property


class Chain(BaseStruct):
    def __init__(self,iterable=None,**kwargs):
        super(Chain,self).__init__(iterable,**kwargs)

    @Property
    def chainid():
        doc = "The chainid property."
        def fget(self):
            for atom in self:
                return atom.chainid
        def fset(self, value):
            for atom in self:
                atom.chainid = value
        return locals()

    @Property
    def addr():
        doc = "The addr property."
        def fget(self):
            return  '%s' % self.chainid
        return locals()

    def iter_seg(self,**kwargs):
        """
        For looping over segments.
        kwargs:
            segTypes
        >>> thisChain.iter_seg(segTypes=['pro','dna','rna','nuc','good','bad'])
        """
        # kwargs
        kwargs.pop('chainids',None)
        segTypes = kwargs.pop('segTypes',None)
        if kwargs.keys():
            raise TypeError('Unprocessed kwargs(%r)' % kwargs.keys())
        #
        if segTypes is None:
            segTypes = list(set(( atom.segType for atom in self )))
            segTypes.sort()
        # Do Work
        for segType in segTypes:
            iterator = ( atom for atom in self if atom.segType == segType )
            result = Seg(iterable=iterator,code=self.code,autoFix=False)
            if result:
                yield result

    def iter_res(self,**kwargs):
        """
        For looping over residues.
        kwargs:
            segTypes
        >>> thisChain.iter_res(segTypes=['pro','good'])
        """
        for seg in self.iter_seg(**kwargs):
            for res in seg.iter_res():
                yield res
