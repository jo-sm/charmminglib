"""
DOCME
"""
# fcp
# 10/27/2010


from charmming.lib.baseatom import BaseAtom
from charmming.lib.metaatom import AtomError
from charmming.tools import Property


class CGAtom(BaseAtom):
    """
    The standard CHARMMing implementation of the coarse-grained
    resolution atom-like object.

    This object does not vary much from the canonical all-atom `Atom`
    object, save for removal of the element Property, and some
    compliance methods.

    Class Attributes
        `_autoInFormat`
        `_sortSegType`
        `_tagMap`
    Parse Definition
        `parse`
    Properties
        `bFactor`
        `tag`
        `weight`
    Public Methods
        `is_ktgo`       TODO
        `is_bln`        TODO
        `Print`
    """

    _autoInFormat = 'charmm'

    _sortSegType = {'ktgo':1, 'bln':2}

    _tagMap = {'ktgo':'atom', 'bln':'atom'}

    def __init__(self,text=None,**kwargs):
        """
        DOCME
        """
        super(CGAtom, self).__init__(text, **kwargs)

####################
# Parse Definition #
####################

    def parse(self, inFormat):
        """
        Parses _text into the following instance variables:
            `atomNum`
            `atomType`
            `resName`
            `chainid`
            `resid`
            `cart`
            `resIndex`

        Initial values for: `chainid`, `segType`, `resid` and `atomNum`
        are saved as well as `chainid0`, etc.

        Recognized `inFormat` values are:
            `"pdborg"`
            `"charmm"`
        """
        if inFormat == 'pdborg':
            self.atomNum = self._text[6:11]
            self.atomType = self._text[12:16]
            self.resName = self._text[16:20]
            self.segType = 'auto'
            self.chainid = self._text[21]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.resIndex = self.resid
        elif inFormat == 'charmm':
            self.atomNum = self._text[6:11]
            self.atomType = self._text[12:16]
            self.resName = self._text[17:21]
            self.segType = 'auto'
            self.chainid = self._text[72:76]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.resIndex = self.resid
        else:
            raise AtomError('parse: unknown `inFormat`: %s' % inFormat)
        # Save initial properties
        self._chainid0 = self.chainid
        self._segType0 = self.segType
        self._resid0 = self.resid
        self._atomNum0 = self.atomNum

##############
# Properties #
##############

    @Property
    def bFactor():
        doc =\
        """
        Always returns a float(0), as `bFactor`s are not currently
        implemented in coarse grained models.
        """
        def fget(self):
            return 0.
        return locals()

@Property
def segType():
    doc =\
    """
    DOCME
    """
    def fget(self):
        return self._segType
    def fset(self, value):
        value = str(value).strip().lower()
        if value == 'auto':
            if self.is_ktgo():
                self._segType = 'ktgo'
            elif self.is_bln():
                self._segType = 'bln'
            else:
                self._segType = 'bad'
        else:
            self._segType = value
    return locals()

    @Property
    def weight():
        doc =\
        """
        Always returns a float(0), as `weight`s are not currently
        implemented in coarse grained models.
        """
        def fget(self):
            return 0.
        return locals()

##################
# Public Methods #
##################

    def is_ktgo(self):
        """
        Logic for determining if an atom was constructed using the KT
        Go model, returns a bool.

        TODO
        """
        return self.atomType in ['   B','   S']

    def is_bln(self):
        """
        Logic for determining if an atom was constructed using the BLN
        (hydrophobic/hydrophillic/apolar) model, returns a bool.

        TODO
        """
        pass

    def Print(self, **kwargs):
        """
        DOCME
        """
        pass
