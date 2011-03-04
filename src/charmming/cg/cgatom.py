"""
DOCME
"""
# fcp
# 10/27/2010


from charmming.const.etc import alpha2num
from charmming.tools import Property, lowerKeys
from charmming.lib.metaatom import AtomError
from charmming.lib.baseatom import BaseAtom
from charmming.lib.atom import Atom


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
        `segType`
        `weight`
    Public Methods
        `is_ktgo`
        `is_bln`        TODO
        `Print`
    """

    _autoInFormat = 'charmm'
    """
    A string that defines the default input formatting for all class
    instances.
    """

    _properties = dict(BaseAtom._properties)
    _properties.update({
        'domain': 'auto',
        'structure': None
        })
    del _properties['bFactor']
    del _properties['weight']
    """
    A dictionary which tells the class which properties (keys) it is
    associated with.  The corresponding values serve as default values
    for said properties.
    """
    _sortSegType = {'ktgo':1, 'bln':2, 'bad':10}

    _tagMap = {'ktgo':'atom', 'bln':'atom', 'bad':'hetatm'}

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
            `segType`

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
            self.chainid = self._text[21]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.resIndex = self.resid
            self.domain = 'auto'
            self.segType = 'auto'
        elif inFormat == 'charmm':
            self.atomNum = self._text[6:11]
            self.atomType = self._text[12:16]
            self.resName = self._text[17:21]
            self.chainid = self._text[72:76]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.resIndex = self.resid
            self.domain = 'auto'
            self.segType = 'auto'
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
    def domain():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._domain
        def fset(self, value):
            value = str(value).strip().lower()
            if value == 'auto':
                self._domain = alpha2num[self.chainid]
            else:
                self._domain = value
        return locals()

    @Property
    def prmString():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return '%s%d%s' % (self.chainid, self.resid, self.atomType)
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
        """
        return '%s%d' % (self.chainid, self.resid) == self.resName

    def is_bln(self):
        """
        Logic for determining if an atom was constructed using the BLN
        (hydrophobic/hydrophillic/apolar) model, returns a bool.

        TODO
        """
        pass

    def Print(self, **kwargs):
        """
        Returns a string representation of the Atom object. Formatting
        defaults to the `_autoInFormat` formatting.

        kwargs:
            `outformat`     ["charmm","debug","xdebug","crd","xcrd"]
            `old_chainid`   [False,True]
            `old_segtype`   [False,True]
            `old_resid`     [False,True]
            `old_atomnum`   [False,True]

        >>> thisAtom.Print(outformat="charmm",old_resid=True)
        """
        # Make kwargs case insensitive
        kwargs = lowerKeys(kwargs)
        # kwargs
        outFormat = kwargs.pop('outformat', self.__class__._autoInFormat)
        old_chainid = kwargs.pop('old_chainid', False)
        old_segType = kwargs.pop('old_segtype', False)
        old_resid = kwargs.pop('old_resid', False)
        old_atomNum = kwargs.pop('old_atomnum', False)
        # Localize variables
        x, y, z = self.cart
        if old_chainid:
            chainid = self.chainid0
        else:
            chainid = self.chainid
        if old_segType:
            segType = self.segType0
        else:
            segType = self.segType
        if old_resid:
            resid = self.resid0
        else:
            resid = self.resid
        if old_atomNum:
            atomNum = self.atomNum0
        else:
            atomNum = self.atomNum

        # Set formatting
        if outFormat == 'charmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s' % \
                (self.tag, atomNum, self.atomType, self.resName, resid,
                x, y, z, self.weight, self.bFactor, chainid)
        elif outFormat == 'debug':
            taco = '%-6s%5i %4s %4s %1s%4i    %8.3f%8.3f%8.3f' % \
                (segType, atomNum, self.atomType, self.resName, chainid, resid,
                x, y, z)
        elif outFormat == 'xdebug':
            taco = '%15s > %15s: %5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f' % \
                (self.addr0, self.addr, atomNum, self.atomType, self.resName,
                chainid, resid, x, y, z)
        elif outFormat == 'repr':
            taco = '%15s :: %4s %4s    %8.3f %8.3f %8.3f' % \
                (self.addr, self.atomType, self.resName, x, y, z)
        elif outFormat in ['crd', 'cor', 'card', 'short', 'shortcard']:
            taco = '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f' % \
                (atomNum, self.resIndex, self.resName, self.atomType,
                x, y, z, chainid, resid, self.weight)
        elif outFormat in ['xcrd', 'xcor', 'xcard', 'long', 'longcard']:
            taco = '%10i%10i  %-4s      %-4s    %20.10f%20.10f%20.10f  %-5s    %-4i    %20.10f' % \
                (atomNum, self.resIndex, self.resName, self.atomType,
                x, y, z, chainid + segType, resid, self.weight)
        else:
            raise AtomError('Print: unknown inFormat %s' % inFormat)
        return taco.upper()

