"""
:Author: fcp
:Date: 10/22/2010
"""


from charmming.const.bio import atomMass, backbone, charmm2pdbAtomNames, good, \
        pro, nuc
from charmming.tools import Property, lowerKeys
from charmming.lib.metaatom import AtomError
from charmming.lib.baseatom import BaseAtom


class Atom(BaseAtom):
    """
    The standard CHARMMing implementation of the all-atom resolution
    atom-like object.

    A longer
    multi-line
    description.

    Class Attributes
        `_autoInFormat`
        `_sortSegType`
        `_tagMap`
    Parse Definition
        `parse`
    Properties
        `element`
        `segType`
    Public Methods
        `derive_element`
        `is_backbone`
        `is_good`
        `is_pro`
        `is_nuc`
        `Print`
    Private Methods
        `_compliance_resName`
        `_compliance_atomType`
        `_init_null`
    """

    _autoInFormat = 'pdborg'
    """
    A string that defines the default input formatting for all class
    instances.
    """

    _properties = dict(BaseAtom._properties)
    _properties.update({
        'element': 'x'
        })
    """
    A dictionary which tells the class which properties (keys) it is
    associated with.  The corresponding values serve as default values
    for said properties.
    """

    _sortSegType = {'pro':1,'dna':2,'rna':3,'nuc':4,'good':5,'bad':10}
    """
    This is a dictionary which defines the recognized `segTypes` (keys)
    for the `_sort` scoring method, and their weighting (values).
    """

    _tagMap = {'nuc':'atom','pro':'atom','good':'atom','bad':'hetatm',
            'dna':'atom','rna':'atom','bln':'atom'}
    """
    Map `segType` -> `tag`

    This is a dictionary which defines the recognized `segTypes` (keys)
    for the `tag` property.
    """

    def __init__(self, text=None, **kwargs):
        """
        DOCME
        """
        super(Atom, self).__init__(text, **kwargs)

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
            `weight`
            `bFactor`
            `element`
            `resIndex`

        Initial values for: `chainid`, `segType`, `resid` and `atomNum`
        are saved as well as `chainid0`, etc.

        Recognized `inFormat` values are:
            `"pdborg"`
            `"charmm"`
            `"shortcard"`
            `"longcard"`
            `"amber"`   TODO
        """
        if inFormat == 'pdborg':
            self.atomNum = self._text[6:11]
            self.atomType = self._text[12:16]
            self.resName = self._text[16:20]
            self.chainid = self._text[21]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.weight = self._text[55:60]
            self.bFactor = self._text[61:66]
            self.element = self._text[66:]
            self.resIndex = self.resid
            self.segType = 'auto'
        elif inFormat == 'charmm':
            self.atomNum = self._text[6:11]
            self.atomType = self._text[12:16]
            self.resName = self._text[17:21]
            self.chainid = self._text[72:76]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.weight = self._text[55:60]
            self.bFactor = self._text[61:66]
            self.resIndex = self.resid
            self.segType = 'auto'
            self.element = 'auto'   # this line must come after `self.segType`
        elif inFormat in ['crd', 'cor', 'card', 'short', 'shortcard']:
            self.atomNum = self._text[0:5]
            self.resid = self._text[5:10]
            self.resName = self._text[11:15]
            self.atomType = self._text[16:20]
            self.cart = (self._text[20:30], self._text[30:40], self._text[40:50])
            self.chainid = self._text[51:55]
            self.weight = 0
            self.bFactor = 1.
            self.resIndex = self.resid
            self.segType = 'auto'
            self.element = 'auto'   # this line must come after `self.segType`
        # TODO: verify these against the format statements in CHARMM
        elif inFormat in ['xcrd', 'xcor', 'xcard', 'long', 'longcard']:
            self.atomNum = self._text[0:10]
            self.resid = self._text[10:20]
            self.resName = self._text[22:27]
            self.atomType = self._text[32:36] # NB, I am leaving this @ 4 chars for now for format compatibility
            self.cart = (self._text[40:60], self._text[60:80], self._text[80:100])
            self.chainid = self._text[102:107]
            self.weight = 0
            self.bFactor = 1.
            self.resIndex = self.resid
            self.segType = 'auto'
            self.element = 'auto'   # this line must come after `self.segType`
        elif inFormat == 'amber':
            raise NotImplementedError
        else:
            raise AtomError('parse: unknown `inFormat`: %s' % inFormat)
        #
        try:
            if self.element == 'x':
                self.element = self.derive_element()
        except AttributeError:
            pass
        # Save initial properties
        self._chainid0 = self.chainid
        self._segType0 = self.segType
        self._resid0 = self.resid
        self._atomNum0 = self.atomNum

##############
# Properties #
##############

    @Property
    def element():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._element
        def fset(self, value):
            value = str(value).strip().lower()
            if value == 'auto':
                self._element = self.derive_element()
                self.mass = atomMass[self._element]
            else:
                if value in atomMass.keys():
                    self._element = value
                    self.mass = atomMass[value]
                elif value in charmm2pdbAtomNames.keys():
                    self._element = charmm2pdbAtomNames[value]
                    self.mass = atomMass[self._element]
                elif self._autoFix:
                    self._element = 'x'
                    self.mass = 0
                else:
                    raise AtomError('element %s: %s is invalid' %
                                    (self.addr0, value))
        def fdel(self):
            del self._element
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
                if self.is_pro():
                    self._segType = 'pro'
                elif self.is_nuc():
                    self._segType = 'nuc'
                elif self.is_good():
                    self._segType = 'good'
                else:
                    self._segType = 'bad'
            else:
                self._segType = value
        return locals()

##################
# Public Methods #
##################

    def derive_element(self):
        """
        DOCME
        """
        if self.segType in ['pro', 'nuc', 'rna', 'dna']:
            return self._atomType.strip()[0]
        if self.segType == 'good' and self.atomType in charmm2pdbAtomNames.keys():
            return charmm2pdbAtomNames[self.atomType]
        else:
            return self.atomType[:2].strip()

    def is_backbone(self):
        """
        Determines if the atom is part of the backbone structure using
        the `atomType` property.  Returns a bool.
        """
        return self.atomType in backbone

    def is_good(self):
        """
        Determines if the atom is a "good" hetero atom using the
        `resName` property.  Returns a bool.
        """
        return self.resName in good

    def is_pro(self):
        """
        Determines if the atom is part of a protein structure using the
        `resName` property.  Returns a bool.
        """
        return self.resName[-3:] in pro

    def is_nuc(self):
        """
        Determines if the atom is part of a nucleic acid structure
        using the `resName` property.  Returns a bool.
        """
        return self.resName[-3:] in nuc

    def Print(self, **kwargs):
        """
        Returns a string representation of the Atom object. Formatting
        defaults to the `_autoInFormat` formatting.

        kwargs:
            `outformat`     ["pdborg","charmm","debug","xdebug","crd","xcrd"]
            `old_chainid`   [False,True]
            `old_segtype`   [False,True]
            `old_resid`     [False,True]
            `old_atomnum`   [False,True]

        >>> thisAtom.Print(outformat="pdborg",old_resid=True)
        """
        # Make kwargs case insensitive
        kwargs = lowerKeys(kwargs)
        # kwargs
        outFormat = kwargs.get('outformat', self.__class__._autoInFormat)
        old_chainid = kwargs.get('old_chainid', False)
        old_segType = kwargs.get('old_segtype', False)
        old_resid = kwargs.get('old_resid', False)
        old_atomNum = kwargs.get('old_atomnum', False)
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
        if outFormat == 'pdborg':
            taco = '%-6s%5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s' % \
                (self.tag, atomNum, self.atomType, self.resName, chainid, resid,
                x, y, z, self.weight, self.bFactor, self.element)
        elif outFormat == 'charmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s' % \
                (self.tag, atomNum, self.atomType, self.resName, resid,
                x, y, z, self.weight, self.bFactor, chainid)
        elif outFormat == 'debug':
            taco = '%-6s%5i %4s %4s %1s%4i    %8.3f%8.3f%8.3f %3s' % \
                (segType, atomNum, self.atomType, self.resName, chainid, resid,
                x, y, z, self.element)
        elif outFormat == 'xdebug':
            taco = '%15s > %15s: %5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s' % \
                (self.addr0, self.addr, atomNum, self.atomType, self.resName,
                chainid, resid, x, y, z, self.weight, self.bFactor, self.element)
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

###################
# Private Methods #
###################

    def _compliance_resName(self):
        """
        Re-label `resName` string to be CHARMM compliant.
        """
        if self.resName in ['hoh','tip3']:
           self.resName = 'tip3'
        elif self.resName == 'his':
           self.resName = 'hsd'
        elif self.resName in ['a','da']:
           self.resName = 'ade'
        elif self.resName in ['t','dt']:
           self.resName = 'thy'
        elif self.resName in ['c','dc']:
           self.resName = 'cyt'
        elif self.resName in ['g','dg']:
           self.resName = 'gua'
        elif self.resName in ['u','du']:
           self.resName = 'ura'
        elif (self.atomType == 'zn  ' and self.resName == 'zn'):
           self.resName = 'zn2'
        elif (self.atomType == 'na  ' and self.resName == 'na'):
           self.resName = 'sod'
        elif (self.atomType == 'cs  ' and self.resName == 'cs'):
           self.resName = 'ces'
        elif (self.atomType == 'cl  ' and self.resName == 'cl'):
           self.resName = 'cla'
        elif (self.atomType == 'ca  ' and self.resName == 'ca'):
           self.resName = 'cal'
        elif (self.atomType == ' k  ' and self.resName == 'k'):
           self.resName = 'pot'

    def _compliance_atomType(self):
        """
        Re-label `atomType` string to be CHARMM compliant.
        """
        if self.resName in ['hoh','tip3']:
           self.atomType = ' oh2'
        elif (self.resName == 'ile' and self.atomType == ' cd1'):
           self.atomType = ' cd '
        elif (self.atomType == 'na  ' and self.resName == 'na'):
           self.atomType = 'sod '
        elif (self.atomType == 'cs  ' and self.resName == 'cs'):
           self.atomType = 'ces '
        elif (self.atomType == 'cl  ' and self.resName == 'cl'):
           self.atomType = 'cla '
        elif (self.atomType == 'ca  ' and self.resName == 'ca'):
           self.atomType = 'cal '
        elif (self.atomType == ' k  ' and self.resName == 'k'):
           self.atomType = 'pot '
