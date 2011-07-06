"""
:Author: fcp
:Date: 10/22/2010
"""


from pychm.const.bio import atomMass, backbone, charmm2pdbAtomNames, good, \
        pro, nuc
from pychm.tools import Property, lowerKeys
from pychm.lib.metaatom import AtomError
from pychm.lib.baseatom import BaseAtom


class Atom(BaseAtom):
    """
    :Note:  This class is derived from :mod:`pychm.lib.baseatom`,
        please familiarize yourself with that documentation before
        proceeding with this article.


    The standard CHARMMing implementation of the all-atom resolution
    atom-like object.

    All **Core Data** defined by parent classes is also implicitly
    included.

    **Core Data:**
        | ``element``

    Both of these private methods are used exclusively by
    :mod:`pychm.scripts.parse`.

    **Private Methods:**
        | ``_compliance_resName``
        | ``_compliance_atomType``

    :TODO:
        | # Input and output *AMBER* style atomic data.
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
        self._element = 'auto'
        self._segType = 'auto'
        super(Atom, self).__init__(text, **kwargs)

####################
# Parse Definition #
####################

    def parse(self, inFormat):
        """
        Parses ``self._text`` into the following instance variables:
            | ``atomNum``
            | ``atomType``
            | ``resName``
            | ``chainid``
            | ``resid``
            | ``cart``
            | ``weight``
            | ``bFactor``
            | ``element``
            | ``resIndex``

        Initial values for addressing are stored in:
            | ``chainid0``
            | ``segType0``
            | ``resid0``
            | ``atomNum0``

        Recognized *inFormat* values are:
            | ``"pdborg"``
            | ``"charmm"``
            | ``"shortcard"``
            | ``"longcard"``
            | ``"amber"``   **TODO**
        """
        if inFormat == 'pdborg':
            self.atomNum = self._text[6:11]
            self.resName = self._text[16:20]
            self.atomType = self._text[12:16]
            self.chainid = self._text[21]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.weight = self._text[55:60]
            self.bFactor = self._text[61:66]
            self.resIndex = self.resid
        elif inFormat == 'charmm':
            self.atomNum = self._text[6:11]
            self.resName = self._text[17:21]
            self.atomType = self._text[12:16]
            self.chainid = self._text[72:76]
            self.resid = self._text[22:26]
            self.cart = (self._text[30:38], self._text[38:46], self._text[46:54])
            self.weight = self._text[55:60]
            self.bFactor = self._text[61:66]
            self.resIndex = self.resid
        elif inFormat in ['crd', 'cor', 'card', 'short', 'shortcard']:
            self.atomNum = self._text[0:5]
            self.resid = self._text[5:10]
            # null [10:11]
            self.resName = self._text[11:15]
            # null [15:16]
            self.atomType = self._text[16:20]
            self.cart = (self._text[20:30], self._text[30:40], self._text[40:50])
            # null [50:51]
            self.chainid = self._text[51:55]
            # null [55:56]
            self.resid = self._text[56:60]
            self.weight = self._text[60:70]
            self.bFactor = 1.
        elif inFormat in ['xcrd', 'xcor', 'xcard', 'long', 'longcard']:
            self.atomNum = self._text[0:10]
            self.resIndex = self._text[10:20]
            # null [20:22]
            self.resName = self._text[22:30]
            # null [30:32]
            self.atomType = self._text[32:40]
            self.cart = (self._text[40:60], self._text[60:80], self._text[80:100])
            # null [100:102]
            self.chainid = self._text[102:110]
            # null [110:112]
            self.resid = self._text[112:120]
            self.weight = self._text[120:140]
            self.bFactor = 1.
        elif inFormat == 'amber':
            raise NotImplementedError
        else:
            raise AtomError('parse: unknown `inFormat`: %s' % inFormat)
        # fix white space padding on atomTypes
        self.fix_atomType()
        # Save initial properties
        self._chainid0 = self.chainid
        self._segType0 = self.segType
        self._resid0 = self.resid
        self._atomNum0 = self.atomNum
        self._atomType0 = self.atomType
        self._resName0 = self.resName

##############
# Properties #
##############

    @Property
    def element():
        doc =\
        """
        A string which defines the cannonical 2 character element
        abbreviation.  If set to *'auto'*, it polls the ``atomType``
        property for this information.  Otherwise it is set manually.
        """
        def fget(self):
            if self._element == 'auto':
                return self._atomType[:2].strip()
            else:
                return self._element
        def fset(self, value):
            self._element = str(value).strip().lower()
        return locals()

    @Property
    def mass():
        doc =\
        """
        Polls the ``element`` property, and returns the mass in AMU of
        the :class:`Atom`.  Read only.
        """
        def fget(self):
            return atomMass[self.element]
        return locals()

    @Property
    def segType():
        doc =\
        """
        A string which defines the biochemical nature of the current
        segment.  Valid ``segTypes`` for the :class:`Atom` are:
        *"pro" "nuc" "rna" "dna" "good" "bad"*.  If set to *'auto'*,
        it polls the ``resName`` for this information.  Otherwise it
        is set manually.
        """
        def fget(self):
            if self._segType == 'auto':
                if self.is_pro():
                    return 'pro'
                if self.is_nuc():
                    return 'nuc'
                if self.is_good():
                    return 'good'
                return 'bad'
            else:
                return self._segType
            return self._segType
        def fset(self, value):
            self._segType = str(value).strip().lower()
        return locals()

##################
# Public Methods #
##################

    def is_backbone(self):
        """
        Determines if the atom is part of the backbone structure using
        the ``atomType`` property.  Returns a :class:`bool`.
        """
        return self.atomType in backbone

    def is_good(self):
        """
        Determines if the atom is a "good" hetero atom using the
        ``resName`` property.  Returns a :class:`bool`.
        """
        return self.resName in good

    def is_pro(self):
        """
        Determines if the atom is part of a protein structure using the
        ``resName`` property.  Returns a :class:`bool`.
        """
        return self.resName[-3:] in pro

    def is_nuc(self):
        """
        Determines if the atom is part of a nucleic acid structure
        using the ``resName`` property.  Returns a :class:`bool`.
        """
        return self.resName[-3:] in nuc

    def fix_atomType(self):
        """
        Fixes the white space padding on ``atomType`` strings to
        conform to the 4 character standard.

        Where the first 2 characters represent the element, and the
        last 2 characters the chemical environment.  It accounts for
        the common exception to this rule where 1 character elements
        may use 3 characters for their chemical environment.
        """
        if self.segType == 'good':
            tmp = self._atomType.strip()
            try:
                self._atomType = charmm2pdbAtomNames[tmp]
            except KeyError:
                if len(tmp) == 4:
                    self._atomType = ' %s' % tmp
                elif len(tmp) in [2, 3]:
                    tmp = '%s  ' % tmp
                    self._atomType = tmp[:4]
                else:
                    tmp = ' %s  ' % tmp
                    self._atomType = tmp[:4]
        elif self.segType == 'pro':
            tmp = self._atomType.strip()
            if len(tmp) == 4:
                self._atomType = ' %s' % tmp
            else:
                tmp = ' %s  ' % tmp
                self._atomType = tmp[:4]

    def Print(self, **kwargs):
        """
        Returns a string representation of the :class:`Atom`. Formatting
        defaults to `self.__class__._autoInFormat`.

        **kwargs:**
            | ``outformat``     ["pdborg","charmm","debug","xdebug","crd","xcrd"]
            | ``old_chainid``   [False,True]
            | ``old_segtype``   [False,True]
            | ``old_resid``     [False,True]
            | ``old_atomnum``   [False,True]

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
        # Unretrofit whitespace padding for long (hydrogen) atomType...
        if len(self.atomType) > 4:
            tmpAtomType = self.atomType[-4:]
        else:
            tmpAtomType = self.atomType
        # Set formatting
        if outFormat == 'pdborg':
            taco = '%-6s%5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s' % \
                (self.tag, atomNum, tmpAtomType, self.resName, chainid, resid,
                x, y, z, self.weight, self.bFactor, self.element)
        elif outFormat == 'charmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s' % \
                (self.tag, atomNum, tmpAtomType, self.resName, resid,
                x, y, z, self.weight, self.bFactor, chainid)
        elif outFormat == 'debug':
            taco = '%-6s%5i %4s %4s %1s%4i    %8.3f%8.3f%8.3f' % \
                (segType, atomNum, tmpAtomType, self.resName, chainid, resid,
                x, y, z)
        elif outFormat == 'xdebug':
            taco = '%15s > %15s: %5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f' % \
                (self.addr0, self.addr, atomNum, tmpAtomType, self.resName,
                chainid, resid, x, y, z, self.weight, self.bFactor)
        elif outFormat == 'repr':
            taco = '%15s :: %4s %4s    %8.3f %8.3f %8.3f' % \
                (self.addr, tmpAtomType, self.resName, x, y, z)
        elif outFormat in ['crd', 'cor', 'card', 'short', 'shortcard']:
            taco = '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f' % \
                (atomNum, self.resIndex, self.resName, tmpAtomType,
                x, y, z, chainid, resid, self.weight)
        elif outFormat in ['xcrd', 'xcor', 'xcard', 'long', 'longcard']:
            taco = '%10i%10i  %-8s  %-8s    %20.10f%20.10f%20.10f  %-8s  %-8i    %20.10f' % \
                (atomNum, self.resIndex, self.resName, tmpAtomType,
                x, y, z, chainid + segType, resid, self.weight)
        else:
            raise AtomError('Print: unknown inFormat %s' % inFormat)
        return taco.upper()

###################
# Private Methods #
###################

    def _compliance_resName(self):
        """
        Re-label ``resName`` string to be CHARMM compliant.
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
        Re-label ``atomType`` string to be CHARMM compliant.
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
