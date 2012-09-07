"""
:Author: fcp
:Date: 10/22/2010
"""


from pychm.const import alphanum
from pychm.lib.metaatom import AtomError, MetaAtom
from pychm.tools import Property


class BaseAtom(MetaAtom):
    """
    :Note:  This class is derived from :mod:`pychm.lib.metaatom`,
        please familiarize yourself with that documentation before
        proceeding with this article.


    This is the base class for all "CHARMM-like" atoms.  It includes
    property definitions and methods useful for classical mechanics
    based molecular models.  However, like :mod:`pychm.lib.metaatom`,
    :mod:`baseatom` is still abstract, and thus not useful without
    deriving and defining the methods marked **STUB** herein.

    All **Core Data** defined by parent classes is also implicitly
    included.

    **Core Data:**
        | ``atomNum``
        | ``atomType``
        | ``bFactor``
        | ``chainid``
        | ``resid``
        | ``resIndex``
        | ``resName``
        | ``segType``
        | ``weight``

    As some **STUBS** are defined in this class, and others are added,
    the following list explicitly states all **STUBS** that must be
    defined by child classes.

    **STUBS:**
        | ``_autoInFormat``
        | ``_sortSegType``
        | ``_tagMap``
        | ``parse``
        | ``Print``

    This private method has been documented as a **Special Method** because
    it defines the scoring system used by the rich comparison operators.

    The weighting for this scoring functions is:
        ``chainid`` > ``segType`` > ``resid`` > ``atomNum``

    Where ``chainid`` is sorted from `a` to `9`, ``segType`` is sorted
    according to the :class:`dict` defined by ``self.__class__._sortSegType``
    and ``resid`` and ``atomNum`` are straightforward.

    **Special Methods:**
        | ``_sort``
    """


    _properties = dict(MetaAtom._properties)
    _properties.update(
        {
        'atomNum': 0,
        'atomType': 'unkt',
        'bFactor': 0,
        'chainid': '?',
        'resid': 0,
        'resIndex': 0,
        'resName': 'unkr',
        'segType': 'bad',
        'weight': 0
        }
    )
    """
    A dictionary which tells the class which properties (keys) it is
    associated with.  The corresponding values serve as default values
    for said properties.
    """


    _
sortSegType = None
    """
    This is a dictionary which defines the recognized ``segTypes`` (keys)
    for the ``_sort`` scoring method, and their weighting (values).

    **STUB**
    """

    _tagMap = None
    """
    Map ``segType`` -> ``tag``

    This is a dictionary which defines the recognized `segTypes` (keys)
    for the `tag` property.

    **STUB**
    """

    def __init__(self, text=None, **kwargs):
        super(BaseAtom, self).__init__(text, **kwargs)

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        The ``addr`` property provides a human readable unique string
        representation for each `BaseAtom` instance.

        The default is: "``chainid.segType.resid.atomNum``"
        """
        def fget(self):
            return '%s.%s.%d.%d' % (self.chainid, self.segType, self.resid,
                                    self.atomNum)
        return locals()

    @Property
    def atomNum():
        doc =\
        """
        An atomic index which spans all atoms in a single model or
        :class:`Mol` object.
        """
        def fget(self):
            return self._atomNum
        def fset(self, value):
            value = int(value)
            if value < 0:
                if self._autoFix:
                    value = 0
                else:
                    raise AtomError('atomNum %s: %d is less than 0' %
                                    (self.addr0, value))
            elif value > 10000:
                if self._autoFix:
                    value = 10000
                else:
                    raise AtomError('atomNum %s: %d is less than 0' %
                                    (self.addr0, value))
            self._atomNum = value
        return locals()

    @Property
    def atomNum0():
        doc =\
        """
        The value of the ``atomNum`` property at instantization, read only.
        """
        def fget(self):
            return self._atomNum0
        return locals()

    @Property
    def atomType():
        doc =\
        """
        The four character string which defines the atomic element,
        and neighboring chemical environment.  Typically the first
        two characters define the element, and the final two define
        the chemical environment.  However many sources do not adhere
        to this rule.

        **Note:** ``atomType`` should always be 4 characters long,
        and whitespace should never be stripped.
        Remember kids, ``' CA ' != 'CA  '``.

        **Note:** The one exception to the above rule, that is now
        recognized by pychm, is for elements which have a single
        character abbreviation, and might have 3 characters denoting
        their type, for example ``HD22``.  Other single character
        elements should also work, however this has not been tested.
        """
        def fget(self):
            return self._atomType
        def fset(self, value):
            value = str(value).lower()
            if len(value) > 5:
                if self._autoFix:
                    value = value[:5]
                else:
                    raise AtomError('atomType %s: %s is longer than 5 characters' %
                                    (self.addr0, value))
            self._atomType = value
        return locals()

    @Property
    def atomType0():
        doc =\
        """
        The value of the ``atomType`` property at instantization, read only.
        """
        def fget(self):
            return self._atomType0
        return locals()

    @Property
    def bFactor():
        doc =\
        """
        A float in the range [0,100].  It is a proxy for the quality of the
        structural model.  Lower values are better.
        """
        def fget(self):
            return self._bFactor
        def fset(self, value):
            try:
                value = float(value)
            except ValueError:
                if self._autoFix:
                    value = 100.0
                else:
                    raise AtomError('bFactor %s: bFactor value is missing' %
                                    self.addr0)
            if value > 100.:
                if self._autoFix:
                    value = 100.
                else:
                    raise AtomError('bFactor %s: %5.2f is greater than 100.0' %
                                    (self.addr0, value))
            elif value < 0.:
                if self._autoFix:
                    value = 0.
                else:
                    raise AtomError('bFactor %s: %5.2f is less than 0.0' %
                                    (self.addr0, value))
            self._bFactor = value
        return locals()

    @Property
    def chainid():
        doc =\
        """
        A single character string indicating the domain of the structure
        where the atom is located.  It is in the range of [`'a'`,`9`].
        """
        def fget(self):
            return self._chainid
        def fset(self, value):
            value = str(value).strip().lower()
            if len(value) > 1:
                if self._autoFix:
                    value = value[0]
                else:
                    raise AtomError('chainid %s: %s is longer than 1 character' %
                                    (self.addr0, value))
            self._chainid = value
        return locals()

    @Property
    def chainid0():
        doc =\
        """
        The value of the `chainid` property at instantization, read only.
        """
        def fget(self):
            return self._chainid0
        return locals()

    @Property
    def resid():
        doc =\
        """
        An integer in the range of [-1000,10000].  It denotes the residue
        by number, and to facilitate CHARMM jobs, it is typically reindexed
        to start at one for each segment.  To retrieve the cannonical
        resid refer to ``resid0``.
        """
        def fget(self):
            return self._resid
        def fset(self, value):
            value = int(value)
            if value < -1000:
                if self._autoFix:
                    value = -1000
                else:
                    raise AtomError('resid %s: %d is less than -1000' %
                                    (self.addr0, value))
            elif value > 10000:
                if self._autoFix:
                    value = 10000
                else:
                    raise AtomError('resid %s: %d is greater than 10000' %
                                    (self.addr0, value))
            self._resid = value
        return locals()

    @Property
    def resid0():
        doc =\
        """
        The value of the ``resid`` property at instantization, read only.
        """
        def fget(self):
            return self._resid0
        return locals()

    @Property
    def resIndex():
        doc =\
        """
        An integer in the range of [-1000,10000].  It denotes the residue
        by number and is indexed from 1.  Unlike the ``resid`` it indexes
        the residues present in the entire molecule, not just in the
        current segment.
        """
        def fget(self):
            return self._resIndex
        def fset(self, value):
            value = int(value)
            if value < -1000:
                if self._autoFix:
                    value = -1000
                else:
                    raise AtomError('resIndex %s: %d is less than -1000' %
                                    (self.addr0, value))
            elif value > 10000:
                if self._autoFix:
                    value = 10000
                else:
                    raise AtomError('resIndex %s: %d is greater than 10000' %
                                    (self.addr0, value))
            self._resIndex = value
        return locals()

    @Property
    def resName():
        doc =\
        """
        A 4 character string denoting the name of the residue.  As per
        the PDB specification the first character is reserved for multi-
        model sections of PDB files, and the last 3 characters are used
        for the actual residue name itself.  This convention is loosely
        adhered to.
        """
        def fget(self):
            return self._resName
        def fset(self, value):
            value = str(value).strip().lower()
            if len(value) > 4:
                if self._autoFix:
                    value = value[:4]
                else:
                   raise AtomError('resName %s: %s is longer than 4 characters' %
                                (self.addr0, value))
            self._resName = value
        return locals()

    @Property
    def resName0():
        doc =\
        """
        The value of the `resName` property at instantization, read only.
        """
        def fget(self):
            return self._resName0
        return locals()

    @Property
    def segid():
        doc =\
        """
        A unique string generated from the segment's parent chainid and
        the segment's type.  Read only.
        """
        def fget(self):
            return '%s-%s' % (self.chainid, self.segType)
        return locals()

    @Property
    def segType():
        doc =\
        """
        A string which defines the biochemical nature of the current
        segment.  Valid ``segTypes`` should be registered as keys in
        ``self.__class__._sortSegType``.
        """
        def fget(self):
            return self._segType
        def fset(self, value):
            self._segType = str(value).strip().lower()
        return locals()

    @Property
    def segType0():
        doc =\
        """
        The value of the `segType` property at instantization, read only.
        """
        def fget(self):
            return self._segType0
        return locals()

    @Property
    def tag():
        doc =\
        """
        The first field which appears in a PDB line, usually 'ATOM' or
        'HETATM'.  A mapping between ``segType`` and ``tag`` should be
        defined by a dictionary ``self.__class__._tagMap``.  Read only.
        """
        def fget(self):
            return self.__class__._tagMap[self.segType]
        return locals()

    @Property
    def weight():
        doc =\
        """
        In multi-model sections of .pdb files the ``weight`` property
        defines the experimental confidence of the atomic coordinates.
        When stripping out multi-model sections of a .pdb, the ``weight``
        property is used to determine the most likely atomic coordinates
        and removes the others.
        """
        def fget(self):
            return self._weight
        def fset(self, value):
            try:
                value = float(value)
            except ValueError:
                if self._autoFix:
                    value = 0.
                else:
                    raise AtomError('weight %s: weight value is missing' %
                                    self.addr0)
            if value > 1.:
                if self._autoFix:
                    value = 1.
                else:
                    raise AtomError('weight %s: %5.2f is greater than 1.0' %
                                    (self.addr0, value))
            elif value < 0.:
                if self._autoFix:
                    value = 0.
                else:
                    raise AtomError('weight %s: %5.2f is less than 0.0' %
                                    (self.addr0, value))
            self._weight = value
        return locals()

###################
# Private Methods #
###################

    def _sort(self):
        """
        Scoring method that determines sorting order, for rich
        comparisons and containerClass.sort() methods.

        The weighting for this scoring functions is as follows:
            ``chainid`` > ``segType`` > ``resid`` > ``atomNum``
        """
        sortChainid = dict(((char,i) for i,char in enumerate(alphanum)))
        return sortChainid[self.chainid] * 1e10 + \
                self.__class__._sortSegType[self.segType] * 1e9 + \
                self.resid * 1e5 + self.atomNum
