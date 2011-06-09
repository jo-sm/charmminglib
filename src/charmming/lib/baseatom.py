"""
DOCME
"""
# fcp
# 10/22/2010


from charmming.const import alphanum
from charmming.lib.metaatom import AtomError, MetaAtom
from charmming.tools import Property


class BaseAtom(MetaAtom):
    """
    This is the base class for all "CHARMM-like" atoms.  It includes
    property definitions and methods useful for classical mechanics
    based molecular models.

    Class Attributes
        `_autoInFormat`         STUB
        `_segTypeSortScore`     STUB
        `_tagMap`               STUB
    Parse Definition
        `parse`                 STUB
    Properties
        `addr`
        `atomNum`
        `atomNum0`
        `atomType`
        `bFactor`
        `chainid`
        `chainid0`
        `resid`
        `resid0`
        `resIndex`
        `resName`
        `segid`
        `segType`
        `segType0`
        `tag`
        `weight`
    Public Methods
        `Print`                 STUB
    Private Methods
        `_sort`
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


    _sortSegType = None
    """
    This is a dictionary which defines the recognized `segTypes` (keys)
    for the `_sort` scoring method, and their weighting (values).

    STUB
    """

    _tagMap = None
    """
    Map `segType` -> `tag`

    This is a dictionary which defines the recognized `segTypes` (keys)
    for the `tag` property.

    STUB
    """

    def __init__(self, text=None, **kwargs):
        """
        DOCME
        """
        super(BaseAtom, self).__init__(text, **kwargs)

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        The `addr` property provides a human readable unique string
        representation for each `BaseAtom` instance.

        The default is: "`chainid`.`segType`.`resid`.`atomNum`"
        """
        def fget(self):
            return '%s.%s.%d.%d' % (self.chainid, self.segType, self.resid,
                                    self.atomNum)
        return locals()

    @Property
    def atomNum():
        doc =\
        """
        DOCME
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
        The value of the `atomNum` property at instantization, read only.
        """
        def fget(self):
            return self._atomNum0
        return locals()

    @Property
    def atomType():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._atomType
        def fset(self, value):
            value = str(value).lower()
            if len(value) > 4:
                if self._autoFix:
                    value = value[:4]
                else:
                    raise AtomError('atomType %s: %s is longer than 4 characters' %
                                    (self.addr0, value))
            self._atomType = value
        return locals()

    @Property
    def bFactor():
        doc =\
        """
        DOCME
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
        DOCME
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
        DOCME
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
        The value of the `resid` property at instantization, read only.
        """
        def fget(self):
            return self._resid0
        return locals()

    @Property
    def resIndex():
        doc =\
        """
        DOCME
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
        DOCME
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
    def segid():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return '%s-%s' % (self.chainid, self.segType)
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
        DOCME
        """
        def fget(self):
            return self.__class__._tagMap[self.segType]
        return locals()

    @Property
    def weight():
        doc =\
        """
        DOCME
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
            chainid > segType > resid > atomNum
        """
        sortChainid = dict(((char,i) for i,char in enumerate(alphanum)))
        return sortChainid[self.chainid] * 1e10 + \
                self.__class__._sortSegType[self.segType] * 1e9 + \
                self.resid * 1e5 + self.atomNum
