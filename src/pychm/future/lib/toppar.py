"""This module contains the data objects for the chemical information found in
topology and parameter files. The classes are "pure" chemical data holders, and
do not contain (much) implementation specific formatting (ie CHARMM). To get
this chemical data into a useful form, the classes contained in this module
must be used in conjunction with their appropriate I/O class found in
:mod:`pychm.future.io`. The general workflow is as follows: Formatted Text Data -> I/O
Class -> Toppar Classes -> {some useful work} -> I/O Class -> Formatted Text
Data. Detailed, concrete examples can be found in the documentation for
:mod:`pychm.future.io.prm` and :mod:`pychm.future.io.rtf`.

:TODO:
    :class:`Residue` -- Complete residue and patch class logic, for psf
    writing.
    :method:`cull` -- Write logic for "culling" toppar, based upon atomtypes
    present in a :class:`Mol` object
    :exec`CMAP_Exception` -- replace with something less specific
"""

__author__ = ("Frank C. Pickard <frank.pickard@nih.gov>")
__all__ = ['Toppar']

from abc import ABCMeta, abstractmethod, abstractproperty
from copy import deepcopy
import warnings


class CMAP_Exception(Exception):
    """This is a temporary exception, and will be removed/changed in a future
    release. It is only included now for compatibility with the current
    CHARMMing release.
    """
    pass


# Convenience functions ############################
def _myfloat(k):
    """Casts a variable to a float, but doesn't barf exceptions on
    uninitialized (`None`) values.
    """
    try:
        return float(k)
    except TypeError:
        if k is None:
            return None
        else:
            raise

def _myint(k):
    """Casts a variable to an int, but doesn't barf exceptions on
    uninitialized (`None`) values.
    """
    try:
        return int(k)
    except TypeError:
        if k is None:
            return None
        else:
            raise


# Merging functions #################################
def _unique(this):
    """Accepts a :class:`list` and returns a new list containing only the unique
    elements, while preserving the initial order.
    """
    if this is None:    # short circuit for Nonetype
        return None
    if this == []:      # short circuit empty lists
        return None
    this = deepcopy(this)
    if len(this) == len(set(this)): # short circuit for already unique lists
        return this
    tmp = []
    dupes = []
    for data in this:
        if data in tmp:
            dupes.append(data)
        else:
            tmp.append(data)
    if dupes:
        warnings.warn("Duplicate data: %r" % dupes)
    return tmp

def _merge_section(this, that):
    """Function for merging two :class:`list`s containing :class:`PRM`-like
    objects. This function should not be called on :class:`Mass` or
    :class:`CmapPRM` objects, they have their own specialized logic.
    """
    this = _unique(this)
    that = _unique(that)
    if this is None:
        return that
    if that is None:
        return this
    if this == that:
        return this
    #
    for data in that:
        if data not in this:
            this.append(data)
    return this

def _merge_mass(this, that):
    """Function for merging two :class:`list`s containing :class:`Mass` objects.

    First attempts to naively merge the two lists by simply adding them
    (preserving ordering and indexing). Failing that, the lists are combined,
    the ordering is lost, and new :attr:`id`s are assigned.
    """
    this = _unique(this)
    that = _unique(that)
    if this is None:
        return that
    if that is None:
        return this
    if this == that:
        return this
    # try merging naively, then check for overlap atom and id
    tmp = this + that
    tmpid = [ x.id for x in tmp ]
    if len(this) + len(that) == len(set(tmp)) and len(this) + len(that) == len(set(tmpid)):
        # naive merge works!
        return sorted(tmp, key=lambda x: x.id)
    # naive merge fails
    tmp = sorted(list(set(tmp)))
    for i in range(len(tmp)):
        tmp[i].id = i+1
    warnings.warn("Original mass id's could not be retained!")
    return tmp

def _merge_cmap(this, that):
    """Function for merging two :class:`CmapPRM` objects. It is currently
    stupid.

    If two different CMAP objects are given, a :exec:`CMAP_Exception` is raised.
    """
    this = _unique(this)
    that = _unique(that)
    if this is None:
        return that
    if that is None:
        return this
    if this == that:
        return this
    raise CMAP_Exception("CMAPs are different between toppar")

def _merge_command(this, that):
    """Function for merging two command_strings. It is currently stupid.

    If two different strings are given, the first one is merely returned, and a
    warning is thrown.
    """
    if this is None:
        return that
    if that is None:
        return this
    if this == that:
        return this
    warnings.warn("Conflicting command strings:\n\n%s\n\n%s\n\nusing first string" % (this,that))
    return this


###############################################################################
# Class Definitions ###########################################################
###############################################################################
class Toppar(object):
    """The class containing chemical data related to molecular topology and
    bonding parameters.
    """
    data_sections = ['bond', 'angle', 'dihedral', 'improper', 'cmap',
                    'nonbond', 'nbfix', 'hbond', 'mass', 'residue', 'patch']
    command_sections = ['declare', 'autogen', 'default']

    def __init__(self):
        self.sections = self.data_sections + self.command_sections
        # data ################################################################
        self.bond = None
        self.angle = None
        self.dihedral = None
        self.improper = None
        self.cmap = None
        self.nonbond = None
        self.nbfix = None
        self.hbond = None
        self.mass = None
        self.residue = None
        self.patch = None
        # commands ############################################################
        #### a dictionary where we stuff unformatted text data, this aspect of
        #### the implementation may change as we add non-charmm formats
        self._init_commands() # sets attr:`self.commands`

    def merge(self, other):
        """Creates a new Toppar object, where rtf/prm objects from `self` are
        given priority over rtf/prm objects from `other`.
        """
        tmp = Toppar()
        # merge data attributes
        for section in self.sections:
            try:
                if section == 'mass':
                    merged = _merge_mass(self.mass, other.mass)
                elif section == 'cmap':
                    merged = _merge_cmap(self.cmap, other.cmap)
                else:
                    merged = _merge_section(getattr(self, section),
                                            getattr(other, section))
                setattr(tmp, section, merged)
            except AttributeError:
                pass
        # merge command attributes
        for section in self.sections:
            merged = _merge_command(self.commands[section], other.commands[section])
            if merged is not None:
                tmp.commands[section] = merged
        return tmp

    def unique(self):
        """Checks each :attr:`data_section` in the current :class:`Toppar` object
        for redundencies, and removes them.
        """
        for section in self.sections:
            try:
                setattr(self, section, _unique(getattr(self, section)))
            except AttributeError:
                pass

    def __add__(self, other):
        return self.merge(other)

    def _init_commands(self):
        """Initialize all possible command values to `None`."""
        tmp = dict()
        for section in self.sections:
            tmp[section] = None
        self.commands = tmp

class PRM(ABCMeta):
    """A dummy :class:`abc.ABCMeta`-class, serves as a registrar and API
    definition.
    """
    # API Definition ################################################
    @abstractproperty
    def _sortkey(self):
        raise NotImplementedError

    @abstractmethod
    def _validate(self):
        raise NotImplementedError


class BasePRM(object):
    """Foundation class for PRM objects. Can not be instantiated because it
    lacks `abstract` definitions from class:`PRM`.

    All methods defined here depend upon :attr:`_sortkey`, a property defining
    how PRM objects are hashed. This hash is in turn used for sorting and
    determining uniqueness.
    """
    __metaclass__ = PRM

    # Magic Methods #################################################
    def __hash__(self):
        return hash(self._sortkey)

    def __repr__(self):
        return '%s%r' % (self.__class__.__name__, self._sortkey)

    #### Comparison Methods #########################################
    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self._sortkey == other._sortkey

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self._sortkey < other._sortkey

    def __le__(self, other):
        return self._sortkey <= other._sortkey

    def __gt__(self, other):
        return self._sortkey > other._sortkey

    def __ge__(self, other):
        return self._sortkey >= other._sortkey


class BondPRM(BasePRM):
    """Object representing bonding parameter constants.

    The first two arguments are strings representing atomtypes. During
    initialization these strings are normalized, such that 'atom0' is always
    the "lesser" string, allowing BondPRM objects to be hashed correctly.
    Changing the 'atom' attributes after initialization is *not* recommended
    because uniqueness is not guaranteed. For example, Bond(AB) == Bond(BA) if
    set at initialization, but not equal if set afterwards. Changing the spring
    constant, 'k' and the equilibrium bond length 'eq' is fine at any time.
    """
    __slots__ = ['atom0', 'atom1', 'k', 'eq']

    def __init__(self, atom0=None, atom1=None, k=None, eq=None):
        self.atom0, self.atom1 = sorted((atom0, atom1))
        self.k = _myfloat(k)
        self.eq = _myfloat(eq)

    @property
    def _sortkey(self):
        return (self.atom0, self.atom1)

    def _validate(self):
        for attrname in self.__slots__:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))


class AnglePRM(BasePRM):
    """Object representing angle bending parameter constants.

    The first three arguments are strings representing atomtypes. During
    initialization these strings are normalized, such that 'atom0' is always
    the "lesser" string, allowing AnglePRM objects to be hashed correctly.
    Changing the 'atom' attributes after initialization is *not* recommended
    because uniqueness is not guaranteed. For example, Angle(ABC) == Angle(CBA)
    if set at initialization, but not equal if set afterwards.  Changing the
    spring constant, 'k' and the equilibrium bond length 'eq' is fine at any
    time. The 1-3 interaction parameters may be optionally set.
    """
    __slots__ = ['atom0', 'atom1', 'atom2', 'k', 'eq', 'k13', 'eq13']

    def __init__(self, atom0=None, atom1=None, atom2=None, k=None, eq=None, k13=None, eq13=None):
        self.atom0, self.atom2 = sorted((atom0, atom2))
        self.atom1 = atom1
        self.k = _myfloat(k)
        self.eq = _myfloat(eq)
        self.k13 = _myfloat(k13)
        self.eq13 = _myfloat(eq13)

    @property
    def _sortkey(self):
        return (self.atom0, self.atom1, self.atom2)

    def _validate(self):
        for attrname in self.__slots__[:5]:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))
        if self.k13 is None and self.eq13 is not None:
            warnings.warn("%r has uninitialized attr: k13" % self)
        if self.eq13 is None and self.k13 is not None:
            warnings.warn("%r has uninitialized attr: eq13" % self)


class DihedralPRM(BasePRM):
    """Object representing dihedral angle parameter constants.

    The first four arguments are strings representing atomtypes. During
    initialization these strings are normalized, such that 'atom0' is always
    the "lesser" string, allowing DihedralPRM objects to be hashed correctly.
    Changing the 'atom' attributes after initialization is *not* recommended
    because uniqueness is not guaranteed. For example, Dihedral(ABCD) ==
    Dihedral(DCBA) if set at initialization, but not equal if set afterwards.
    Furthermore, the multiplicity of the dihedral, 'mult' is also used when
    determining parameter uniqueness. Changing the spring constant, 'k' and the
    equilibrium bond length 'eq' is fine at any time.
    """
    __slots__ = ['atom0', 'atom1', 'atom2', 'atom3' 'k', 'mult', 'eq']

    def __init__(self, atom0=None, atom1=None, atom2=None, atom3=None, k=None, mult=None, eq=None):
        self.atom0, self.atom3 = sorted((atom0, atom3))
        if self.atom0 == atom0:
            self.atom1, self.atom2 = (atom1, atom2)
        else:
            self.atom1, self.atom2 = (atom2, atom1)
        self.k = _myfloat(k)
        self.mult = _myint(mult)
        self.eq = _myfloat(eq)

    @property
    def _sortkey(self):
        return (self.atom0, self.atom1, self.atom2, self.atom3, self.mult)

    def _validate(self):
        for attrname in self.__slots__:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))


class ImproperPRM(DihedralPRM):
    """Object representing improper angle parameter constants.

    The first four arguments are strings representing atomtypes. During
    initialization these strings are normalized, such that 'atom0' is always
    the "lesser" string, allowing DihedralPRM objects to be hashed correctly.
    Changing the 'atom' attributes after initialization is *not* recommended
    because uniqueness is not guaranteed. For example, Improper(ABCD) ==
    Improper(DCBA) if set at initialization, but not equal if set afterwards.
    Furthermore, the multiplicity of the Improper, 'mult' is also used when
    determining parameter uniqueness. Changing the spring constant, 'k' and the
    equilibrium bond length 'eq' is fine at any time.
    """
    __slots__ = ['atom0', 'atom1', 'atom2', 'atom3' 'k', 'mult', 'eq']


class CmapPRM(BasePRM):
    """Object representing CHARMM CMAP parameters. It is currently a black box.
    """
    __slots__ = ['text']

    def __init__(self, arg=None):
        self.text = arg

    @property
    def _sortkey(self):
        return ('CMAP_BLOCK',)

    def _validate(self):
        pass


class NonbondPRM(BasePRM):
    """Object representing unary nonbonding VDW parameter constants.

    The first argument is a string representing the atomtype. The next three
    arguments, 'ig', 'k' and 'eq' are ignored, the LJ 6-12 well depth and the
    equilibrium interaction distance upon two, respectively.  The last three
    arguments are optional, and only used for calculating 1-4 nonbonded
    interactions.
    """
    __slots__ = ['atom', 'ig', 'k', 'eq', 'ig14', 'k14', 'eq14']

    def __init__(self, atom=None, ig=None, k=None, eq=None, ig14=None, k14=None, eq14=None):
        self.atom = atom
        self.ig = _myfloat(ig)     # ignored
        self.k = _myfloat(k)       # epsilon
        self.eq = _myfloat(eq)     # r_min/2
        self.ig14 = _myfloat(ig14)
        self.k14 = _myfloat(k14)
        self.eq14 = _myfloat(eq14)

    @property
    def _sortkey(self):
        return (self.atom,)

    def _validate(self):
        for attrname in self.__slots__[:4]:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))
        if self.ig14 is None and (self.k14 is not None or self.eq14 is not None):
            warnings.warn("%r has uninitialized attr: ig14" % self)
        if self.k14 is None and (self.ig14 is not None or self.eq14 is not None):
            warnings.warn("%r has uninitialized attr: k14" % self)
        if self.eq14 is None and (self.ig14 is not None or self.k14 is not None):
            warnings.warn("%r has uninitialized attr: eq14" % self)


class NBFixPRM(BasePRM):
    """Object representing binary nonbonding VDW parameter constant overrides.
    These parameters explicitly override pairwise parameters procedurally
    generated bu Nonbond parameters.

    The first two arguments are strings representing atomtypes. During
    initialization these strings are normalized, such that 'atom0' is always
    the "lesser" string, allowing NBFixPRM objects to be hashed correctly.
    Changing the 'atom' attributes after initialization is *not* recommended
    because uniqueness is not guaranteed. For example, NBFix(AB) == NBFix(BA)
    if set at initialization, but not equal if set afterwards. The arguments
    'k' and 'eq' are the LJ 6-12 well depth and equilibrium distances,
    respectively. The last two arguments are optional, and are used for 1-4
    nonbonded interactions.
    """
    __slots__ = ['atom0', 'atom1', 'k', 'eq', 'k14', 'eq14']

    def __init__(self, atom0=None, atom1=None, k=None, eq=None, k14=None, eq14=None):
        self.atom0, self.atom1 = sorted((atom0, atom1))
        self.k = _myfloat(k)
        self.eq = _myfloat(eq)
        self.k14 = _myfloat(k14)
        self.eq14 = _myfloat(eq14)

    @property
    def _sortkey(self):
        return (self.atom0, self.atom1)

    def _validate(self):
        for attrname in self.__slots__[:4]:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))
        if self.ig14 is None and (self.k14 is not None or self.eq14 is not None):
            warnings.warn("%r has uninitialized attr: ig14" % self)
        if self.k14 is None and (self.ig14 is not None or self.eq14 is not None):
            warnings.warn("%r has uninitialized attr: k14" % self)
        if self.eq14 is None and (self.ig14 is not None or self.k14 is not None):
            warnings.warn("%r has uninitialized attr: eq14" % self)


class HBondPRM(BasePRM):
    """Object representing hydrogen bonding interactions, their use in CHARMM
    is now depricated.

    The first two arguments are strings representing atomtypes. During
    initialization these strings are normalized, such that 'atom0' is always
    the "lesser" string, allowing NBFixPRM objects to be hashed correctly.
    Changing the 'atom' attributes after initialization is *not* recommended
    because uniqueness is not guaranteed. For example, HBond(AB) == HBond(BA)
    if set at initialization, but not equal if set afterwards. The arguments
    'k' and 'eq' are the well depth and equilibrium distances, respectively.
    """
    __slots__ = ['atom0', 'atom1', 'k', 'eq']

    def __init__(self, atom0=None, atom1=None, k=None, eq=None):
        self.atom0, self.atom1 = sorted((atom0, atom1))
        self.k = _myfloat(k)
        self.eq = _myfloat(eq)

    @property
    def _sortkey(self):
        return (self.atom0, self.atom1)

    def _validate(self):
        for attrname in self.__slots__[:4]:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))


class Mass(BasePRM):
    """Object representing particle mass data.

    The first argument is a string representing the atomtype. The second
    argument (this might be changed to the last arg, and be made optional) is
    the unique integer ID used by CHARMM. The third argument is a float
    representing the mass in AMU. The fourth argument is an optional string
    representing the elemental symbol.
    """
    __slots__ = ['atom', 'id', 'mass', 'element']

    def __init__(self, id=None, atom=None, mass=None, element=None):
        self.id = int(id)
        self.atom = atom
        self.mass = _myfloat(mass)
        self.element = element

    @property
    def _sortkey(self):
        return (self.atom,)

    def _validate(self):
        for attrname in self.__slots__[:3]:
            if getattr(self, attrname) is None:
                warnings.warn("%r has uninitialized attr: %s" % (self, attrname))


class Residue(BasePRM):
    """Object representing residue information. Mostly a black box.

    The first two arguments are the residue name and its total charge. The remaining
    information is stored as a black box in the third argument. This will likely change
    in the near future and a fully featured residue object will be implemented.
    """
    def __init__(self, name=None, charge=None, body=None):
        self.name = name
        self.charge = _myfloat(charge)
        self.body = body

    @property
    def _sortkey(self):
        return "a%r" % self.name

    def _validate(self):
        pass

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, self.name)


class Patch(Residue):
    """Object representing patch residue information. Mostly a black box.

    The first two arguments are the residue name and its total charge. The
    remaining information is stored as a black box in the third argument. This
    will likely change in the near future and a fully featured residue object
    will be implemented.
    """
    @property
    def _sortkey(self):
        return "b%r" % self.name
