

__all__ = ['Toppar']

from abc import ABCMeta, abstractmethod, abstractproperty
from copy import deepcopy
import warnings


# Convenience functions ############################
def _myfloat(k):
    try:
        return float(k)
    except TypeError:
        return k

def _myint(k):
    try:
        return int(k)
    except TypeError:
        return k
#####################################################


class CMAP_Exception(Exception):
    pass


class Toppar(object):
    sections = ('bond', 'angle', 'dihedral', 'improper', 'cmap', 'nonbond',
                'nbfix', 'hbond')

    def __init__(self):
        # prm
        self.bond = None
        self.angle = None
        self.dihedral = None
        self.improper = None
        self.cmap = None
        self.nonbond = None
        self.nbfix = None
        self.hbond = None
        self.bond_opts = None
        self.angle_opts = None
        self.dihedral_opts = None
        self.improper_opts = None
        self.cmap_opts = None
        self.nonbond_opts = None
        self.nbfix_opts = None
        self.hbond_opts = None
        # rtf
        self.mass = None
        self.residue = None
        self.patch = None

    def merge(self, other):
        """Creates a new Toppar object, where rtf/prm objects from `self` are
        given priority over rtf/prm objects from `other`.
        """
        tmp_self = deepcopy(self)
        tmp_self.make_unique()

        tmp_other = deepcopy(other)
        tmp_other.make_unique()

        for section in self.sections:
            self_section = getattr(tmp_self, section)
            other_section = getattr(tmp_other, section)
            if self_section is None and other_section is None:
                setattr(tmp_self, section, None)
            elif self_section is None and other_section is not None:
                setattr(tmp_self, section, other_section)
            elif self_section is not None and other_section is None:
                pass
            else:
                if section == 'cmap':
                    if self_section != other_section:
                        raise CMAP_Exception("CMAPs are different between toppar")
                    setattr(tmp_self, section, getattr(tmp_self, section))
                else:
                    setattr(tmp_self, section, getattr(tmp_self, section) + getattr(tmp_other, section))
        return tmp_self

    def __add__(self, other):
        return self.merge(other)

    def make_unique(self):
        """Currently destroys original ordering."""
        for section in self.sections:
            if section == 'cmap':
                continue
            if getattr(self, section) is None:
                continue
            tmp = sorted(list(set(getattr(self, section))))

class PRM(ABCMeta):
    """Dummy metaclass, serves as a registrar and API definition."""

    # API Definition ################################################
    @abstractproperty
    def _sortkey(self):
        raise NotImplementedError

    @abstractmethod
    def _validate(self):
        raise NotImplementedError


class BasePRM(object):
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

class ImproperPRM(BasePRM):
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


class CmapPRM(BasePRM):
    __slots__ = ['text']

    def __init__(self, arg=None):
        self.text = arg

    @property
    def _sortkey(self):
        return ('CMAP_BLOCK',)

    def _validate(self):
        pass


class NonbondPRM(BasePRM):
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
    __slots__ = ['atom0', 'atom1', 'k', 'eq', 'ig14', 'k14', 'eq14']

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


class Residue(object):
    def __init__(self, name=None, charge=None, body=None):
        self.name = name
        self.charge = _myfloat(charge)
        self.body = body

    @property
    def _sortkey(self):
        return "a%r" % self.name

    # Magic Methods #################################################
    def __hash__(self):
        return hash(self._sortkey)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, self.name)

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

class Patch(Residue):
    @property
    def _sortkey(self):
        return "b%r" % self.name
