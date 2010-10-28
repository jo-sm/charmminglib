"""
DOCME
"""
# fcp
# 10/22/2010


from numpy import subtract, dot, cross, arccos, array
from numpy.linalg import norm
from charmming.const.units import RAD2DEG
from charmming.tools import Property, lowerKeys


class AtomError(Exception):
    """
    The exception to raise when errors occur involving BaseAtom, or
    derived classes.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class MetaAtom(object):
    """
    This isn't actually a metaclass per se, instead it is a framework
    from which all "atom-like" classes should be derived.  Common
    elements are predefined, other elements are left stubbed out, to be
    defined by subsequent subclasses.

    Methods, attributes and properties which are present here, but are
    left undefined have a "STUB" listing next to them.

    Class Attributes
        `_autoInFormat`     STUB
    Private Attributes
        `_autoFix`
        `_hash`
        `_index`
        `_text`
    Parse Definition
        `parse`             STUB
    Properties
        `addr`              STUB
        `addr0`
        `mass`
        `cart`
    Public Methods
        `Print`             STUB
        `calc_length`
        `calc_angle`
        `calc_dihedral`
        `calc_signedDihedral`
    Private Methods
        `_sort`             STUB
        `_set_hash`
    """

    _autoIndex = 0

    _autoInFormat = None
    """
    A string that defines the default input formatting for all class
    instances.

    STUB
    """

    def __init__(self, text=None, **kwargs):
        """
        DOCME
        """
        super(MetaAtom, self).__init__()
        # kwargs
        kwargs = lowerKeys(kwargs)
        commentChar = kwargs.pop('commentchar', '#')
        inFormat = kwargs.pop('informat', self.__class__._autoInFormat)
        self._index = kwargs.pop('index', MetaAtom._autoIndex)
        self._autoFix = kwargs.pop('autofix', True)
        # Main
        if text is None:
            self.cart = (0., 0., 0.)
        else:
            self._text  = text.lower().split(commentChar)[0].strip()
            self.parse(inFormat)
            self._addr0 = self.addr
        # Finally
        self._set_hash()
        MetaAtom._autoIndex += 1

####################
# Parse Definition #
####################

    def parse(self, inFormat):
        """
        A method that parses the `_text` attribute into instance
        variables.

        STUB
        """
        raise NotImplementedError

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        The `addr` property provides a human readable unique string
        representation for each `MetaAtom` instance.

        STUB
        """
        def fget(self):
            raise NotImplementedError
        return locals()

    @Property
    def addr0():
        doc =\
        """
        The `addr0` property gives the value of `addr` at instantization.

        This is a handy attribute to have for debugging, as `BaseAtom`
        instances are mutable, and thus their `addr` property can change
        over time, however `addr0` should not.
        """
        def fget(self):
            try:
                return self._addr0
            except AttributeError:
                return 'Index: %d' % self._index
        return locals()

    @Property
    def cart():
        doc =\
        """
        The `cart` property is a `numpy.array` which contains the
        instances' cartesian data.

        This entry is given in units of Angstrom, and it should be
        within the range (-10000,10000) due to limitations of CHARMM.
        """
        def fget(self):
            return self._cart
        def fset(self, value):
            value = map(float, value)
            for i,taco in enumerate(value):
                if taco < -10000.:
                    if self._autoFix:
                        value[i] = -10000.
                    else:
                        raise AtomError('cart %s: %9.2f is less than -10000.0' %
                                        (self.addr0, taco))
                elif taco > 10000.:
                    if self._autoFix:
                        value[i] = 10000.
                    else:
                        raise AtomError('cart %s: %9.2f is greater than 10000.0' %
                                        (self.addr0, taco))
            self._cart = array(value)
        return locals()

    @Property
    def mass():
        doc =\
        """
        The `mass` property is a float representing the mass in units
        of AMU.
        """
        def fget(self):
            return self._mass
        def fset(self, value):
            value = float(value)
            if value <= 0:
                if self._autoFix:
                    value = 0
                else:
                    raise AtomError('mass %s: %7.2f is invalid' %
                                    (self.addr0, value))
            self._mass = value
        return locals()

##################
# Public Methods #
##################

    def Print(self, **kwargs):
        """
        The method responsible for formatting the instance variables
        into strings of text to be written to external files.

        At a bare minimum, the kwarg `outformat` should be defined.

        STUB
        """
        raise NotImplementedError

    def calc_length(self, other):
        """
        Returns the cartesian distance between two `BaseAtom` objects.
        """
        # Localize data
        i = self.cart
        j = other.cart
        # Do work
        return norm(subtract(i, j))

    def calc_angle(self, other1, other2, units='deg'):
        """
        Returns the cartesian angle between three `BaseAtom` objects.

        By default the angle is output in radians, however with the
        `units="deg"` keyword, it will return degrees.
        """
        # Localize data
        i = self.cart
        j = other1.cart
        k = other2.cart
        # Do work
        Rji = subtract(j, i)
        Rji /= norm(Rji)
        Rjk = subtract(j, k)
        Rjk /= norm(Rjk)
        result = arccos(dot(Rji, Rjk))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

    def calc_dihedral(self, other1, other2, other3, units='deg'):
        """
        Returns the dihedral angle between four `BaseAtom` objects.

        By default the angle is output in radians, however with the
        `units="deg"` keyword, it will return degrees.
        """
        # Localize data
        i = self.cart
        j = other1.cart
        k = other2.cart
        l = other3.cart
        # Do work
        Rij = subtract(i, j)
        Rjk = subtract(j, k)
        Rkl = subtract(k, l)
        Nijk = cross(Rij, Rjk)
        Nijk /= norm(Nijk)
        Njkl = cross(Rjk, Rkl)
        Njkl /= norm(Njkl)
        result = arccos(dot(Nijk, Njkl))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

    def calc_signedDihedral(self,other1,other2,other3,units='deg'):
        """
        Returns the phase corrected dihedral angle between four
        `BaseAtom` objects.

        By default the angle is output in radians, however with the
        `units="deg"` keyword, it will return degrees.
        """
        # Localize data
        i = self.cart
        j = other1.cart
        k = other2.cart
        l = other3.cart
        # Do work
        Rij = subtract(i, j)
        Rjk = subtract(j, k)
        Rkl = subtract(k, l)
        Nijk = cross(Rij, Rjk)
        Nijk /= norm(Nijk)
        Njkl = cross(Rjk, Rkl)
        Njkl /= norm(Njkl)
        # Phase
        sign = -1 if dot(cross(Nijk, Njkl)) > 0 else 1
        result = sign * arccos(dot(Nijk, Njkl))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

###################
# Private Methods #
###################

    def _sort(self):
        """
        A scoring method that determines sorting order, for rich
        comparisons and containerClass.sort() methods.

        STUB
        """
        raise NotImplementedError

    def _set_hash(self):
        """
        Method for setting the instance's `_hash` value, which is in
        turn used with the built-in `hash` function.
        """
        self._hash = sum(map(hash,self.cart))

###################
# Special Methods #
###################

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,self.Print(outFormat='repr'))

    def __str__(self):
        return self.addr

    def __hash__(self):
        return self._hash

    def __eq__(self,other):
        return self._hash == other._hash and self.__class__.__name__ == \
                other.__class__.__name__

    def __ne__(self,other):
        return self._hash != other._hash

    def __lt__(self,other):
        return self._sort() < other._sort()

    def __le__(self,other):
        return self._sort() <= other._sort()

    def __gt__(self,other):
        return self._sort() > other._sort()

    def __ge__(self,other):
        return self._sort() >= other._sort()
