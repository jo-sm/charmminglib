"""
:Author: fcp
:Date: 10/22/2010
"""


from numpy import subtract, dot, cross, arccos, array, cos, sin
from numpy.linalg import norm
from charmming.const.units import RAD2DEG, DEG2RAD
from charmming.tools import Property, lowerKeys


class AtomError(Exception):
    """
    The exception to raise when errors occur involving ``MetaAtom``, or
    derived classes.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class MetaAtom(object):
    """
    A class which, by itself, doesn't actually do anything!

    This isn't actually a metaclass, instead it is a general framework
    for building ``atom``-like classes.  Accordingly, it should probably
    renamed `PseudoClass`, or something fancy like that.  Common
    elements are defined herein.  Others, which must be defined, yet
    which are disimilar, are left stubbed out, to be defined by subsequent
    subclassing.  Appropriately, this is denoted with a **STUB**.

    Because this class is explicitly designed for hacking and tweaking,
    and because the abstraction is a bit egregious, the documentation will
    attempt be more thorough than other sections, and will discuss features
    that might not be of interest to many developers.


    **Class Attributes:**
        | ``_autoIndex``
        | ``_autoInFormat``     **STUB**
        | ``_properties``

    **Private Attributes:**
        | ``_autoFix``
        | ``_hash``
        | ``_index``
        | ``_text``

    **Parse Definition:**
        | ``parse``             **STUB**

    **Properties:**
        | ``addr``              **STUB**
        | ``addr0``
        | ``mass``
        | ``cart``

    **Public Methods:**
        | ``Print``             **STUB**

    **Private Methods:**
        | ``_init_null``
        | ``_sort``             **STUB**
        | ``_set_hash``

    **kwargs:**
        | ``commentchar`` :: Defaults to '#'
        | ``informat`` :: Defaults to ``self.__class__._autoInFormat``
            which may set as a ``MetaAtom`` subclass class variable.
        | ``index`` :: A potentially useful way of differentiating
            ``Atom``-like objects during initialization.  For example,
            when generating error and debugging messages of *atoms*
            that fail to initialize properly.  This value defaults to
            ``MetaAtom._autoIndex``, which itself is just a counter of
            ``MetaAtom`` objects.
        | ``atomfix`` :: Defaults to ``True``
    """

    _autoIndex = 0
    """
    Keeps track of how many class instances exist, and serves as the
    default indexing value.
    """

    _autoInFormat = None
    """
    A string that defines the default input formatting for all class
    instances.

    **STUB**
    """

    _properties = {
        'cart': array((0., 0., 0.)),
        'mass': 0
    }
    """
    A dictionary which tells the class which properties (keys) it is
    associated with.  The corresponding values serve as default values
    for said properties.
    """

    def __init__(self, text=None, **kwargs):
        super(MetaAtom, self).__init__()
        # kwargs
        kwargs = lowerKeys(kwargs)
        commentChar = kwargs.get('commentchar', '#')
        inFormat = kwargs.get('informat', self.__class__._autoInFormat)
        self._index = kwargs.get('index', MetaAtom._autoIndex)
        self._autoFix = kwargs.get('autofix', True)
        # Main
        if text is None:
            self._init_null()
        elif isinstance(text, str):
            # card format files are whitespace sensitive
            if inFormat == 'shortcard' or inFormat == 'longcard':
                self._text  = text.lower()
            else:
                self._text  = text.lower().split(commentChar)[0].strip()

            self.parse(inFormat)
            self._addr0 = self.addr
        elif isinstance(text, MetaAtom):
            selfProp = set(self.__class__._properties.iterkeys())
            otherProp = set(text._properties.iterkeys())
            for key in selfProp.intersection(otherProp):
                setattr(self, key, getattr(text, key))
            for key in selfProp - otherProp:
                setattr(self, key, self.__class__._properties[key])
        else:
            raise TypeError("""Invalid input, initialization requires `None`,
                            `str` or `MetaAtom` type.""")
        # Finally
        self._set_hash()
        MetaAtom._autoIndex += 1

####################
# Parse Definition #
####################

    def parse(self, inFormat):
        """
        A method that parses the ``_text`` attribute into instance
        variables.

        **STUB**
        """
        raise NotImplementedError

##############
# Properties #
##############

    @Property
    def addr():
        doc =\
        """
        A ``property`` that provides a human readable string
        representation unique to each ``MetaAtom`` instance.

        **STUB**
        """
        def fget(self):
            raise NotImplementedError
        return locals()

    @Property
    def addr0():
        doc =\
        """
        A ``property`` that gives the value of ``addr`` at
        instantization.

        This is a handy attribute to have for debugging, as ``MetaAtom``
        instances are mutable, so ``addr`` can change over time,
        however ``addr0`` should not.
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
        A ``property`` that returns a :class:`numpy.array` which contains the
        instances' cartesian data in ``array([x, y, z])`` format.

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
        A ``property`` that returns a float representing the atomic
        mass in units of AMU.
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
        into strings of text to be written to external files or read
        in the interactive terminal.

        At a bare minimum, the kwarg ``outformat`` should be defined.

        **STUB**
        """
        raise NotImplementedError

    def calc_length(self, other):
        """
        Returns the cartesian distance between two ``MetaAtom`` objects.
        """
        # Localize data
        i = self.cart
        j = other.cart
        # Do work
        return norm(subtract(i, j))

    def calc_angle(self, other1, other2, units='rad'):
        """
        Returns the valence angle between three ``MetaAtom`` objects.

        By default the angle is output in radians, however with the
        ``units="deg"`` keyword, it will return degrees.
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

    def calc_dihedral(self, other1, other2, other3, units='rad'):
        """
        Returns the dihedral angle between four ``MetaAtom`` objects.

        By default the angle is output in radians, however with the
        ``units="deg"`` keyword, it will return degrees.
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

    def calc_signedDihedral(self, other1, other2, other3, units='rad'):
        """
        Returns the phase corrected dihedral angle between four
        ``MetaAtom`` objects.

        By default the angle is output in radians, however with the
        ``units="deg"`` keyword, it will return degrees.
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
        sign = -1 if dot(cross(Nijk, Njkl), Rjk) > 0 else 1
        result = sign * arccos(dot(Nijk, Njkl))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

    def rotate(self, axis, angle, units='rad'):
        """
        Rotate a ``MetaAtom`` through an `axis` specified by a
        cartesian vector and an `angle`.

        By default the angle is intput in radians, however with the
        ``units="deg"`` keyword, it will accept degrees.
        """
        # axis
        assert len(axis) == 3
        axis = array(axis)
        axis /= norm(axis)
        x, y, z = axis
        # angle
        if units == 'deg':
            t = angle * DEG2RAD
        else:
            t = angle
        #
        ct = cos(t)
        ct1 = 1 - cos(t)
        st = sin(t)
        #
        R = array([
            [ct+x*x*ct1,   x*y*ct1+z*st, x*z*ct1-y*st],
            [y*x*ct1-z*st, ct+y*y*ct1,   y*z*ct1+x*st],
            [z*x*ct1+y*st, z*y*ct1-x*st, ct+z*z*ct1  ]
            ])
        self.cart = dot(self.cart, R)

    def translate(self, i, j, k):
        pass

###################
# Private Methods #
###################

    def _init_null(self):
        """
        If the constructor is called without a `text` argument, this
        function defines the initial variable values of the Atom like
        object.
        """
        for key, value in self.__class__._properties.iteritems():
            setattr(self, key, value)

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
        for key, value in self.__class__._properties.iteritems():
            try:
                self._hash += hash(value)
            except TypeError:
                pass

###################
# Special Methods #
###################

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,self.Print(outFormat='repr'))

    def __str__(self):
        return self.addr

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        for i in xrange(3):
            if self._cart[i] != other._cart[i]:
                return False
        for key in self._properties.iterkeys():
            try:
                if getattr(self,key) != getattr(other,key):
                    return False
            except ValueError:
                pass
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self._sort() < other._sort()

    def __le__(self, other):
        return self._sort() <= other._sort()

    def __gt__(self, other):
        return self._sort() > other._sort()

    def __ge__(self, other):
        return self._sort() >= other._sort()
