from numpy import subtract,dot,cross,arccos,array
from numpy.linalg import norm
from charmming.const.units import RAD2DEG
from charmming.tools import Property


class AtomError(Exception):
    """
    Exception to raise when errors occur involving the Atom class.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class BaseAtom(object):
    """
    In order for a class to inherit properly from this one, it must implement
    everything that is indicated as a stub.

    Private Attributes
        _text           string
        _index          [0,int]
        _autoFix        [True,False]
        _hash
    Parse Definition
        parse           STUB
    Properties
        addr            STUB
        addr0           STUB
        mass
        cart
    Public Methods
        Print           STUB
        calc_length
        calc_angle
        calc_dihedral
        calc_signedDihedral
    Private Methods
        _sort           STUB
    """
    def __init__(self,text=None,**kwargs):
        super(BaseAtom, self).__init__()
        # kwargs
        cc              = kwargs.pop('cc','#')
        self._index     = kwargs.pop('index',0)
        self._autoFix   = kwargs.pop('autoFix',True)
        if kwargs.keys():
            raise TypeError('Unprocessed kwargs(%r)' % kwargs.keys())
        if text is None:
            self.cart = (0.,0.,0.)
        else:
            self._text  = text.lower().split(cc)[0].strip()

####################
# Parse Definition #
####################

# STUB
    def parse(self):
        raise NotImplementedError

##############
# Properties #
##############

# STUB
    @Property
    def addr():
        doc = "The addr property."
        def fget(self):
            raise NotImplementedError
        return locals()

# STUB
    @Property
    def addr0():
        doc = "The addr0 property."
        def fget(self):
            return self._index
        return locals()

    @Property
    def cart():
        doc = "The cart property, define it as a 3-tuple."
        def fget(self):
            return self._cart
        def fset(self, value):
            value = map(float,value)
            for i,taco in enumerate(value):
                if taco < -10000.:
                    if self._autoFix:
                        value[i] = -10000.
                    else:
                        raise AtomError('cart %s: %9.2f is less than -10000.0' % (self.addr0,taco))
                elif taco > 10000.:
                    if self._autoFix:
                        value[i] = 10000.
                    else:
                        raise AtomError('cart %s: %9.2f is greater than 10000.0' % (self.addr0,taco))
            self._cart  = array(value)
        return locals()

    @Property
    def mass():
        doc = "The mass property."
        def fget(self):
            return self._mass
        def fset(self, value):
            value = float(value)
            if value < 0.:
                if self._autoFix:
                    value = 'badMass'
                else:
                    raise AtomError('mass %s: %7.2f is less than 0.0' % (self.addr0,value))
            elif value > 1000.:
                if self._autoFix:
                    value = 'badMass'
                else:
                    raise AtomError('mass %s: %7.2f is greater than 1000.0' % (self.addr0,value))
            self._mass = value
        return locals()

##################
# Public Methods #
##################

# STUB
    def Print(self):
       raise NotImplementedError

    def calc_length(self,other):
        """
        Returns the cartesian distance between two BaseAtom objects.
        """
        # Localize data
        i = self.cart
        j = other.cart
        # Do work
        return nor(sub(i,j))

    def calc_angle(self,other1,other2,units='deg'):
        """
        Returns the bond calc_angle between three BaseAtom objects.  You can optionally
        specify the units as the 3rd argument with 'rad', default is 'deg'.
        """
        # Localize data
        i = self.cart
        j = other1.cart
        k = other2.cart
        # Do work
        Rji = subtract(j,i)
        Rji /= norm(Rji)
        Rjk = subtract(j,k)
        Rjk /= norm(Rjk)
        result = arccos(dot(Rji,Rjk))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

    def calc_dihedral(self,other1,other2,other3,units='deg'):
        """
        Returns the torsion calc_angle between four BaseAtom objects.  You can optionally
        specify the units as the 4th argument with 'rad', default is 'deg'.
        """
        # Localize data
        i = self.cart
        j = other1.cart
        k = other2.cart
        l = other3.cart
        # Do work
        Rij = subtract(i,j)
        Rjk = subtract(j,k)
        Rkl = subtract(k,l)
        Nijk = cross(Rij,Rjk)
        Nijk /= norm(Nijk)
        Njkl = cross(Rjk,Rkl)
        Njkl /= norm(Njkl)
        result = arccos(dot(Nijk,Njkl))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

    def calc_signedDihedral(self,other1,other2,other3,units='deg'):
        """
        Returns the signed torsion calc_angle between four BaseAtom objects.  You can optionally
        specify the units as the 4th argument with 'rad', default is 'deg'.
        """
        # Localize data
        i = self.cart
        j = other1.cart
        k = other2.cart
        l = other3.cart
        # Do work
        Rij = subtract(i,j)
        Rjk = subtract(j,k)
        Rkl = subtract(k,l)
        Nijk = cross(Rij,Rjk)
        Nijk /= norm(Nijk)
        Njkl = cross(Rjk,Rkl)
        Njkl /= norm(Njkl)
        result = ( -1 if dot(cross(Nijk,Njkl)) > 0 else 1 ) * arccos(dot(Nijk,Njkl))
        # Units
        if units == 'deg':
            return result * RAD2DEG
        else:
            return result

###################
# Private Methods #
###################

# STUB
    def _sort(self):
        raise NotImplementedError

    def _set_hash(self):
        """
        This should be called at the very end of derived classes __init__.
        """
        self._hash = sum(map(hash,self.cart))

###################
# Special Methods #
###################

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,self.Print(format='repr'))

    def __str__(self):
        return self.addr

    def __hash__(self):
        return self._hash

    def __eq__(self,other):
        return self._hash == other._hash and self.__class__.__name__ == other.__class__.__name__

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
