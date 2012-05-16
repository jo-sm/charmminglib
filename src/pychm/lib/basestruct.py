"""
:Author: fcp
:Date: 10/26/2010
"""


from copy import deepcopy
from itertools import tee
from numpy import array, fromiter, float, dot, sin, cos, ones
from numpy.linalg import eig, norm
from pychm.const.units import DEG2RAD
from pychm.tools import Property, expandPath, lowerKeys
from pychm.lib.metaatom import MetaAtom


class StructError(Exception):
    """
    Exception to raise when errors occur involving the :class:`BaseStruct` class.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class BaseStruct(list):
    """
    This is the base class for all classes that are containers for
    "atom-like" objects.

    The current implementation uses :class:`list` as a base class, this is not
    finalized, other candidates to derive from include:

        | :class:`numpy.array`
        | :class:`OrderedSet`
        | :class:`OrderedDict`
        | :class:`blist`

    The tradeoffs between the various classes are many and varied, but they
    basically boil down to performance versus flexibility.  A list based
    solution gives quite a bit of flexibility, the one potential drawbackc
    is accessing data members without using an *index*.  The :meth:`find`
    method attempts to address this short-coming.

    The most fundamental purpose of all :class:`BaseStruct`-like objects
    is to provide a clean interface for creating iterators over arbitrary
    :class:`Atom`-like objects.  The :class:`BaseStruct` class itself provides
    methods for atom selection :meth:`find`, :meth:`findByDistance`,
    :meth:`__add__` and :meth:`__sub__`.

    :class:`BaseStruct` objects may be manipulated using most of the
    standard :class:`list` based methods.  However note that the
    :meth:`__add__` method has been changed, and the :meth:`__mul__` and
    :meth:`__setslice__` has been explicitly removed.  Other methods such
    as :meth:`list.sort` may behave differently than expected, please see
    :mod:`pychm.lib.baseatom` for relevant documentation.  To assist
    in atom deletion, the cumberson :meth:`BaseStruct.del_atoms` has been
    added.

    Coordinate manipulations are also possible using :class:`BaseStruct`
    methods.  One can perform arbitrary translations and rotations on a
    given atom selection, or simply choose to use canonical molecule
    fixed axies using the :meth:`orient` method.

    Finally, molecular data can be written using the aptly named :meth:`write`
    method.  Any output format which is supplied to the :meth:`Print` method
    for the container's :class:`Atom`-like objects is valid.

    Default values for *kwargs* are listed first.

    **kwargs:**
        | ``autofix``       [True,False]
        | ``code`` :: pdbcode
        | ``name``

    **Special Methods:**

        | *Removed:*

            | ``__imul__``
            | ``__mul__``
            | ``__setslice__``
            | ``__isub__``
            | ``__iadd__``

        | *Changed:*

            | ``__sub__``
            | ``__add__``

    :TODO:
        | ``rotateByEuler``
    """
    def __init__(self, iterable=None, **kwargs):
        # kwargs
        kwargs = lowerKeys(kwargs)
        self._autoFix = kwargs.get('autofix', True)
        self._code = kwargs.get('code', 'None')
        self._name = kwargs.get('name', 'None')
        #
        self._dict = {}
        # Gatekeeper
        if iterable is None:
            super(BaseStruct, self).__init__()
        else:
            if self._autoFix:
                iterable = ( key for key in iterable if isinstance(key, MetaAtom) )
            else:
                iterable,iterchk = tee(iterable, 2)
                for key in iterchk:
                    if not isinstance(key, MetaAtom):
                        raise StructError('Only objects derived from `MetaAtom`\
                                        class may be added to a BaseStruct object')
            super(BaseStruct, self).__init__()
            self.extend(iterable)

##############
# Properties #
##############

    @Property
    def code():
        doc =\
        """
        A property for the string representing the Protein Data Bank
        code from which this ``BaseStruct`` object was derived.
        """
        def fget(self):
            return self._code
        def fset(self, value):
            self._code = value
        return locals()

    @Property
    def com():
        doc =\
        """
        The property for the center of mass of the collection of
        atoms, in Angstroms.  Read only.
        """
        def fget(self):
            result = array([ atom.mass * atom.cart for atom in self ])
            result = result.sum(axis=0)
            return result / self.mass
        return locals()

    @Property
    def mass():
        doc =\
        """
        The property for the total mass of the collection of atoms,
        in AMU.  Read only.
        """
        def fget(self):
            return sum( ( atom.mass for atom in self ) )
        return locals()

    @Property
    def name():
        doc =\
        """
        The property for an arbitrary string label for the
        :class:`BaseStruct`.
        """
        def fget(self):
            return self._name
        def fset(self, value):
            self._name = value
        return locals()

##################
# Public Methods #
##################

    def center(self):
        """
        A method to translate the collection of atoms, such that their
        new collective center of mass is at the origin.
        """
        self.translate(-1*self.com)

    def del_atoms(self, iterable):
        """
        Deletes, in place, one or more :class:`Atom`-like objects as
        specified by `iterable`.  `Iterable` should iterate over the
        atom-like objects to be deleted.

        For example...

        >>> taco = BaseStruct(some_iterable_of_atoms)
        >>> taco.del_atoms(taco[3:6])

        will remove the 4th, 5th and 6th atoms from `taco`.
        """
        for key in iterable:
            if not isinstance(key, MetaAtom):
                raise StructError('del_atoms: only objects from \
                                `MetaAtom` class may be deleted')
        del_atoms = set(iterable)
        del_indicies = [ i for i, atom in enumerate(self) if atom in del_atoms ]
        for index in reversed(del_indicies):
            self.pop(index)

    def find(self, **kwargs):
        """
        Returns a :class:`BaseStruct` object containing all "atom-like"
        objects which match the specified search criteria.

        **Note:**
            This method will always return a :class:`BaseStruct` object
            even if it is called by a child class.  If you want to then use
            methods from the child class on this selection, you need to
            then reinstantize that child class using this result as an
            iterator.  Sub-optimal for sure.

        **kwargs:**
            | ``chainid``
            | ``segtype``
            | ``resid``
            | ``atomnum``
            | ``atomtype``
            | ``resName``

        >>> self.find(chainid='a', segtype='pro', resid=3)

        """
        if kwargs:
            kwargs = lowerKeys(kwargs)
            chainid = kwargs.get('chainid', None)
            segtype = kwargs.get('segtype', None)
            resid = kwargs.get('resid', None)
            atomnum = kwargs.get('atomnum', None)
            atomtype = kwargs.get('atomtype', None)
            resname = kwargs.get('resname', None)
            chainid0 = kwargs.get('chainid0', None)
            segtype0 = kwargs.get('segtype0', None)
            resid0 = kwargs.get('resid0', None)
            atomnum0 = kwargs.get('atomnum0', None)
            atomtype0 = kwargs.get('atomtype0', None)
            resname0 = kwargs.get('resname0', None)
        else:
            return BaseStruct([], autofix=False)
        #
        try:
            resid = int(resid)
        except (ValueError, TypeError):
            pass
        try:
            atomnum = int(atomnum)
        except (ValueError, TypeError):
            pass
        try:
            resid0 = int(resid0)
        except (ValueError, TypeError):
            pass
        try:
            atomnum0 = int(atomnum0)
        except (ValueError, TypeError):
            pass
        #
        iterator = ( atom for atom in self)
        if chainid:
            iterator = ( atom for atom in iterator if atom.chainid == chainid )
        if segtype:
            iterator = ( atom for atom in iterator if atom.segType == segtype )
        if resid:
            iterator = ( atom for atom in iterator if atom.resid == resid )
        if atomnum:
            iterator = ( atom for atom in iterator if atom.atomNum == atomnum )
        if atomtype:
            iterator = ( atom for atom in iterator if atom.atomType == atomtype )
        if resname:
            iterator = ( atom for atom in iterator if atom.resName == resname )
        if chainid0:
            iterator = ( atom for atom in iterator if atom.chainid0 == chainid0 )
        if segtype0:
            iterator = ( atom for atom in iterator if atom.segType0 == segtype0 )
        if resid0:
            iterator = ( atom for atom in iterator if atom.resid0 == resid0 )
        if atomnum0:
            iterator = ( atom for atom in iterator if atom.atomNum0 == atomnum0 )
        if atomtype0:
            iterator = ( atom for atom in iterator if atom.atomType0 == atomtype0 )
        if resname0:
            iterator = ( atom for atom in iterator if atom.resName0 == resname0 )
        return BaseStruct(iterator, autofix=False)

    def find_byDistance(self, selection, distance):
        """
        Returns a :class:`BaseStruct` which is a subset of the
        :class:`BaseStruct` instance upon which this method is called.

        By specifying a `selection` (a :class:`BaseStruct`) and a
        `distance` (a :class:`float`), any atoms in `self` within
        `distance` of `selection` are returned in a new :class:`BaseStruct`.

        >>> taco.find_byDistance(taco[0:1],2)
        [Atom(      A.PRO.1.1 ::  N    LEU      10.451   -1.168   -1.716),
         Atom(      A.PRO.1.2 ::  CA   LEU       9.226   -0.710   -1.004),
         Atom(      A.PRO.1.9 ::  H    LEU      11.269   -1.350   -1.211)]

        This returns all atoms in `taco` that are within 2 Angstroms of
        the first two atoms in `taco`.
        """
        def proximal(atom, selection, distance):
            for atom_i in selection:
                if atom.calc_length(atom_i) <= distance:
                    return True
            return False
        iterator = ( atom for atom in self if proximal(atom, selection, distance) )
        return BaseStruct(iterator, autofix=False)

    def get_inertiaTensor(self, eigen=False):
        """
        Returns a 3 by 3 :class:`numpy.array` corresponding to the
        inertia tensor of the atom selection.  Alternatively, if
        `eigen` is set to ``True``, :meth:`numpy.linalg.eig` is called
        upon the inertia tensor, and a tuple containing the eigen values
        and eigen vectors are returned directly.

        >>> taco.get_inertiaTensor(eigen=True)
        (array([ 250566.8631789 ,  429622.74865092,  377518.3392818 ]),
         array([[-0.99381953,  0.08612146, -0.07004171],
                [ 0.0908166 ,  0.26793577, -0.9591469 ],
                [ 0.06383645,  0.95957987,  0.27410106]]))

        This returns the eigen values are the first element in the tuple,
        and the eigen vectors, as column vectors are the second element
        in the tuple.
        """
        xx, yy, zz, xy, xz, yz = (0., 0., 0., 0., 0., 0.)
        for atom in self:
            x, y, z = atom.cart
            m = atom.mass
            #
            xx += m*(y*y+z*z)
            yy += m*(x*x+z*z)
            zz += m*(x*x+y*y)
            xy += m*x*y
            xz += m*x*z
            yz += m*y*z
        #
        I = array([
            [ xx, -xy, -xz],
            [-xy,  yy, -yz],
            [-xz, -yz,  zz]
            ])
        if eigen:
            return eig(I)
        else:
            return I

    def get_rmsd(self, other, orient=False, mass=False):
        """
        Get the root mean squared deviation between two struct objects.
        Flags:
            'orient'    -- Defaults to False, if set orients the two structures
                        to each other before calculating RMSD.
            'mass'      -- Defaults to False, if set calculates a mass weighted
                        RMSD instead of a non-weighted RMSD.
        """
        # validate data
        assert len(self) == len(other)
        if mass:
            assert abs(self.mass - other.mass) < 0.001
        # make copies so we dont change original objects
        tmp_self = deepcopy(self)
        tmp_other = deepcopy(other)
        if orient:
            tmp_self.orient()
            tmp_other.orient()
        #
        iterator = ( crd for atom in tmp_self for crd in atom.cart )
        self_crd = fromiter(iterator, float)
        self_crd.resize((len(self), 3))
        iterator = ( crd for atom in tmp_other for crd in atom.cart )
        other_crd = fromiter(iterator, float)
        other_crd.resize((len(other), 3))
        # weighting
        if mass:
            iterator = ( atom.mass for atom in tmp_self )
            weight = fromiter(iterator, float)
            weight = array((weight, weight, weight)).T
        else:
            weight = ones(self_crd.shape)
        # calc rms
        diff_crd = (self_crd - other_crd) * weight
        return ((diff_crd**2).mean())**0.5

    def get_span(self):
        """
        Returns a 3-tuple which represents the span (max - min)for
        the x, y and z coordinates.
        """
        iterator = ( crd for atom in self for crd in atom.cart )
        tmp = fromiter(iterator, float)
        tmp.resize((len(self), 3))
        x = tmp[:,0]
        y = tmp[:,1]
        z = tmp[:,2]
        return (x.max() - x.min(), y.max() - y.min(), z.max() - z.min())

    def orient(self):
        """
        This method centers the atom selection about the origin, and
        rotates it such that the priciples axies of rotation (A, B, C)
        are now coincident with (z, y, x).
        """
        self.center()
        self.rotateByMatrix(self.get_inertiaTensor(eigen=True)[1].transpose())

    def rotate(self, rotVector, angle, units='deg'):
        """
        Rotate an atom selection about an axis defined by a cartesian
        vector `rotVector` and about an `angle`.  By default, angle
        is specified in degrees, you may also use radians.
        """
        # axis
        assert len(rotVector) == 3
        rotVector = array(rotVector)
        rotVector /= norm(rotVector)
        x, y, z = rotVector
        # angle
        if units == 'deg':
            t = angle * DEG2RAD
        else:
            t = angle
        # rotation matrix
        ct = cos(t)
        ct1 = 1 - cos(t)
        st = sin(t)
        #
        R = array([
            [ct+x*x*ct1,   x*y*ct1-z*st, x*z*ct1+y*st],
            [y*x*ct1+z*st, ct+y*y*ct1,   y*z*ct1-x*st],
            [z*x*ct1-y*st, z*y*ct1+x*st, ct+z*z*ct1  ]
            ])
        #
        self.rotateByMatrix(R)

    def rotateByEuler(self, alpha, beta, gamma):
        """
        Rotate an atom selection using Euler angles `alpha`, `beta` and
        `gamma`.

        **TODO**
        """


    def rotateByMatrix(self, rotMatrix):
        """
        Rotate an atom selection by an arbitrary rotation matrix.
        """
        # create coordinate matrix
        iterator = ( crd for atom in self for crd in atom.cart )
        tmp = fromiter(iterator, float)
        tmp.resize((len(self), 3))
        # rotate
        tmp = dot(tmp, rotMatrix.transpose())
        # unpack results
        for i, atom in enumerate(self):
            atom.cart = tmp[i]

    def translate(self, transVector):
        """
        Translate an atom selection by an arbitrary translation vector.
        """
        assert len(transVector) == 3
        transVector = array(transVector)
        for atom in self:
            atom.cart += transVector

    def write(self, filename, **kwargs):
        """
        Writes a file containing the molecular information contained in
        :class:`BaseStruct` object.

        **kwargs:**
            | ``outformat``     ["charmm","pdborg","debug","xdebug","crd","xcrd"]
            | ``old_chainid``   [False,True]
            | ``old_segType``   [False,True]
            | ``old_resid``     [False,True]
            | ``old_atomNum``   [False,True]
            | ``ter``           [False,True]
            | ``end``           True if `outformat` in ["pdborg","charmm"]
            | ``append``        [False,True]

        >>> thisSeg.write('~/1yjp.pdb',outformat='charmm',old_resid=True)
        """
        # kwargs
        kwargs = lowerKeys(kwargs)
        outFormat = kwargs.get('outformat', 'charmm')
        end = kwargs.get('end', None)
        ter = kwargs.get('ter', None)
        append = kwargs.get('append', False)
       
        #
        filename = expandPath(filename)
        writeMe = []
        if outFormat in ['pdborg', 'charmm']:
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
                if ter is None:
                    ter = True
                if end is None:
                    end = True
        elif outFormat in ['debug', 'xdebug']:
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
        elif outFormat in ['crd', 'cor', 'card', 'short', 'shortcard',
                        'xcrd', 'xcor', 'xcard', 'long', 'longcard']:
            writeMe.append('*')
            writeMe.append('   %d' % len(self))
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
        # TER/END
        if ter:
            writeMe.append('TER')
        if end:
            writeMe.append('END\n')
        # Write file
        writeMe = '\n'.join(writeMe)
        if append:
            writeTo = open(filename,'a')
        else:
            writeTo = open(filename,'w')
        writeTo.write(writeMe)
        writeTo.close()

###################
# Special Methods #
###################

    def __imul__(self):
        raise NotImplementedError

    def __mul__(self):
        raise NotImplementedError

    def __setslice__(self):
        raise NotImplementedError

    def __isub__(self):
        raise NotImplementedError

    def __sub__(self, other, **kwargs):
        if not isinstance(other, BaseStruct):
            raise TypeError("unsupported operand type for '%s' and '%s'" %
                            (type(self), type(other)))
        iterator = ( atom for atom in self if atom not in other )
        return BaseStruct(iterator, **kwargs)

    def __iadd__(self):
        raise NotImplementedError

    def __add__(self, other):
        if not isinstance(other, BaseStruct):
            raise TypeError("unsupported operand type for '%s' and '%s'" %
                            (type(self), type(other)))
        return super(BaseStruct, self).__add__(other.__sub__(self))

    def __getitem__(self, key):
        try:
            return super(BaseStruct, self).__getitem__(key)
        except TypeError:
            for atom in self:
                if atom.addr == key:
                    return atom
            raise KeyError('No atom with the specified addr value: "%s"' % key)
