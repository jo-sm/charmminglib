"""
DOCME
"""
# fcp
# 10/26/2010


from IPython.external.pretty import pretty
from itertools import tee
from numpy import array, fromiter, float, dot, sin, cos
from numpy.linalg import eig, norm
from charmming.const.units import DEG2RAD
from charmming.tools import Property, expandPath, lowerKeys
from charmming.lib.metaatom import MetaAtom


class StructError(Exception):
    """
    Exception to raise when errors occur involving the BaseStruct class.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class BaseStruct(list):
    """
    This is the base class for all classes that are containers for
    "atom-like" objects.

    The current implementation uses `list` as a base class, this is not
    finalized, other candidates to derive from include:
        `np.array`
        `OrderedSet`
        `OrderedDict`
        `blist`

    Private Attributes
        `_autoFix`
    Properties
        `addr`              STUB
        `code`
        `com`
        `mass`
        `name`
    Public Methods
        `del_atoms`
        `rotate`            TODO
        `translate`         TODO
        `write`
    """
    def __init__(self, iterable=None, **kwargs):
        """
        kwargs:
            `autofix`       [True,False]
            `code`          string          # pdbcode
            `name`          string
        """
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
    def addr():
        doc =\
        """
        The `addr` property provides a human readable unique string
        representation for each `BaseStruct` instance.

        STUB
        """
        def fget(self):
            raise NotImplementedError
        return locals()

    @Property
    def code():
        doc =\
        """
        A string representing the Protein Data Bank code from which
        this `BaseStruct` object was derived.
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
        The center of mass in Angstrom.
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
        The total mass, in AMU.
        """
        def fget(self):
            return sum( ( atom.mass for atom in self ) )
        return locals()

    @Property
    def name():
        doc =\
        """
        DOCME
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
        self.translate(-1*self.com)

    def del_atoms(self, iterable):
        """
        Deletes, in place, one or more `Atoms` as specified by `iterable`.
        The `Atoms` themselves should go in the iterable which indicates
        the atoms to delete.

        >>> BaseStruct.del_atoms(BaseStruct[3:6])
        """
        for key in iterable:
            if not isinstance(key, MetaAtom):
                raise StructError('del_atoms: only objects from \
                                `MetaAtom` class may be deleted')
        del_atoms = set(iterable)
        del_indicies = [ i for i, atom in enumerate(self) if atom in del_atoms ]
        for index in reversed(del_indicies):
            self.pop(index)

    def find(self, addr=None, **kwargs):
        """
        Returns a BaseStruct object containing all "atom-like" objects which
        match the specified search criteria.

        Searches may be performed using an `addr` value as a string, or using
        the following kwargs:
            `chainid`
            `segtype`
            `resid`
            `atomnum`
            `atomtype`

        >>> self.find('a.pro.3')

        >>> self.find(chainid='a', segtype='pro', resid=3)
        """
        if addr is not None:
            tmp = addr.split('.')
            chainid = segtype = resid = atomnum = atomtype = None
            try:
                chainid = tmp[0]
                segtype = tmp[1]
                resid = tmp[2]
                atomnum = tmp[3]
                atomtype = tmp[4]
            except IndexError:
                pass
        elif kwargs:
            chainid = kwargs.get('chainid', None)
            segtype = kwargs.get('segtype', None)
            resid = kwargs.get('resid', None)
            atomnum = kwargs.get('atomnum', None)
            atomtype = kwargs.get('atomtype', None)
        else:
            return BaseStruct([], autofix=False)
        #
        try:
            resid = int(resid)
        except ValueError:
            pass
        except TypeError:
            pass
        try:
            atomnum = int(atomnum)
        except ValueError:
            pass
        except TypeError:
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
            iterator = ( atom for atom in iterator if atom.atomNum == resid )
        if atomtype:
            iterator = ( atom for atom in iterator if atom.atomType == atomtype )
        return BaseStruct(iterator, autofix=False)

    def find_byDistance(self, selection, distance):
        """
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
            [ xx, -yz, -xz],
            [-yz,  yy, -yz],
            [-xz, -yz,  zz]
            ])
        if eigen:
            return eig(I)
        else:
            return I

    def orient(self):
        """
        """
        self.center()
        self.rotateByMatrix(self.get_inertiaTensor(eigen=True)[1].transpose())

    def rotate(self, rotVector, angle, units='deg'):
        """
        TODO
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

    def rotateByMatrix(self, rotMatrix):
        """
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
        TODO
        """
        transVector = array(transVector)
        for atom in self:
            atom.cart += transVector

    def write(self, filename, **kwargs):
        """
        Writes a file containing the molecular information contained in the
        `BaseStruct` object.

        kwargs:
            `outformat`     ["charmm","pdborg","debug","xdebug","crd","xcrd"]
            `old_chainid`   [False,True]
            `old_segType`   [False,True]
            `old_resid`     [False,True]
            `old_atomNum`   [False,True]
            `ter`           [False,True]
            `end`           True if `outformat` in ["pdborg","charmm"]

        >>> thisSeg.write('~/1yjp.pdb',outformat='charmm',old_resid=True)
        """
        # kwargs
        kwargs = lowerKeys(kwargs)
        outFormat = kwargs.get('outformat', 'charmm')
        end = kwargs.get('end', None)
        ter = kwargs.get('ter', None)
        #
        filename = expandPath(filename)
        writeMe = []
        if outFormat in ['pdborg', 'charmm']:
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
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
            writeMe.append('END')
        # Write file
        writeMe = '\n'.join(writeMe)
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
