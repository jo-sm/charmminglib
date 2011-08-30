"""
"""


import numpy
import scipy.ndimage.filters as filters
from copy import deepcopy
from pychm.tools import Property
import pychm.emap.ext
import itertools



def pairwise_add(iterable):
    def addTwo(iterable):
        def key(x, s=2, a=[-1]):
            r = a[0] = a[0] + 1
            return r // s
        for k, g in itertools.groupby(iterable, key):
            yield sum(g)
    iterable = list(iterable)
    while len(iterable) > 1:
        iterable = list(addTwo(iterable))
    return iterable[0]


class EMap(object):
    """
    """

    dx = 3.
    dy = 3.
    dz = 3.
    diff = numpy.array([dx, dy, dz])

    def __init__(self, mol=None):
        super(EMap, self).__init__()
        #
        self.basis = numpy.identity(3)
        self._cartArray = None
        self.chargeArray = None
        self.coreArray = None
        self.laplacianArray = None
        # Operations
        self.translations = []
        self.rotations = []

    @Property
    def cartArray():
        doc =\
        """
        DOCME
        """
        def fget(self):
            trans, rot = self.state
            tmp = numpy.dot(self._cartArray, rot.transpose())
            return tmp + trans
        return locals()

    @Property
    def state():
        doc =\
        """
        DOCME
        """
        def fget(self):
            trans0 = numpy.zeros(3)
            for trans in self.translations:
                trans0 += trans
            rot0 = numpy.identity(3)
            for rot in self.rotations:
                rot0 = numpy.dot(rot0, rot)
            return (trans0, rot0)
        return locals()

########################################
# Create Pixel x/y/z -> Cart x/y/z Map #
########################################

    def get_charge(self, i, j, k):
        """
        """
        if i > self.lx - 1:
            print "Warning, i index is too big ( > %d )" % (self.lx -1)
        if j > self.ly - 1:
            print "Warning, j index is too big ( > %d )" % (self.ly -1)
        if k > self.lz - 1:
            print "Warning, k index is too big ( > %d )" % (self.lz -1)

        return self.chargeArray[i + self.lx * j + self.lx * self.ly * k]

    def get_cart(self, i, j, k):
        """
        """
        if i > self.lx - 1:
            print "Warning, i index is too big ( > %d )" % (self.lx -1)
        if j > self.ly - 1:
            print "Warning, j index is too big ( > %d )" % (self.ly -1)
        if k > self.lz - 1:
            print "Warning, k index is too big ( > %d )" % (self.lz -1)

        return self.cartArray[i + self.lx * j + self.lx * self.ly * k]

    def init_cart(self):
        # Center and translate
        chargeArray = self.chargeArray.reshape((len(self.chargeArray), 1))
        tmp = (self.cartArray * chargeArray)
        coq = tmp.sum(axis=0)/chargeArray.sum()
        self.translations.append(coq)
        self._cartArray -= coq
        # Orient and rotate
        evec = self.get_chargeTensor(eigen=True)[1]
        self.rotate_byMatrix(evec.transpose())
        self.rotations.append(evec)

    def get_chargeTensor(self, eigen=False, array=None):
        tmp = pychm.emap.ext.get_Tensor(self._cartArray, self.chargeArray)
        #
        if eigen:
            return numpy.linalg.eigh(tmp)
        else:
            return tmp

    def rotate_byMatrix(self, matrix):
        """
        """
        self._cartArray = numpy.dot(self._cartArray, matrix.transpose())

    def build_map(self):
        # create coordinate matrix
        iterator = ( crd for atom in self.mol for crd in atom.cart )
        crdMatrix = numpy.fromiter(iterator, float)
        crdMatrix.resize((len(self.mol), 3))
        if self.transVector is not None:
            crdMatrix += self.transVector
        if self.rotMatrix is not None:
            crdMatrix = dot(crdMatrix, self.rotMatrix.transpose())

        # translate minimum to origin
        xMin = crdMatrix[:,0].min()
        yMin = crdMatrix[:,1].min()
        zMin = crdMatrix[:,2].min()
        crdMatrix -= numpy.array([xMin, yMin, zMin])

        # pad minimum borders
        padding = 2 * self.resolution / self.diff + 1
        crdMatrix += padding

        # pad maximum borders
        xMax = int(crdMatrix[:,0].max() + padding[0])
        yMax = int(crdMatrix[:,1].max() + padding[1])
        zMax = int(crdMatrix[:,2].max() + padding[2])

        # initialize
        tmp = numpy.zeros( xMax * yMax * zMax )
        tmp.resize((xMax, yMax, zMax))
        tmp1 = deepcopy(tmp)

        # build kernel
        for i, atom in enumerate(self.mol):
            # for each atom we pull the charge apart and assign it to one
            # of 8 verticies on a unit cube using a linear scheme
            #
            # c -- atom charge
            # x0, y0, z0 -- the lesser index of the cube vertex
            # x1, y1, z1 -- the greater index of the cube vertex
            # x, y, z -- the cartesian position inside the unit cube

            c = atom.charge
            x0, y0, z0 = map(int, crdMatrix[i] // self.resolution)
            x1, y1, z1 = (x0 + 1, y0 + 1, z0 + 1)
            x, y, z = crdMatrix[i] % self.resolution / self.resolution

            # populate the 8 cube verticies
            tmp[x0][y0][z0] += c * (1 - x) * (1 - y) * (1 - z)
            tmp[x1][y0][z0] += c * x * (1 - y) * (1 - z)
            tmp[x0][y1][z0] += c * (1 - x) * y * (1 - z)
            tmp[x0][y0][z1] += c * (1 - x) * (1 - y) * z
            tmp[x1][y1][z0] += c * x * y * (1 - z)
            tmp[x1][y0][z1] += c * x * (1 - y) * z
            tmp[x0][y1][z1] += c * (1 - x) * y * z
            tmp[x1][y1][z1] += c * x * y * z

        # blur kernel
        totalCharge = tmp.sum()
        filters.gaussian_filter(tmp, self.resolution/2., order=0, output=tmp1,
                                mode='constant', cval=0.)
        tmp1 *= totalCharge / tmp1.sum()

        # return
        return (tmp, tmp1)
