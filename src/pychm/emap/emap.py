"""
"""


import numpy
import scipy.ndimage.filters as filters
from copy import deepcopy
from pychm.tools import Property
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
        # Operations
        self.translations = []
        self.rotations = []
        ## Parameters
        #self.resolution = 15. # In Angstrom
        ## State Variables
        #self.transVector = None
        #self.rotMatrix = None
        #if arg is not None:
        #    self.mol = mol
        #    self.mol.orient()

    @Property
    def cartArray():
        doc =\
        """
        DOCME
        """
        def fget(self):
            tmpTrans = numpy.zeros(3)
            for trans in self.translations:
                tmp += trans
            return self._cartArray + tmp
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

    def init_cart2(self):
        # Orient and rotate
        tensor0 = self.get_chargeTensor(eigen=True)[1]
        self.rotate_byMatrix(tensor0.transpose())
        tensor1 = self.get_chargeTensor(eigen=True)[1]
        self.rotate_byMatrix(tensor1.transpose())
        self.rotations.append(numpy.dot(tensor0, tensor1))

    def get_chargeTensor(self, eigen=False, array=None):
        if array is None:
            array = self._cartArray
        xx, yy, zz, xy, xz, yz = (0., 0., 0., 0., 0., 0.)
        for i, cart in enumerate(array):
            x, y, z = cart
            q = self.chargeArray[i]
            xx += q*(y*y+z*z)
            yy += q*(x*x+z*z)
            zz += q*(x*x+y*y)
            xy += q*x*y
            xz += q*x*z
            yz += q*y*z
        #
        tmp = numpy.array([
            [ xx, -yz, -xz],
            [-yz,  yy, -yz],
            [-xz, -yz,  zz]
            ])
        #
        if eigen:
            return numpy.linalg.eigh(tmp)
        else:
            return tmp

    def rotate_byMatrix(self, matrix):
        """
        """
        self._cartArray = numpy.dot(self._cartArray, matrix.transpose())

    def orient(self):
        self.rotate_byMatrix(self.get_chargeTensor(eigen=True)[1].transpose())

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

    def get_chargeTensorAddTwo(self, eigen=False):
        xx = []
        yy = []
        zz = []
        xy = []
        xz = []
        yz = []
        for i, cart in enumerate(self._cartArray):
            x, y, z = cart
            q = self.chargeArray[i]
            xx.append(q*(y*y+z*z))
            yy.append(q*(x*x+z*z))
            zz.append(q*(x*x+y*y))
            xy.append(q*x*y)
            xz.append(q*x*z)
            yz.append(q*y*z)
        xx = pairwise_add(xx)
        yy = pairwise_add(yy)
        zz = pairwise_add(zz)
        xy = pairwise_add(xy)
        xz = pairwise_add(xz)
        yz = pairwise_add(yz)
        #
        tmp = numpy.array([
            [ xx, -yz, -xz],
            [-yz,  yy, -yz],
            [-xz, -yz,  zz]
            ])
        #
        if eigen:
            return numpy.linalg.eigh(tmp)
        else:
            return tmp
    #def set_cartArrayFromPixelArray(self):
    #    """
    #    """

    #    def pixel2cart_map(x, y, z):
    #        """
    #        I think this is wrong because nc/nr/ns-start are not properly permuted.
    #        """
    #        xcart = (self.ncstart + x) * (self.xlen/self.nx)
    #        ycart = (self.nrstart + y) * (self.ylen/self.ny)
    #        zcart = (self.nsstart + z) * (self.zlen/self.nz)
    #        return numpy.array([xcart, ycart, zcart])

    #    def fixed_pixel2cart_map(x, y, z):
    #        """
    #        Properly permuted starting indicies.

    #        TODO -- Include crystal skews in mapping.
    #        """
    #        # Do work
    #        xcart = (self.nxstart + x) * (self.xlen/self.nx)
    #        ycart = (self.nystart + y) * (self.ylen/self.ny)
    #        zcart = (self.nzstart + z) * (self.zlen/self.nz)
    #        return numpy.array([xcart, ycart, zcart])

    #    # make 3-D array
    #    tmp = numpy.zeros((self.lx, self.ly, self.lz, 3))
    #    iterable = ( (i, j, k) for i in xrange(self.lx) for j in xrange(self.ly)
    #                for k in xrange(self.lz) )
    #    for i, j, k in iterable:
    #        tmp[i][j][k] = pixel2cart_map(i, j, k)
    #    self.cartArray = tmp
    #    # store flattened cart array
    #    self.flatCartArray = self.get_flatCartArray()

