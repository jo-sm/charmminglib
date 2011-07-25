"""
"""


import numpy
import scipy.ndimage.filters as filters
from copy import deepcopy


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
        self.cartArray = None
        self.chargeArray = None
        ## Parameters
        #self.resolution = 15. # In Angstrom
        ## State Variables
        #self.transVector = None
        #self.rotMatrix = None
        #if arg is not None:
        #    self.mol = mol
        #    self.mol.orient()

########################################
# Create Pixel x/y/z -> Cart x/y/z Map #
########################################

    def set_cartArrayFromPixelArray(self):
        """
        """

        def pixel2cart_map(x, y, z):
            """
            I think this is wrong because nc/nr/ns-start are not properly permuted.
            """
            xcart = (self.ncstart + x) * (self.xlen/self.nx)
            ycart = (self.nrstart + y) * (self.ylen/self.ny)
            zcart = (self.nsstart + z) * (self.zlen/self.nz)
            return numpy.array([xcart, ycart, zcart])

        def fixed_pixel2cart_map(x, y, z):
            """
            Properly permuted starting indicies.

            TODO -- Include crystal skews in mapping.
            """
            # Do work
            xcart = (self.nxstart + x) * (self.xlen/self.nx)
            ycart = (self.nystart + y) * (self.ylen/self.ny)
            zcart = (self.nzstart + z) * (self.zlen/self.nz)
            return numpy.array([xcart, ycart, zcart])

        # make 3-D array
        tmp = numpy.zeros((self.lx, self.ly, self.lz, 3))
        iterable = ( (i, j, k) for i in xrange(self.lx) for j in xrange(self.ly)
                    for k in xrange(self.lz) )
        for i, j, k in iterable:
            tmp[i][j][k] = pixel2cart_map(i, j, k)
        self.cartArray = tmp
        # store flattened cart array
        self.flatCartArray = self.get_flatCartArray()

    def get_flatChargeArray(self):
        """
        """
        iterable = ( self.chargeArray[i][j][k] for k in xrange(self.lz)
                    for j in xrange(self.ly) for i in xrange(self.lx) )
        tmp = numpy.fromiter(iterable, numpy.float)
        tmp.resize(len(tmp), 1)
        return tmp

    def get_flatCartArray(self):
        """
        """
        tmp = numpy.zeros((self.lx * self.ly * self.lz, 3))
        iterable = ( self.cartArray[i][j][k] for k in xrange(self.lz)
                    for j in xrange(self.ly) for i in xrange(self.lx) )
        for n in xrange(self.lx * self.ly * self.lz):
            tmp[n] = iterable.next()
        return tmp

    def get_coq(self):
        return (self.flatCartArray * self.flatChargeArray).sum(axis=0)/self.chargeArray.sum()

    def get_basis(self):
        tmpCart = self.flatCartArray - self.get_coq()
        tmpCharge = self.flatChargeArray.flatten()
        xx, yy, zz, xy, xz, yz = numpy.zeros(6)
        for i, cart in enumerate(tmpCart):
            x, y, z = cart
            q = tmpCharge[i]
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
        return numpy.linalg.eig(tmp)[1]

    def rotate_byMatrix(self, matrix):
        """
        """
        # Centered Cartesian Array
        tmp = self.cartArray - self.get_coq()
        # Unpack cartesian data into n x 3 matrix
        tmp.resize((self.lx * self.ly * self.lz, 3))
        # Rotate
        tmp = numpy.dot(tmp, matrix)
        # Repack
        tmp.resize((self.lx, self.ly, self.lz, 3))
        return tmp

    def rotate_2(self, matrix):
        """
        """
        tmp = deepcopy(self.flatCartArray)
        tmp -= self.get_coq()
        tmp = numpy.dot(tmp, matrix)
        tmp += self.get_coq()
        self.intermediate = tmp
        returnMe = numpy.zeros(self.lx * self.ly * self.lz * 3)
        returnMe.resize((self.lx, self.ly, self.lz, 3))
        iterator = ( (i, j, k) for k in xrange(self.lz) for j in xrange(self.ly)
                    for i in xrange(self.lx) )
        for n, (i, j, k) in enumerate(iterator):
            returnMe[i][j][k] = tmp[n]
        return returnMe

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
