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
        self.pixelArray = None
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

        tmp = numpy.zeros((self.lx, self.ly, self.lz, 3))
        iterable = ( (i, j, k) for i in xrange(self.lx) for j in xrange(self.ly)
                    for k in xrange(self.lz) )
        for i, j, k in iterable:
            tmp[i][j][k] = pixel2cart_map(i, j, k)
        #
        self.cartArray = tmp

    def get_coq(self):
        """
        """
        iterable = ( (i, j, k) for i in xrange(self.lx) for j in xrange(self.ly)
                    for k in xrange(self.lz) )
        tmp = numpy.zeros(3)
        for i, j, k in iterable:
            tmp += self.cartArray[i][j][k] * self.pixelArray[i][j][k]
        return tmp / self.pixelArray.sum()

    def get_basis(self):
        """
        """
        # Centered Cartesian Array
        tmp = self.cartArray - self.get_coq()
        # Compute Moment of Charge Tensor
        xx, yy, zz, xy, xz, yz = numpy.zeros(6)
        iterable = ( (i, j, k) for i in xrange(self.lx) for j in xrange(self.ly)
                    for k in xrange(self.lz) )
        for i, j, k in iterable:
            x, y, z = tmp[i][j][k]
            q = self.pixelArray[i][j][k]
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
        self.basis = numpy.linalg.eig(tmp)[1].transpose()

    def rotate_byMatrix(self, matrix):
        """
        """
        # Unpack cartesian data into n x 3 matrix
        tmp = numpy.reshape(self.cartArray, (self.lx * self.ly * self.lz, 3))
        # Rotate
        tmp = numpy.dot(tmp, matrix.transpose())
        # Repack
        tmp.resize((self.lx, self.ly, self.lz, 3))
        return tmp

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
