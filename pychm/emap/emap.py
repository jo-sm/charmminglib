"""
"""


import numpy as np
import scipy.ndimage.filters as filters
from copy import deepcopy
from pychm.tools import Property
import pychm.emap.ext as ext
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

    _defaults = {
        'resolution' : np.float64(15),
        'rhocut' : np.float64(0),
        'dx' : np.float64(3),
        'dy' : np.float64(3),
        'dz' : np.float64(3)
    }

    _metaData = [
        'nc', 'nr', 'ns',
        'mode',
        'ncstart', 'nrstart', 'nsstart',
        'nx', 'ny', 'nz',
        'xlen', 'ylen', 'zlen',
        'alpha', 'beta', 'gamma',
        'mapc', 'mapr', 'maps',
        'amin', 'amax', 'amean',
        #'ispg', 'nsymbt', 'lskflg', 'skwmat', 'skwtrn',
        #'extra', 'maplabel', 'macst', 'arms', 'nlabl', 'label'
        ]

    def __init__(self, mol=None):
        super(EMap, self).__init__()
        #
        self._cartArray = None
        self.chargeArray = None
        self.coreArray = None
        self.laplacianArray = None
        # Operations
        self.translations = []
        self.rotations = []
        # default parameters
        for k, v in self.__class__._defaults.iteritems():
            setattr(self, k, v)

    @Property
    def shape():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return (self.lx, self.ly, self.lz)
        return locals()

    @Property
    def res():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return (self.xres, self.yres, self.zres)
        return locals()

    @Property
    def xres():
        doc =\
        """
        """
        def fget(self):
            return self.xlen / self.nx
        def fset(self, value):
            self.xlen = self.nx * np.float64(value)
        return locals()

    @Property
    def yres():
        doc =\
        """
        """
        def fget(self):
            return self.ylen / self.ny
        def fset(self, value):
            self.ylen = self.ny * np.float64(value)
        return locals()

    @Property
    def zres():
        doc =\
        """
        """
        def fget(self):
            return self.zlen / self.nz
        def fset(self, value):
            self.zlen = self.nz * np.float64(value)
        return locals()

    @Property
    def cartArray3d():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return np.reshape(self.cartArray, (self.lx, self.ly, self.lz, 3), order='F')
        return locals()

    @Property
    def chargeArray3d():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return np.reshape(self.chargeArray, (self.lx, self.ly, self.lz), order='F')
        return locals()

    @Property
    def coreArray3d():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return np.reshape(self.coreArray, (self.lx, self.ly, self.lz), order='F')
        return locals()

    @Property
    def laplacianArray3d():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return np.reshape(self.laplacianArray, (self.lx, self.ly, self.lz), order='F')
        return locals()

    @Property
    def cartArray():
        doc =\
        """
        DOCME
        """
        def fget(self):
            trans, rot = self.state
            tmp = np.dot(self._cartArray, rot.transpose())
            return tmp + trans
        def fset(self, value):
            trans, rot = self.state
            value -= trans
            self._cartArray = np.dot(value, rot)
        return locals()

    @Property
    def state():
        doc =\
        """
        DOCME
        """
        def fget(self):
            trans0 = np.zeros(3)
            for trans in self.translations:
                trans0 += trans
            rot0 = np.identity(3)
            for rot in self.rotations:
                rot0 = np.dot(rot0, rot)
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
        tmp = ext.get_Tensor(self._cartArray, self.chargeArray)
        #
        if eigen:
            return np.linalg.eigh(tmp)
        else:
            return tmp

    def get_laplacianArray(self):
        return ext.get_laplacian(self.chargeArray3d)

    def get_coreArray(self):
        return ext.get_core(self.chargeArray3d, self.laplacianArray3d, self.rhocut)

    def rotate_byMatrix(self, matrix):
        """
        """
        self._cartArray = np.dot(self._cartArray, matrix.transpose())

    def densityCorrelation(self, other):
        """
        DC_{mn} = \frac{    \bar{\rho_m \rho_n} -\bar{\rho_m}\bar{\rho_n}}
                       {    \delta(\rho_m)\delta(\rho_n)}
        """
        assert isinstance(other, EMap)
        sRho = self.chargeArray
        oRho = other.chargeArray
        rhoSquared = np.fromiter(itertools.imap(lambda (x,y): x*y,
                                itertools.izip(sRho, oRho)), dtype=np.float64)
        return ((rhoSquared.mean() - sRho.mean() * oRho.mean()) /
                (sRho.std() * oRho.std()))

    def laplacianCorrelation(self, other):
        """
        """
        assert isinstance(other, EMap)
        sLap = self.laplacianArray
        oLap = other.laplacianArray
        lapSquared = np.fromiter(itertools.imap(lambda (x,y): x*y,
                                itertools.izip(sLap, oLap)), dtype=np.float64)
        return ((lapSquared.mean() - sLap.mean() * oLap.mean()) /
                (sLap.std() * oLap.std()))

    def pad_mapByPixel(self, lessI, moreI, lessJ, moreJ, lessK, moreK):
        lessI = int(lessI)
        moreI = int(moreI)
        lessJ = int(lessJ)
        moreJ = int(moreJ)
        lessK = int(lessK)
        moreK = int(moreK)
        assert lessI >= 0
        assert moreI >= 0
        assert lessJ >= 0
        assert moreJ >= 0
        assert lessK >= 0
        assert moreK >= 0
        #
        tmp = self._new()
        tmp.nxstart = self.nxstart - (lessI + moreI - 1)
        tmp.nystart = self.nystart - (lessJ + moreJ - 1)
        tmp.nzstart = self.nzstart - (lessK + moreK - 1)
        tmp.ncstart = tmp.nxstart
        tmp.nrstart = tmp.nystart
        tmp.nsstart = tmp.nzstart
        #
        chargeArray3d = ext.padZero(self.chargeArray3d, lessI, moreI, lessJ, moreJ, lessK, moreK)
        tmp.lx, tmp.ly, tmp.lz = chargeArray3d.shape
        tmp.nc, tmp.nr, tmp.ns = chargeArray3d.shape
        tmp.chargeArray = chargeArray3d.flatten('F')
        #
        tmp.cartArray = ext.build_cart(tmp.lx, tmp.ly, tmp.lz,
                                        tmp.nxstart, tmp.nystart, tmp.nzstart,
                                        tmp.xres, tmp.yres, tmp.zres)
        #
        return tmp

    def pad_mapByCart(self, lessX, moreX, lessY, moreY, lessZ, moreZ):
        assert lessX >= 0.
        assert moreX >= 0.
        assert lessY >= 0.
        assert moreY >= 0.
        assert lessZ >= 0.
        assert moreZ >= 0.
        #
        lessI = int(lessX / self.xres) + 1
        moreI = int(moreX / self.xres) + 1
        lessJ = int(lessY / self.yres) + 1
        moreJ = int(moreY / self.yres) + 1
        lessK = int(lessZ / self.zres) + 1
        moreK = int(moreZ / self.zres) + 1
        #
        return self.pad_mapByPixel(lessI, moreI, lessJ, moreJ, lessK, moreK)

    def combine_maps(self, other):
        pass

    def recast_map(self):
        pass

    def write_map(self):
        pass

    def _new(self):
        tmp = EMap()
        tmp.nx = self.nx
        tmp.ny = self.ny
        tmp.nz = self.nz
        tmp.mode = 2
        tmp.xlen = self.xlen
        tmp.ylen = self.ylen
        tmp.zlen = self.zlen
        tmp.alpha = self.alpha
        tmp.beta = self.beta
        tmp.gamma = self.gamma
        tmp.mapc = 1
        tmp.mapr = 2
        tmp.maps = 3
        tmp.rotations = deepcopy(self.rotations)
        tmp.translations = deepcopy(self.translations)
        return tmp

    def build_map(self):
        # create coordinate matrix
        iterator = ( crd for atom in self.mol for crd in atom.cart )
        crdMatrix = np.fromiter(iterator, float)
        crdMatrix.resize((len(self.mol), 3))
        if self.transVector is not None:
            crdMatrix += self.transVector
        if self.rotMatrix is not None:
            crdMatrix = dot(crdMatrix, self.rotMatrix.transpose())

        # translate minimum to origin
        xMin = crdMatrix[:,0].min()
        yMin = crdMatrix[:,1].min()
        zMin = crdMatrix[:,2].min()
        crdMatrix -= np.array([xMin, yMin, zMin])

        # pad minimum borders
        padding = 2 * self.resolution / self.diff + 1
        crdMatrix += padding

        # pad maximum borders
        xMax = int(crdMatrix[:,0].max() + padding[0])
        yMax = int(crdMatrix[:,1].max() + padding[1])
        zMax = int(crdMatrix[:,2].max() + padding[2])

        # initialize
        tmp = np.zeros( xMax * yMax * zMax )
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

