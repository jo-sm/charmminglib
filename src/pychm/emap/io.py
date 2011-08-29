"""
"""


import struct
import numpy as np
from emap import EMap


class BaseBinaryIO(object):
    """
    """
    def read(self):
        NotImplementedError

    def write(self):
        NotImplementedError


class CCP4File(BaseBinaryIO):
    """
    ccp4 file specification as per --
    http://www.ccp4.ac.uk/html/maplib.html#description

    # HEADER VARS

    My Var                  WXW Var                 Comment
    ----------------        ----------------        ------------------------
    nc, nr, ns              lx, ly, lz              # of Col/Row/Sect in array
    mode                    mode                    should be `2` -> Real
    ncstart, nrstart, ...   nsx, nsy, nsz           number of first col in map
    nx, ny, nz              mx, my, mz              number of interval in unit cell
    xlen, ylen, zlen        xlen, ylen, zlen        unit cell dim, in Ang
    alpha, beta, gamma      alpha, beta, gamma      unit cell angles, in Deg
    mapc, mapr, maps        mapc, mapr, maps        col/row/sec -> x/y/z mapping
    amin, amax, amean       amin, amax, amean       eDensity statistics
    ispg                    ispg                    space group number
    nsymbt                  nsymbt                  number of bytes for storing sym op
    lskflg                  lskflg                  flag for skew transformation (=0 none,
                                                    =1 if it follows in `symmetry`
    skwmat                  skwmat                  skew matrix S (in order S11, S12, S13,
                                                    S21, etc) if `lskflg` != 0
    skwtrn                  skwtrn                  skew translation t if `lskflg` != 0
    extra                   extra                   UNUSED
    maplabel                maplabel                the string 'MAP '
    machst                  machst                  machine stamp - indicates what wrote file
    arms                    arms                    rms dev of map from mean density
    nlabl                                           # of labels being used
    label                                           a list of 10 labels, 80 characters each
    symmetry                                        TODO

    # NON HEADER VARS
    chargeArray              A `numpy.array` of dimension lx, ly, lz, whose values are
                            electron density in coulomb
    transposeTuple          a tuple of mapc, mapr, maps -- only zero indexed
    lx, ly, lz              # of x, y, z in array -- it has been mapped from nc, nr, ns
                            and correctly permutated using mapc, mapr, maps
    nxstart, nystart, ...   Starting pixel index in the x, y, z directions -- it has been
                            mapped from ncstart, nrstart, nsstart and correctly permuted
    """
    @classmethod
    def read(cls, filename):
        cls.map = EMap()
        cls.filename = filename
        cls._read_binData()
        cls._map_pixelVars2CartVars()
        cls._processDensityData()
        cls._processCartData()
        cls.map.init_cart()
        return cls.map

    @classmethod
    def _read_binData(cls):
        with open(cls.filename, mode='rb') as fp:
            nc, nr, ns = struct.unpack('iii', fp.read(12))
            mode = struct.unpack('i', fp.read(4))[0]
            ncstart, nrstart, nsstart = struct.unpack('iii', fp.read(12))
            nx, ny, nz = map(np.float64, struct.unpack('iii', fp.read(12)))
            xlen, ylen, zlen = map(np.float64, struct.unpack('fff', fp.read(12)))
            alpha, beta, gamma = map(np.radians, struct.unpack('fff', fp.read(12)))
            mapc, mapr, maps = struct.unpack('iii', fp.read(12))
            amin, amax, amean = map(np.float64, struct.unpack('fff', fp.read(12)))
            ispg = struct.unpack('i', fp.read(4))[0]
            nsymbt = struct.unpack('i', fp.read(4))[0]
            lskflg = struct.unpack('i', fp.read(4))[0]
            skwmat = struct.unpack('f'*9, fp.read(36))
            skwtrn = map(np.float64, struct.unpack('fff', fp.read(12)))
            extra = struct.unpack('i'*15, fp.read(60))
            maplabel = ''.join(struct.unpack('cccc', fp.read(4)))
            machst = struct.unpack('cccc', fp.read(4))
            arms = np.float64(struct.unpack('f', fp.read(4))[0])
            nlabl = struct.unpack('i', fp.read(4))[0]
            label = []
            for i in range(10):
                label.append(''.join(struct.unpack('c'*80, fp.read(80))))

        # TODO -- This needs a testcase to debug it!
            symmetry = []
            for i in range(nsymbt/80):
                symmetry.append(''.join(struct.unpack('c'*80, fp.read(80))))
        # TODO

            rawArray = np.array(struct.unpack('f'*nc*nr*ns,
                                            fp.read(4*nc*nr*ns)), dtype=np.float64)
        cls.map.__dict__.update(locals())
        del cls.map.fp
        del cls.map.cls
        del cls.map.i

    @classmethod
    def _map_pixelVars2CartVars(cls):
        tmp = cls.map

        if tmp.mapc == 1 and tmp.mapr == 2 and tmp.maps == 3:
            tmp.lx = tmp.nc
            tmp.ly = tmp.nr
            tmp.lz = tmp.ns
            tmp.nxstart = tmp.ncstart
            tmp.nystart = tmp.nrstart
            tmp.nzstart = tmp.nsstart
            tmp.transposeTuple = (0, 1, 2)
        elif tmp.mapc == 1 and tmp.mapr == 3 and tmp.maps == 2:
            tmp.lx = tmp.nc
            tmp.ly = tmp.ns
            tmp.lz = tmp.nr
            tmp.nxstart = tmp.ncstart
            tmp.nystart = tmp.nsstart
            tmp.nzstart = tmp.nrstart
            tmp.transposeTuple = (0, 2, 1)
        elif tmp.mapc == 2 and tmp.mapr == 1 and tmp.maps == 3:
            tmp.lx = tmp.nr
            tmp.ly = tmp.nc
            tmp.lz = tmp.ns
            tmp.nxstart = tmp.nrstart
            tmp.nystart = tmp.ncstart
            tmp.nzstart = tmp.nsstart
            tmp.transposeTuple = (1, 0, 2)
        elif tmp.mapc == 2 and tmp.mapr == 3 and tmp.maps == 1:
            tmp.lx = tmp.ns
            tmp.ly = tmp.nc
            tmp.lz = tmp.nr
            tmp.nxstart = tmp.nsstart
            tmp.nystart = tmp.ncstart
            tmp.nzstart = tmp.nrstart
            tmp.transposeTuple = (1, 2, 0)
        elif tmp.mapc == 3 and tmp.mapr == 1 and tmp.maps == 2:
            tmp.lx = tmp.nr
            tmp.ly = tmp.ns
            tmp.lz = tmp.nc
            tmp.nxstart = tmp.nrstart
            tmp.nystart = tmp.nsstart
            tmp.nzstart = tmp.ncstart
            tmp.transposeTuple = (2, 0, 1)
        elif tmp.mapc == 3 and tmp.mapr == 2 and tmp.maps == 1:
            tmp.lx = tmp.ns
            tmp.ly = tmp.nr
            tmp.lz = tmp.nc
            tmp.nxstart = tmp.nsstart
            tmp.nystart = tmp.nrstart
            tmp.nzstart = tmp.ncstart
            tmp.transposeTuple = (2, 1, 0)
        else:
            raise IOError("NULL: Error with Axis Remapping Data.")

    @classmethod
    def _processDensityData(cls):
        tmp = cls.map
        rawArray = tmp.rawArray

        # Roll raw electron density array data into a 3D matrix
        rawMatrix = np.zeros((tmp.nc, tmp.nr, tmp.ns))
        n = 0
        for k in xrange(tmp.ns):
            for j in xrange(tmp.nr):
                for i in xrange(tmp.nc):
                    rawMatrix[i, j, k] = rawArray[n]
                    n += 1

        # Remap electron density data
        chargeMatrix = np.transpose(rawMatrix, tmp.transposeTuple)

        # Unroll electron density data into flat array
        iterator = ( chargeMatrix[i, j, k] for k in xrange(tmp.lz)
                    for j in xrange(tmp.ly) for i in xrange(tmp.lx) )
        tmp.chargeArray = np.fromiter(iterator, np.float64)

    @classmethod
    def _processCartData(cls):
        tmp = cls.map
        # build cartesian data
        tmp._cartArray = np.zeros((tmp.lx * tmp.ly * tmp.lz, 3))
        n = 0
        for k in xrange(tmp.nzstart, tmp.nzstart + tmp.lz):
            for j in xrange(tmp.nystart, tmp.nystart + tmp.ly):
                for i in xrange(tmp.nxstart, tmp.nxstart + tmp.lx):
                    tmp._cartArray[n, 0] = i
                    tmp._cartArray[n, 1] = j
                    tmp._cartArray[n, 2] = k
                    n += 1
        # scale
        tmp._cartArray *= np.array([tmp.xlen/tmp.nx, tmp.ylen/tmp.ny,
                                    tmp.zlen/tmp.nz])
        if tmp.alpha != np.pi/2 or tmp.beta != np.pi/2 or tmp.gamma != np.pi/2:
            # skew matrix (un-tested)
            cg = np.cos(tmp.gamma)
            sg = np.sin(tmp.gamma)
            cb = np.cos(tmp.beta)
            sb = np.sin(tmp.beta)
            ca = np.cos(tmp.alpha)
            sa = np.sin(tmp.alpha)
            skew = np.array([(1, cg,      cb),
                             (0, sg, sb * ca),
                             (0,  0, sb * sa)])
            #
            tmp._cartArray = np.dot(tmp._cartArray, skew.T)


if __name__ == '__main__':


    import sys

    print "\nThis currently only supports interactive debugging.\n"
    print "`taco` = `EMap()`"

    #taco = read_CCP4File(sys.argv[1])
    taco = CCP4File.read(sys.argv[1])
    #taco = read_CCP4FileWXW(sys.argv[1])

    #rawArray2 = numpy.zeros((taco.nc, taco.nr, taco.ns))
    #for n, ele in enumerate(taco.rawArray):
    #    k = n // ( taco.nc * taco.nr)
    #    j = (n % ( taco.nc * taco.nr )) // taco.nc
    #    i = ((n % ( taco.nc * taco.nr )) % taco.nc)

    #    k1 = n // (taco.nc * taco.nr)
    #    q1 = n - k1 * taco.nc * taco.nr
    #    j1 = q1 // taco.nc
    #    i1 = q1 - j1 * taco.nc

    #    if not k == k1 and not j == j1 and not i == i1:
    #        print 'ASDF'
    #        break
    #def reMap(n):
    #    k = n // ( taco.nc * taco.nr)
    #    j = (n % ( taco.nc * taco.nr )) // taco.nc
    #    i = ((n % ( taco.nc * taco.nr )) % taco.nc)
    #    return (i, j, k)


def read_CCP4File(filename):
    """
    ccp4 file specification as per --
    http://www.ccp4.ac.uk/html/maplib.html#description

    # HEADER VARS

    My Var                  WXW Var                 Comment
    ----------------        ----------------        ------------------------
    nc, nr, ns              lx, ly, lz              # of Col/Row/Sect in array
    mode                    mode                    should be `2` -> Real
    ncstart, nrstart, ...   nsx, nsy, nsz           number of first col in map
    nx, ny, nz              mx, my, mz              number of interval in unit cell
    xlen, ylen, zlen        xlen, ylen, zlen        unit cell dim, in Ang
    alpha, beta, gamma      alpha, beta, gamma      unit cell angles, in Deg
    mapc, mapr, maps        mapc, mapr, maps        col/row/sec -> x/y/z mapping
    amin, amax, amean       amin, amax, amean       eDensity statistics
    ispg                    ispg                    space group number
    nsymbt                  nsymbt                  number of bytes for storing sym op
    lskflg                  lskflg                  flag for skew transformation (=0 none,
                                                    =1 if it follows in `symmetry`
    skwmat                  skwmat                  skew matrix S (in order S11, S12, S13,
                                                    S21, etc) if `lskflg` != 0
    skwtrn                  skwtrn                  skew translation t if `lskflg` != 0
    extra                   extra                   UNUSED
    maplabel                maplabel                the string 'MAP '
    machst                  machst                  machine stamp - indicates what wrote file
    arms                    arms                    rms dev of map from mean density
    nlabl                                           # of labels being used
    label                                           a list of 10 labels, 80 characters each
    symmetry                                        TODO

    # NON HEADER VARS
    chargeArray              A `numpy.array` of dimension lx, ly, lz, whose values are
                            electron density in coulomb
    transposeTuple          a tuple of mapc, mapr, maps -- only zero indexed
    lx, ly, lz              # of x, y, z in array -- it has been mapped from nc, nr, ns
                            and correctly permutated using mapc, mapr, maps
    nxstart, nystart, ...   Starting pixel index in the x, y, z directions -- it has been
                            mapped from ncstart, nrstart, nsstart and correctly permuted
    """
    import numpy
    # Instantize!
    tmp = EMap()

    with open(filename, mode='rb') as fpointer:
        tmp.nc, tmp.nr, tmp.ns = struct.unpack('iii', fpointer.read(12))
        tmp.mode = struct.unpack('i', fpointer.read(4))[0]
        tmp.ncstart, tmp.nrstart, tmp.nsstart = struct.unpack('iii', fpointer.read(12))
        tmp.nx, tmp.ny, tmp.nz = struct.unpack('iii', fpointer.read(12))
        tmp.xlen, tmp.ylen, tmp.zlen = struct.unpack('fff', fpointer.read(12))
        tmp.alpha, tmp.beta, tmp.gamma = struct.unpack('fff', fpointer.read(12))
        tmp.mapc, tmp.mapr, tmp.maps = struct.unpack('iii', fpointer.read(12))
        tmp.amin, tmp.amax, tmp.amean = struct.unpack('fff', fpointer.read(12))
        tmp.ispg = struct.unpack('i', fpointer.read(4))[0]
        tmp.nsymbt = struct.unpack('i', fpointer.read(4))[0]
        tmp.lskflg = struct.unpack('i', fpointer.read(4))[0]
        tmp.skwmat = struct.unpack('f'*9, fpointer.read(36))
        tmp.skwtrn = struct.unpack('fff', fpointer.read(12))
        tmp.extra = struct.unpack('i'*15, fpointer.read(60))
        tmp.maplabel = ''.join(struct.unpack('cccc', fpointer.read(4)))
        tmp.machst = struct.unpack('cccc', fpointer.read(4))
        tmp.arms = struct.unpack('f', fpointer.read(4))[0]
        tmp.nlabl = struct.unpack('i', fpointer.read(4))[0]
        tmp.label = []
        for i in range(10):
            tmp.label.append(''.join(struct.unpack('c'*80, fpointer.read(80))))


    # TODO -- This needs a testcase to debug it!
        tmp.symmetry = []
        for i in range(tmp.nsymbt/80):
            tmp.symmetry.append(''.join(struct.unpack('c'*80, fpointer.read(80))))
    # TODO


        rawArray = numpy.array(struct.unpack('f'*tmp.nc*tmp.nr*tmp.ns,
                                        fpointer.read(4*tmp.nc*tmp.nr*tmp.ns)))

    ##
    ## done with file reading
    ##

    # Roll raw electron density array data into a 3D matrix
    rawMatrix = numpy.zeros((tmp.nc, tmp.nr, tmp.ns))
    n = 0
    for k in xrange(tmp.ns):
        for j in xrange(tmp.nr):
            for i in xrange(tmp.nc):
                rawMatrix[i, j, k] = rawArray[n]
                n += 1

    # Map from column/row/section space to x/y/z space
    if tmp.mapc == 1 and tmp.mapr == 2 and tmp.maps == 3:
        tmp.lx = tmp.nc
        tmp.ly = tmp.nr
        tmp.lz = tmp.ns
        tmp.nxstart = tmp.ncstart
        tmp.nystart = tmp.nrstart
        tmp.nzstart = tmp.nsstart
        transposeTuple = (0, 1, 2)
    elif tmp.mapc == 1 and tmp.mapr == 3 and tmp.maps == 2:
        tmp.lx = tmp.nc
        tmp.ly = tmp.ns
        tmp.lz = tmp.nr
        tmp.nxstart = tmp.ncstart
        tmp.nystart = tmp.nsstart
        tmp.nzstart = tmp.nrstart
        transposeTuple = (0, 2, 1)
    elif tmp.mapc == 2 and tmp.mapr == 1 and tmp.maps == 3:
        tmp.lx = tmp.nr
        tmp.ly = tmp.nc
        tmp.lz = tmp.ns
        tmp.nxstart = tmp.nrstart
        tmp.nystart = tmp.ncstart
        tmp.nzstart = tmp.nsstart
        transposeTuple = (1, 0, 2)
    elif tmp.mapc == 2 and tmp.mapr == 3 and tmp.maps == 1:
        tmp.lx = tmp.ns
        tmp.ly = tmp.nc
        tmp.lz = tmp.nr
        tmp.nxstart = tmp.nsstart
        tmp.nystart = tmp.ncstart
        tmp.nzstart = tmp.nrstart
        transposeTuple = (1, 2, 0)
    elif tmp.mapc == 3 and tmp.mapr == 1 and tmp.maps == 2:
        tmp.lx = tmp.nr
        tmp.ly = tmp.ns
        tmp.lz = tmp.nc
        tmp.nxstart = tmp.nrstart
        tmp.nystart = tmp.nsstart
        tmp.nzstart = tmp.ncstart
        transposeTuple = (2, 0, 1)
    elif tmp.mapc == 3 and tmp.mapr == 2 and tmp.maps == 1:
        tmp.lx = tmp.ns
        tmp.ly = tmp.nr
        tmp.lz = tmp.nc
        tmp.nxstart = tmp.nsstart
        tmp.nystart = tmp.nrstart
        tmp.nzstart = tmp.ncstart
        transposeTuple = (2, 1, 0)
    else:
        raise IOError("NULL: Error with Axis Remapping Data.")

    # Remap electron density data
    chargeMatrix = numpy.transpose(rawMatrix, transposeTuple)

    # Unroll electron density data into flat array
    tmp.chargeArray = numpy.zeros(tmp.lx * tmp.ly * tmp.lz)
    n = 0
    for k in xrange(tmp.lz):
        for j in xrange(tmp.ly):
            for i in xrange(tmp.lx):
                tmp.chargeArray[n] = chargeMatrix[i, j, k]
                n += 1

    # Generate cartesian data into 'flat' array
    tmp._cartArray = numpy.zeros((tmp.lx * tmp.ly * tmp.lz, 3))
    n = 0
    mz = tmp.zlen / tmp.nz
    my = tmp.ylen / tmp.ny
    mx = tmp.xlen / tmp.nx
    for k in xrange(tmp.nzstart, tmp.nzstart + tmp.lz):
        zcart = mz * k
        for j in xrange(tmp.nystart, tmp.nystart + tmp.ly):
            ycart = my * j
            for i in xrange(tmp.nxstart, tmp.nxstart + tmp.lx):
                xcart = mx * i
                tmp._cartArray[n, 0] = xcart
                tmp._cartArray[n, 1] = ycart
                tmp._cartArray[n, 2] = zcart
                n += 1

    # Fin!
    tmp._cartArray = tmp._cartArray
    tmp.chargeArray = tmp.chargeArray
    #tmp.init_cart()
    return tmp


def read_CCP4FileWXW(filename):
    # Instantize!
    tmp = EMap()

    with open(filename, mode='rb') as fpointer:
        tmp.nc, tmp.nr, tmp.ns = struct.unpack('iii', fpointer.read(12))
        tmp.mode = struct.unpack('i', fpointer.read(4))[0]
        tmp.ncstart, tmp.nrstart, tmp.nsstart = struct.unpack('iii', fpointer.read(12))
        tmp.nx, tmp.ny, tmp.nz = struct.unpack('iii', fpointer.read(12))
        tmp.xlen, tmp.ylen, tmp.zlen = struct.unpack('fff', fpointer.read(12))
        tmp.alpha, tmp.beta, tmp.gamma = struct.unpack('fff', fpointer.read(12))
        tmp.mapc, tmp.mapr, tmp.maps = struct.unpack('iii', fpointer.read(12))
        tmp.amin, tmp.amax, tmp.amean = struct.unpack('fff', fpointer.read(12))
        tmp.ispg = struct.unpack('i', fpointer.read(4))[0]
        tmp.nsymbt = struct.unpack('i', fpointer.read(4))[0]
        tmp.lskflg = struct.unpack('i', fpointer.read(4))[0]
        tmp.skwmat = struct.unpack('f'*9, fpointer.read(36))
        tmp.skwtrn = struct.unpack('fff', fpointer.read(12))
        tmp.extra = struct.unpack('i'*15, fpointer.read(60))
        tmp.maplabel = ''.join(struct.unpack('cccc', fpointer.read(4)))
        tmp.machst = struct.unpack('cccc', fpointer.read(4))
        tmp.arms = struct.unpack('f', fpointer.read(4))[0]
        tmp.nlabl = struct.unpack('i', fpointer.read(4))[0]
        tmp.label = []
        for i in range(10):
            tmp.label.append(''.join(struct.unpack('c'*80, fpointer.read(80))))


    # TODO -- This needs a testcase to debug it!
        tmp.symmetry = []
        for i in range(tmp.nsymbt/80):
            tmp.symmetry.append(''.join(struct.unpack('c'*80, fpointer.read(80))))
    # TODO


        rawArray = numpy.array(struct.unpack('f'*tmp.nc*tmp.nr*tmp.ns,
                                        fpointer.read(4*tmp.nc*tmp.nr*tmp.ns)))

    ##
    ## done with file reading
    ##

    # Map from column/row/section space to x/y/z space
    if tmp.mapc == 1 and tmp.mapr == 2 and tmp.maps == 3:
        tmp.lx = tmp.nc
        tmp.ly = tmp.nr
        tmp.lz = tmp.ns
        tmp.nxstart = tmp.ncstart
        tmp.nystart = tmp.nrstart
        tmp.nzstart = tmp.nsstart
        def transpose(i, j, k):
            return (i, j, k)
        def reMap(i, j, k):
            return i + j * tmp.nc + k * tmp.nc * tmp.nr
    elif tmp.mapc == 1 and tmp.mapr == 3 and tmp.maps == 2:
        tmp.lx = tmp.nc
        tmp.ly = tmp.ns
        tmp.lz = tmp.nr
        tmp.nxstart = tmp.ncstart
        tmp.nystart = tmp.nsstart
        tmp.nzstart = tmp.nrstart
        def transpose(i, j, k):
            return (i, k, j)
        def reMap(i, j, k):
            return i + k * tmp.nc + j * tmp.nc * tmp.ns
    elif tmp.mapc == 2 and tmp.mapr == 1 and tmp.maps == 3:
        tmp.lx = tmp.nr
        tmp.ly = tmp.nc
        tmp.lz = tmp.ns
        tmp.nxstart = tmp.nrstart
        tmp.nystart = tmp.ncstart
        tmp.nzstart = tmp.nsstart
        def transpose(i, j, k):
            return (j, i, k)
        def reMap(i, j, k):
            return j + i * tmp.nr + k * tmp.nr * tmp.nc
    elif tmp.mapc == 2 and tmp.mapr == 3 and tmp.maps == 1:
        tmp.lx = tmp.ns
        tmp.ly = tmp.nc
        tmp.lz = tmp.nr
        tmp.nxstart = tmp.nsstart
        tmp.nystart = tmp.ncstart
        tmp.nzstart = tmp.nrstart
        def transpose(i, j, k):
            return (j, k, i)
        def reMap(i, j, k):
            return j + k * tmp.nr + i * tmp.nr * tmp.ns
    elif tmp.mapc == 3 and tmp.mapr == 1 and tmp.maps == 2:
        tmp.lx = tmp.nr
        tmp.ly = tmp.ns
        tmp.lz = tmp.nc
        tmp.nxstart = tmp.nrstart
        tmp.nystart = tmp.nsstart
        tmp.nzstart = tmp.ncstart
        def transpose(i, j, k):
            return (k, i, j)
        def reMap(i, j, k):
            return k + i * tmp.ns + j * tmp.ns * tmp.nc
    elif tmp.mapc == 3 and tmp.mapr == 2 and tmp.maps == 1:
        tmp.lx = tmp.ns
        tmp.ly = tmp.nr
        tmp.lz = tmp.nc
        tmp.nxstart = tmp.nsstart
        tmp.nystart = tmp.nrstart
        tmp.nzstart = tmp.ncstart
        def transpose(i, j, k):
            return (k, j, i)
        def reMap(i, j, k):
            return k + j * tmp.ns + i * tmp.ns * tmp.nr
    else:
        raise IOError("NULL: Error with Axis Remapping Data.")

    # Unroll electron density data into flat array
    tmp.chargeArray = numpy.zeros(tmp.lx * tmp.ly * tmp.lz)
    for n in range(len(rawArray)):
        k = n // ( taco.nc * taco.nr)
        j = (n % ( taco.nc * taco.nr )) // taco.nc
        i = ((n % ( taco.nc * taco.nr )) % taco.nc)
        tmp.chargeArray[n] = rawArray[reMap(i, j, k)]

    # fin!
    return tmp


