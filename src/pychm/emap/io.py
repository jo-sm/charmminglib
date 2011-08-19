"""
"""


import struct
import numpy
from emap import EMap


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
                rawMatrix[i][j][k] = rawArray[n]
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
                tmp.chargeArray[n] = chargeMatrix[i][j][k]
                n += 1

    # Generate cartesian data into 'flat' array
    def pixel2cart_map(x, y, z):
        xcart = (tmp.nxstart + x) * (tmp.xlen / tmp.nx)
        ycart = (tmp.nystart + y) * (tmp.ylen / tmp.ny)
        zcart = (tmp.nzstart + z) * (tmp.zlen / tmp.nz)
        return numpy.array([xcart, ycart, zcart])

    tmp._cartArray = numpy.zeros((tmp.lx * tmp.ly * tmp.lz, 3))
    n = 0
    for k in xrange(tmp.lz):
        for j in xrange(tmp.ly):
            for i in xrange(tmp.lx):
                tmp._cartArray[n] = pixel2cart_map(i, j, k)
                n += 1

    # Fin!
    tmp._cartArray = tmp._cartArray
    tmp.chargeArray = tmp.chargeArray
    tmp.init_cart()
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


if __name__ == '__main__':


    import sys, itertools

    print "\nThis currently only supports interactive debugging.\n"
    print "`taco` = `EMap()`"

    taco = read_CCP4File(sys.argv[1])

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
    # taco2 = read_CCP4FileWXW(sys.argv[1])

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
