"""
ccp4 file specification as per --
http://www.ccp4.ac.uk/html/maplib.html#description
"""


import struct
import numpy
from emap import EMap


def get_chargeArrayFromCCP4File(filename):
    """
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
    tmp.filename = filename

    # Parse Header
    with open(tmp.filename, mode='rb') as fpointer:
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
        tmp.chargeArray = numpy.array(struct.unpack('f'*tmp.nc*tmp.nr*tmp.ns,
                                        fpointer.read(4*tmp.nc*tmp.nr*tmp.ns)))

    ## Map Column/Row/Section -> x/y/z
    #tmp.transposeTuple = (tmp.mapc - 1, tmp.mapr - 1, tmp.maps - 1)
    #tmp.chargeArray.resize((tmp.nc, tmp.nr, tmp.ns))
    #tmp.chargeArray = numpy.transpose(tmp.chargeArray, tmp.transposeTuple)
    ## Permutation Related Bookkeeping
    #tmp.lx, tmp.ly, tmp.lz = tmp.chargeArray.shape
    #start = [tmp.ncstart, tmp.nrstart, tmp.nsstart]
    #startDict = dict(zip(start, tmp.transposeTuple))
    #start.sort(key=lambda x: startDict[x])
    #tmp.nxstart, tmp.nystart, tmp.nzstart = start
    ## store flattened charge array
    #tmp.flatChargeArray = tmp.get_flatChargeArray()
    # Fin!
    return tmp


if __name__ == '__main__':


    import sys

    #print "\nThis currently only supports interactive debugging.\n"
    #print "`taco` = `EMap()`"
    #taco = get_chargeArrayFromCCP4File(sys.argv[1])

    tmp = EMap()
    tmp.filename = sys.argv[1]

    with open(tmp.filename, mode='rb') as fpointer:
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
        chargeArray = numpy.zeros((tmp.nc, tmp.nr, tmp.ns))
        n = 0
        for k in xrange(tmp.ns):
            for j in xrange(tmp.nr):
                for i in xrange(tmp.nc):
                    chargeArray[i][j][k] = rawArray[n]
                    n += 1
