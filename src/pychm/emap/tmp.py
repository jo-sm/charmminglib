import numpy as np
import pychm.emap.ext as ext

def densityCorrelation(emap1, emap2):
    array1 = emap1.chargeArray
    array2 = emap2.chargeArray
    return ext.correlation(emap1.chargeArray, emap2.chargeArray)

def get_lap(emap):
    rho_1 = ext.padZero(emap.chargeArray3d, 1, 1, 1, 1, 1, 1)
    tmpArray = np.zeros(emap.lx * emap.ly * emap.lz, dtype=np.float64)
    iterator = ( (i, j, k) for k in xrange(1, 1+emap.lz)
                for j in xrange(1, 1+emap.ly) for i in xrange(1, 1+emap.lx) )
    for n, (i, j, k) in enumerate(iterator):
        tmp = np.float64(0)
        tmp += rho_1[i+1, j, k]
        tmp += rho_1[i-1, j, k]
        tmp += rho_1[i, j+1, k]
        tmp += rho_1[i, j-1, k]
        tmp += rho_1[i, j, k+1]
        tmp += rho_1[i, j, k-1]
        tmp -= 6*rho_1[i, j, k]
        tmpArray[n] = tmp
    return tmpArray

def get_core(emap):
    # initialize
    chargeMatrix = emap.chargeArray3d
    lapMatrix = emap.laplacianArray3d
    rhocut = emap.rhocut
    tmpArray = np.ones((emap.lx-2, emap.ly-2, emap.lz-2), dtype=np.int)
    tmpArray = ext.padZero(tmpArray, 1)
    #
    def minF(i, j, k):
        tmp = []
        tmp.append(tmpArray[i+1, j, k])
        tmp.append(tmpArray[i-1, j, k])
        tmp.append(tmpArray[i, j+1, k])
        tmp.append(tmpArray[i, j-1, k])
        tmp.append(tmpArray[i, j, k+1])
        tmp.append(tmpArray[i, j, k-1])
        return min(tmp)
    #
    iteration = 1
    while True:
        print "Iteration: %d" % iteration
        nextArray = tmpArray.copy()
        iterator = ( (i, j, k) for k in xrange(1, emap.lz-1) for j in xrange(1, emap.ly-1) for i in xrange(1, emap.lx-1) )
        for i, j, k in iterator:
            if chargeMatrix[i, j, k] <= rhocut and minF(i, j, k) == 0:
                nextArray[i, j, k] = 0
            elif lapMatrix[i, j, k] > 0 and minF(i, j, k) == 0:
                nextArray[i, j, k] = 0
            else:
                nextArray[i, j, k] = minF(i, j, k) + 1
        #
        if (tmpArray == nextArray).all():
            return tmpArray.flatten('F')
        else:
            iteration += 1
            tmpArray = nextArray

def almost_equal(x, y, eps=6):
    eps = 10**-eps
    for i in xrange(len(x)):
        if not abs(x[i] - y[i]) < eps:
            print "%06d %12.8f %12.8f" % (i, x[i], y[i])

def interpolate(cartMatrix, weightArray, targetArray,
                xres, yres, zres):
    for i in range(len(cartMatrix)):
        q = weightArray[i]
        xCrd = cartMatrix[i, 0]
        yCrd = cartMatrix[i, 1]
        zCrd = cartMatrix[i, 2]
        x0 = int(xCrd // xres)
        y0 = int(yCrd // yres)
        z0 = int(zCrd // zres)
        x1 = x0 + 1
        y1 = y0 + 1
        z1 = z0 + 1
        x = xCrd % xres / xres
        y = yCrd % yres / yres
        z = zCrd % zres / zres
        #
        targetArray[x0, y0, z0] += q * (1 - x) * (1 - y) * (1 - z)
        targetArray[x1, y0, z0] += q * x * (1 - y) * (1 - z)
        targetArray[x0, y1, z0] += q * (1 - x) * y * (1 - z)
        targetArray[x0, y0, z1] += q * (1 - x) * (1 - y) * z
        targetArray[x1, y1, z0] += q * x * y * (1 - z)
        targetArray[x1, y0, z1] += q * x * (1 - y) * z
        targetArray[x0, y1, z1] += q * (1 - x) * y * z
        targetArray[x1, y1, z1] += q * x * y * z
    return targetArray

def trimZeros(emap1):
    # Initialize
    imax = 0
    jmax = 0
    kmax = 0
    imin = emap1.lx
    jmin = emap1.ly
    kmin = emap1.lz
    #
    for k in range(emap1.lz):
        for j in range(emap1.ly):
            for i in range(emap1.lx):
                if emap1.chargeArray3d[i, j, k] == 0.:
                    continue
                if i > imax:
                    imax = i
                elif i < imin:
                    imin = i
                if j > jmax:
                    jmax = j
                elif j < jmin:
                    jmin = j
                if k > kmax:
                    kmax = k
                elif k < kmin:
                    kmin = k
    #
    newmap = emap1._new()
    newmap.nxstart = emap1.nxstart
    newmap.nystart = emap1.nystart
    newmap.nzstart = emap1.nzstart
    newmap.ncstart = emap1.nxstart
    newmap.nrstart = emap1.nystart
    newmap.nsstart = emap1.nzstart
    #
    chargeArray3d = emap1.chargeArray3d[imin:imax+1, jmin:jmax+1, kmin:kmax+1]
    newmap.lx, newmap.ly, newmap.lz = chargeArray3d.shape
    newmap.nc, newmap.nr, newmap.ns = chargeArray3d.shape
    cartArray3d = emap1.cartArray3d[imin:imax+1, jmin:jmax+1, kmin:kmax+1]
    newmap.chargeArray = chargeArray3d.flatten('F')
    newmap.cartArray = np.reshape(cartArray3d, (newmap.lx * newmap.ly * newmap.lz, 3), order='F')
    #
    newmap.nc = newmap.lx
    newmap.nr = newmap.ly
    newmap.ns = newmap.lz
    #
    return newmap

def combine_maps(emap1, emap2):
    """
    Loop over all the grid points of emap1, and interpolate them into grid points of emap2.
    """
    # (1) Find Min/Max cartesian values of emap1 and emap2
    xm1, ym1, zm1 = emap1.cartArray.min(axis=0)
    xm2, ym2, zm2 = emap2.cartArray.min(axis=0)
    #
    xmin = xm1 if xm1 < xm2 else xm2
    ymin = ym1 if ym1 < ym2 else ym2
    zmin = zm1 if zm1 < zm2 else zm2
    #
    xm1, ym1, zm1 = emap1.cartArray.max(axis=0)
    xm2, ym2, zm2 = emap2.cartArray.max(axis=0)
    #
    xmax = xm1 if xm1 > xm2 else xm2
    ymax = ym1 if ym1 > ym2 else ym2
    zmax = zm1 if zm1 > zm2 else zm2

    # (2) create valid metadata
    newmap = emap2._new()

    newmap.nxstart = np.int(np.floor(xmin/emap2.xres)) - 1
    newmap.nystart = np.int(np.floor(ymin/emap2.yres)) - 1
    newmap.nzstart = np.int(np.floor(zmin/emap2.zres)) - 1
    newmap.lx = np.int(np.ceil(xmax/emap2.xres)) - newmap.nxstart + 1
    newmap.ly = np.int(np.ceil(ymax/emap2.yres)) - newmap.nystart + 1
    newmap.lz = np.int(np.ceil(zmax/emap2.zres)) - newmap.nzstart + 1

    newmap.nc = newmap.lx
    newmap.nr = newmap.ly
    newmap.ns = newmap.lz
    newmap.ncstart = newmap.nxstart
    newmap.nrstart = newmap.nystart
    newmap.nsstart = newmap.nzstart

    # (3) Create a new blank 3d array using the min/max cartesian values
    # (4) copy data from emap2
    lessX = emap2.nxstart - newmap.nxstart
    moreX = newmap.lx + newmap.nxstart - (emap2.lx + emap2.nxstart)
    lessY = emap2.nystart - newmap.nystart
    moreY = newmap.ly + newmap.nystart - (emap2.ly + emap2.nystart)
    lessZ = emap2.nzstart - newmap.nzstart
    moreZ = newmap.lz + newmap.nzstart - (emap2.lz + emap2.nzstart)

    chargeArray = ext.padZero(emap2.chargeArray3d,
                            lessX, moreX, lessY, moreY, lessZ, moreZ)
    newmap.cartArray = ext.build_cart(newmap.lx, newmap.ly, newmap.lz,
                                       newmap.nxstart, newmap.nystart, newmap.nzstart,
                                       newmap.xres, newmap.yres, newmap.zres)

    # (5) interpolate data from emap1
    newmap.chargeArray = ext.interpolate(emap1.cartArray, emap1.chargeArray, chargeArray,
                            newmap.xres, newmap.yres, newmap.zres).flatten('F')

    return newmap

def combine(emapX, emapY):
    # (1) Find Min/Max cartesian values of emap1 and emap2
    xmin1, ymin1, zmin1 = emapX.cartArray.min(axis=0)
    xmin2, ymin2, zmin2 = emapY.cartArray.min(axis=0)
    xmax1, ymax1, zmax1 = emapX.cartArray.max(axis=0)
    xmax2, ymax2, zmax2 = emapY.cartArray.max(axis=0)
    xmin = xmin1 if xmin1 < xmin2 else xmin2
    ymin = ymin1 if ymin1 < ymin2 else ymin2
    zmin = zmin1 if zmin1 < zmin2 else zmin2
    xmax = xmax1 if xmax1 > xmax2 else xmax2
    ymax = ymax1 if ymax1 > ymax2 else ymax2
    zmax = zmax1 if zmax1 > zmax2 else zmax2
    #
    lessX = xmin2 - xmin
    moreX = xmax2 - xmax
    lessY = ymin2 - ymin
    moreY = ymax2 - ymax
    lessZ = zmin2 - zmin
    moreZ = zmax2 - zmax
    #
    tmp = emapY.pad_mapByCart(lessX, moreX, lessY, moreY, lessZ, moreZ)
    tmp.chargeArray = ext.interpolate(emapX.cartArray, emapX.chargeArray,
                                    tmp.chargeArray3d, tmp.xres, tmp.yres,
                                    tmp.zres).flatten('F')
    return tmp

if __name__ == '__main__':
    import io
    import copy
    from emap import EMap
    import matplotlib.pyplot as pp


    #emap2 = io.CCP4File.read('/v/bigbox5/home/fpickard/python/package_prep/pychm/src/pychm/emap/1a7n.ccp4')
    #emap2 = io.CCP4File.read('/home/fpickard/chien/map2.ccp4')
    #emap1 = io.CCP4File.read('/home/fpickard/chien/map1.ccp4')
    #emap2 = copy.deepcopy(emap1)
    #emap2.translations.append(np.array((-140., -140., -140.)))

    def plotXY(emap1):
        newCart1 = emap1.cartArray[:, :2]
        newCharge1 = np.zeros((emap1.shape[0] * emap1.shape[1]), dtype=np.float64)
        newChargeLen1 = newCharge1.size
        for i, charge in enumerate(emap1.chargeArray):
            newCharge1[i % newChargeLen1] += charge
        newCharge1 = np.reshape(newCharge1, (emap1.shape[0], emap1.shape[1]), 'F')
        x = np.arange(emap1.cartArray.min(axis=0)[0], emap1.cartArray.max(axis=0)[0] + 0.1 * emap1.xres, emap1.xres)
        y = np.arange(emap1.cartArray.min(axis=0)[1], emap1.cartArray.max(axis=0)[1] + 0.1 * emap1.yres, emap1.yres)
        pp.contour(x, y, newCharge1.T)
        pp.show()

    #newCart2 = emap2.cartArray[:, :2]
    #newCharge2 = np.zeros((emap2.shape[0] * emap2.shape[1]), dtype=np.float64)
    #newChargeLen2 = newCharge2.size
    #for i, charge in enumerate(emap2.chargeArray):
    #    newCharge2[i % newChargeLen2] += charge
    #newCharge2 = np.reshape(newCharge2, (emap2.shape[0], emap2.shape[1]), 'F')

    emapX = io.CCP4File.read('/home/fpickard/chien/map1.ccp4')
    emapY = copy.deepcopy(emapX)
    emapY.translations.append(np.array((-40., -40., -40.)))
    # (1) Find Min/Max cartesian values of emap1 and emap2
    xmin1, ymin1, zmin1 = emapX.cartArray.min(axis=0)
    xmin2, ymin2, zmin2 = emapY.cartArray.min(axis=0)
    xmax1, ymax1, zmax1 = emapX.cartArray.max(axis=0)
    xmax2, ymax2, zmax2 = emapY.cartArray.max(axis=0)
    xmin = xmin1 if xmin1 < xmin2 else xmin2
    ymin = ymin1 if ymin1 < ymin2 else ymin2
    zmin = zmin1 if zmin1 < zmin2 else zmin2
    xmax = xmax1 if xmax1 > xmax2 else xmax2
    ymax = ymax1 if ymax1 > ymax2 else ymax2
    zmax = zmax1 if zmax1 > zmax2 else zmax2
    #
    lessX = xmin2 - xmin
    moreX = xmax - xmax2
    lessY = ymin2 - ymin
    moreY = ymax - ymax2
    lessZ = zmin2 - zmin
    moreZ = zmax - zmax2
    #
    tmp = emapY.pad_mapByCart(lessX, moreX, lessY, moreY, lessZ, moreZ)
    tmp.chargeArray = ext.interpolate(emapX.cartArray, emapX.chargeArray,
                                    tmp.chargeArray3d, tmp.xres, tmp.yres,
                                    tmp.zres).flatten('F')
