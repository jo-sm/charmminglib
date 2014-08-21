import numpy as np
cimport numpy as np
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


def build_cart(int lx, int ly, int lz,
                int x0, int y0, int z0,
                DTYPE_t xres, DTYPE_t yres, DTYPE_t zres):
    """
    x0 -> emap.nxstart
    y0 -> emap.nystart
    z0 -> emap.nzstart
    """
    cdef np.ndarray[DTYPE_t, ndim=2] tmp
    cdef int n, i, j, k
    #
    tmp = np.zeros((lx*ly*lz, 3), dtype=DTYPE)
    n = 0
    for k in range(z0, z0+lz):
        for j in range(y0, y0+ly):
            for i in range(x0, x0+lx):
                tmp[n, 0] = i * xres
                tmp[n, 1] = j * yres
                tmp[n, 2] = k * zres
                n += 1
    return tmp

def correlation(np.ndarray[DTYPE_t, ndim=1] dataArray1 not None,
                np.ndarray[DTYPE_t, ndim=1] dataArray2 not None):
    cdef int len1 = len(dataArray1)
    cdef int len2 = len(dataArray2)
    cdef int n = len1 if len1 < len2 else len2
    cdef int i
    cdef DTYPE_t mean1 = dataArray1.mean()
    cdef DTYPE_t mean2 = dataArray2.mean()
    cdef DTYPE_t std1 = dataArray1.std()
    cdef DTYPE_t std2 = dataArray2.std()
    cdef DTYPE_t tmpMean, tmp
    cdef np.ndarray[DTYPE_t, ndim=1] tmpArray
    #
    tmpArray = np.zeros(n, dtype=DTYPE)
    for i in range(n):
        tmpArray[i] = dataArray1[i] * dataArray2[i]
    #
    tmpMean = tmpArray.mean()
    tmp = (tmpMean - mean1 * mean2) / (std1 * std2)
    return tmp

def get_Tensor(np.ndarray[DTYPE_t, ndim=2] cartMatrix not None,
                np.ndarray[DTYPE_t, ndim=1] weightArray not None):
    cdef DTYPE_t xx = 0.
    cdef DTYPE_t yy = 0.
    cdef DTYPE_t zz = 0.
    cdef DTYPE_t xy = 0.
    cdef DTYPE_t xz = 0.
    cdef DTYPE_t yz = 0.
    cdef DTYPE_t x, y, z, q
    cdef int i
    cdef np.ndarray[DTYPE_t, ndim=2] tmp
    #
    for i in range(<int> len(cartMatrix)):
        x = cartMatrix[i, 0]
        y = cartMatrix[i, 1]
        z = cartMatrix[i, 2]
        q = weightArray[i]
        xx += q*(y*y+z*z)
        yy += q*(x*x+z*z)
        zz += q*(x*x+y*y)
        xy += q*x*y
        xz += q*x*z
        yz += q*y*z
    #
    tmp = np.array((
        ( xx, -xy, -xz),
        (-xy,  yy, -yz),
        (-xz, -yz,  zz)
        ), dtype=DTYPE)
    #
    return tmp

def interpolate(np.ndarray[DTYPE_t, ndim=2] cartMatrix not None,
                np.ndarray[DTYPE_t, ndim=1] weightArray not None,
                np.ndarray[DTYPE_t, ndim=3] targetArray not None,
                DTYPE_t xres, DTYPE_t yres, DTYPE_t zres):
    cdef DTYPE_t xCrd, yCrd, zCrd, x, y, z, q
    cdef int i, x0, x1, y0, y1, z0, z1
    #
    for i in range(<int> len(cartMatrix)):
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

def padZero(dataArray, int lx, int gx, int ly, int gy, int lz, int gz):
    if dataArray.dtype == np.float64:
        return padZeroFloat(dataArray, lx, gx, ly, gy, lz, gz)
    elif dataArray.dtype == np.int:
        return padZeroInt(dataArray, lx, gx, ly, gy, lz, gz)
    else:
        raise AssertionError

cdef padZeroFloat(np.ndarray[DTYPE_t, ndim=3] dataArray,
                int lx, int gx, int ly, int gy, int lz, int gz):
    cdef int x = dataArray.shape[0]
    cdef int y = dataArray.shape[1]
    cdef int z = dataArray.shape[2]
    cdef np.ndarray[DTYPE_t, ndim=3] tmp
    #
    tmp = np.zeros((x+lx+gx, y+ly+gy, z+lz+gz), dtype=DTYPE)
    if gx != 0 and gy != 0 and gz != 0:
        tmp[lx:-gx, ly:-gy, lz:-gz] = dataArray
    elif gx == 0 and gy == 0 and gz == 0:
        tmp[lx:, ly:, lz:] = dataArray
    elif gx != 0 and gy == 0 and gz == 0:
        tmp[lx:-gx, ly:, lz:] = dataArray
    elif gx == 0 and gy != 0 and gz == 0:
        tmp[lx:, ly:-gy, lz:] = dataArray
    elif gx == 0 and gy == 0 and gz != 0:
        tmp[lx:, ly:, lz:-gz] = dataArray
    elif gx != 0 and gy != 0 and gz == 0:
        tmp[lx:-gx, ly:-gy, lz:] = dataArray
    elif gx != 0 and gy == 0 and gz != 0:
        tmp[lx:-gx, ly:, lz:-gz] = dataArray
    elif gx == 0 and gy != 0 and gz != 0:
        tmp[lx:, ly:-gy, lz:-gz] = dataArray
    else:
        raise AssertionError
    return tmp

cdef padZeroInt(np.ndarray[np.int_t, ndim=3] dataArray,
                int lx, int gx, int ly, int gy, int lz, int gz):
    cdef int x = dataArray.shape[0]
    cdef int y = dataArray.shape[1]
    cdef int z = dataArray.shape[2]
    cdef np.ndarray[np.int_t, ndim=3] tmp
    #
    tmp = np.zeros((x+lx+gx, y+ly+gy, z+lz+gz), dtype=np.int)
    if gx != 0 and gy != 0 and gz != 0:
        tmp[lx:-gx, ly:-gy, lz:-gz] = dataArray
    elif gx == 0 and gy == 0 and gz == 0:
        tmp[lx:, ly:, lz:] = dataArray
    elif gx != 0 and gy == 0 and gz == 0:
        tmp[lx:-gx, ly:, lz:] = dataArray
    elif gx == 0 and gy != 0 and gz == 0:
        tmp[lx:, ly:-gy, lz:] = dataArray
    elif gx == 0 and gy == 0 and gz != 0:
        tmp[lx:, ly:, lz:-gz] = dataArray
    elif gx != 0 and gy != 0 and gz == 0:
        tmp[lx:-gx, ly:-gy, lz:] = dataArray
    elif gx != 0 and gy == 0 and gz != 0:
        tmp[lx:-gx, ly:, lz:-gz] = dataArray
    elif gx == 0 and gy != 0 and gz != 0:
        tmp[lx:, ly:-gy, lz:-gz] = dataArray
    else:
        raise AssertionError
    return tmp

def get_laplacian(np.ndarray[DTYPE_t, ndim=3] dataArray not None):
    cdef int x = dataArray.shape[0]
    cdef int y = dataArray.shape[1]
    cdef int z = dataArray.shape[2]
    cdef np.ndarray[DTYPE_t, ndim=3] array3d_1
    cdef np.ndarray[DTYPE_t, ndim=1] tmpArray
    cdef int i, j, k, n
    cdef DTYPE_t tmp
    #
    array3d_1 = padZeroFloat(dataArray, 1, 1, 1, 1, 1, 1)
    tmpArray = np.zeros(x*y*z, dtype=DTYPE)
    n = 0
    for k in range(1, 1+z):
        for j in range(1, 1+y):
            for i in range(1, 1+x):
                tmp = array3d_1[i+1, j, k]
                tmp += array3d_1[i-1, j, k]
                tmp += array3d_1[i, j+1, k]
                tmp += array3d_1[i, j-1, k]
                tmp += array3d_1[i, j, k+1]
                tmp += array3d_1[i, j, k-1]
                tmp -= 6.*array3d_1[i, j, k]
                tmpArray[n] = tmp
                n += 1
    return tmpArray

def get_core(np.ndarray[DTYPE_t, ndim=3] chargeArray not None,
            np.ndarray[DTYPE_t, ndim=3] lapArray not None,
            int rhocut):
    cdef int x = chargeArray.shape[0]
    cdef int y = chargeArray.shape[1]
    cdef int z = chargeArray.shape[2]
    cdef np.ndarray[np.int_t, ndim=3] tmpArray, nextArray
    cdef int n, i, j, k, tmp
    cdef int[6] tmpList
    # initialize
    tmpArray = np.ones((x-2, y-2, z-2), dtype=np.int)
    tmpArray = padZeroInt(tmpArray, 1, 1, 1, 1, 1, 1)
    #
    n = 0
    while True:
        nextArray = tmpArray.copy()
        for k in range(1, z-1):
            for j in range(1, y-1):
                for i in range(1, x-1):
                    if tmpArray[i, j, k] == 0:
                        continue
                    tmpList[0] = tmpArray[i+1, j, k]
                    tmpList[1] = tmpArray[i-1, j, k]
                    tmpList[2] = tmpArray[i, j+1, k]
                    tmpList[3] = tmpArray[i, j-1, k]
                    tmpList[4] = tmpArray[i, j, k+1]
                    tmpList[5] = tmpArray[i, j, k-1]
                    tmp = _min6(tmpList)
                    #
                    if chargeArray[i, j, k] <= rhocut and tmp == 0:
                        nextArray[i, j, k] = 0
                    elif lapArray[i, j, k] > 0 and tmp == 0:
                        nextArray[i, j, k] = 0
                    else:
                        nextArray[i, j, k] = tmp + 1
        #
        if (tmpArray == nextArray).all():
            return tmpArray.flatten('F')
        else:
            tmpArray = nextArray
            n += 1
            if n > 100:
                raise AssertionError("Solution has not converged after 100 iterations.")

cdef inline int _min6(int[6] x):
    cdef int i
    cdef int min = x[0]
    for i in range(1, 6):
        if x[i] < min:
            min = x[i]
    return min
