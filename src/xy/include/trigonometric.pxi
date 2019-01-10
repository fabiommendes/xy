from libc.math cimport cos as rcos, sin as rsin, tan as rtan
from math import pi

cdef double dg = 180 / pi
cdef double rad = pi / 180

cdef double inf = float('inf')
cdef double ninf = -float('inf')

cpdef inline double r360(double x):
    return x % 360

cpdef inline double r180(double x):
    return x % 180

cpdef inline double s180(double x):
    return x % 180

cpdef inline double r90(double x):
    return x % 90

cpdef inline double s90(double x):
    return x % 90


cpdef inline double cos(double x):
    """Stable cosine function with argument in degrees."""
    x = r360(x)
    if x == 0.0:
        return 1.0
    elif x == 90.0:
        return 0.0
    elif x == 180.0:
        return -1.0
    elif x == 270.0:
        return 0.0
    else:
        return rcos(x * rad)


cpdef inline double sin(double x):
    """Stable sine function with argument in degrees."""
    x = r360(x)
    if x == 0.0:
        return 0.0
    elif x == 90.0:
        return 1.0
    elif x == 180.0:
        return 0.0
    elif x == 270.0:
        return -1.0
    else:
        return rsin(x * rad)


cpdef inline double tan(double x):
    """Stable tangent function with argument in degrees."""
    x = s90(x)
    if x == 0.0:
        return 0.0
    elif x == 45.0:
        return 1.0
    elif x == 90.0:
        return inf
    elif x == -45.0:
        return -1.0
    elif x == -90.0:
        return ninf
    else:
        return rtan(x * rad)


cdef inline set_cs(double x, double* c, double* s):
    """Set sin and cos from given angle."""
    x = r360(x)
    cdef double angle
    if x == 0.0:
        c[0] = 1.0
        s[0] = 0.0
    elif x == 90.0:
        c[0] = 0.0
        s[0] = 1.0
    elif x == 180.0:
        c[0] = -1.0
        s[0] = 0.0
    elif x == 270.0:
        c[0] = 0.0
        s[0] = -1.0
    else:
        angle = x * rad
        c[0] = rcos(angle)
        s[0] = rsin(angle)


cdef inline set_rcs(double x, double* c, double* s):
    c[0] = rcos(x)
    s[0] = rsin(x)
