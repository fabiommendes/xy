# distutils: language = c++
# cython: cdivision   = True

include "../include/trigonometric.pxi"
include "../include/text.pxi"
cimport cython
from libc.math cimport sqrt, atan2, fabs
from cpython.object cimport PyObject_TypeCheck


# ------------------------------------------------------------------------------
# VECTOR 2D CLASS
# ------------------------------------------------------------------------------

# noinspection PyMethodParameters
@cython.freelist(8)
cdef class Vec2:
    """
    Simple 2D Vector
    """
    ndim = 2
    size = 2
    dtype = float
    shape = (2,)

    #
    # Properties and accessors
    #
    @property
    def length(self):
        return length(self.data)

    @property
    def length_sqr(self):
        return length2(self.data)

    @property
    def angle(self):
        if self.y == 0:
            return 0.0 if self.data.x >= 0 else 180.0
        if self.x == 0:
            return 90 if self.data.y >= 0 else 270.0
        return atan2(self.data.y, self.data.x) * dg

    @property
    def angle_rad(self):
        return atan2(self.data.y, self.data.x)

    @property
    def direction(self):
        return direction2(self.data)

    @property
    def x(self):
        return self.data.x

    @property
    def y(self):
        return self.data.y

    #
    # Constructors
    #
    @classmethod
    def from_flat(cls, data):
        """Create vector from 2-sequence."""
        x, y = data
        return newvec2(x, y)

    @classmethod
    def from_polar(cls, float radius, float angle):
        """Create vector from polar coordinates."""
        return newvec2(radius * cos(angle), radius * sin(angle))

    @classmethod
    def from_rpolar(cls, float radius, float angle):
        """
        Create vector from polar coordinates.

        Uses radians instead of degrees to define angle.
        """
        return newvec2(radius * rcos(angle), radius * rsin(angle))

    def __cinit__(self, double x, double y):
        self.data = dvec2(x, y)

    #
    # Generic API
    #
    def copy(self, x=None, y=None):
        """
        Creates a new copy of vector overriding either the x, y or both
        components.
        """
        cdef double x_, y_
        set_xy(self, &x_, &y_)
        if x is not None:
            x_ = x
        if y is not None:
            y_ = y
        return newvec2(x_, y_)

    #
    # Magic methods
    #
    def __getstate__(self):
        return self.x, self.y

    def __setstate__(self, st):
        cdef double x, y
        x, y = st
        self.data = dvec2(x, y)

    def __repr__(self):
        return f'{type(self).__name__}({fmt(self.x)}, {fmt(self.y)})'

    def __len__(self):
        return 2

    def __iter__(self):
        yield self.data.x
        yield self.data.y

    def __getitem__(self, idx):
        cdef int i
        try:
            i = idx
        except TypeError:
            if isinstance(idx, slice):
                return [self.data.x, self.data.y][idx]
        else:
            if i == 0 or i == -2:
                return x(self)
            elif i == 1 or i == -1:
                return y(self)
        raise IndexError(idx)

    def __neg__(self):
        return vec2(-self.data)

    def __pos__(self):
        return self

    def __abs__(self):
        return length(self.data)

    def __bool__(self):
        return self.data.x != 0.0 or self.data.y != 0.0

    def __eq__(u, v):
        if isvec(u) and isvec(v):
            return (<Vec2> u).data == (<Vec2> v).data
        try:
            return tovec2(u).data == tovec2(v).data
        except TypeError:
            return NotImplemented

    def __add__(u, v):
        if isvec(u) and isvec(v):
            return vec2((<Vec2> u).data + (<Vec2> v).data)
        try:
            return vec2(tovec2(u).data + tovec2(v).data)
        except TypeError:
            return NotImplemented

    def __sub__(u, v):
        if isvec(u) and isvec(v):
            return vec2((<Vec2> u).data - (<Vec2> v).data)
        try:
            return vec2(tovec2(u).data - tovec2(v).data)
        except TypeError:
            return NotImplemented

    def __mul__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec2((<Vec2> u).data * m)
            else:
                m = u
                return vec2((<Vec2> v).data * m)
        except TypeError:
            return NotImplemented

    def __truediv__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec2((<Vec2> u).data / m)
            else:
                return NotImplemented
        except TypeError:
            return NotImplemented

    def __floordiv__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec2(floor((<Vec2> u).data / m))
            else:
                return NotImplemented
        except TypeError:
            return NotImplemented

    #
    # Exclusive 2D functions
    #
    def rotate(self, double angle):
        """Rotate vector around origin."""

        cdef double cos, sin
        set_cs(angle, &cos, &sin)
        return newvec2(self.x * cos - self.y * sin, self.x * sin + self.y * cos)

    def rotate_around(self, double angle, axis):
        """Rotate vector around given axis."""

        cdef double cos, sin, x0, y0, dx, dy
        set_xy(axis, &x0, &y0)
        set_cs(angle, &cos, &sin)
        dx = self.x - x0
        dy = self.y - y0
        return newvec2(x0 + dx * cos - dy * sin, y0 + dx * sin + dy * cos)

    def rotate_rad(self, double angle):
        """Rotate vector around origin (angle in radians)."""

        cdef double cos = rcos(angle), sin = rsin(angle)
        return newvec2(self.x * cos - self.y * sin, self.x * sin + self.y * cos)

    def rotate_around_rad(self, double angle, axis):
        """Rotate vector around origin (angle in radians)."""

        cdef double cos = rcos(angle), sin = rsin(angle), x0, y0, dx, dy
        set_xy(axis, &x0, &y0)
        dx = self.x - x0
        dy = self.y - y0
        return newvec2(x0 + dx * cos - dy * sin, y0 + dx * sin + dy * cos)

    def cross(self, other):
        """
        The z component of the cross product between two bi-dimensional
        smallvectors.
        """
        return cross(self.data, tovec2(other).data)

    def polar(self):
        """
        Return a tuple with the (radius, angle) polar coordinates.

        Angle is returned in degrees.
        """
        return length(self.data), atan2(self.data.y, self.data.x) * dg

    def polar_rad(self):
        """
        Return a tuple with the (radius, theta) polar coordinates.

        Angle is returned in radians.
        """
        return length(self.data), atan2(self.data.y, self.data.x)

    def perpendicular(self, ccw=True):
        """
        Return the counterclockwise perpendicular vector.

        If ccw is False, do the rotation in the clockwise direction.
        """
        return newvec2(-self.data.y, self.data.x) if ccw else newvec2(self.data.y, -self.data.x)

    #
    # Generic vector functions
    #
    def dot(self, other):
        """Dot product with another vector."""

        return dot(self.data, tovec2(other).data)

    def distance_to(self, other):
        """Compute distance from vector."""

        return distance(self.data, tovec2(other).data)

    def angle_to(self, other):
        """Return the angle to other vector."""

        return angle(self.data, tovec2(other).data) * dg

    def angle_to_rad(self, other):
        """Return the angle to other vector (measured in radians)."""

        return angle(self.data, tovec2(other).data)

    def reflect(self, direction):
        """Reflection of vector around given normal direction."""

        cdef dvec2 *n
        if isdirection(direction):
            n = &(<Vec2> direction).data
            return vec2(self.data - 2 * (self.data - dot(n[0], self.data) * n[0]))
        elif isvec(direction):
            n = &(<Vec2> direction).data
            return vec2(self.data - 2 * (self.data - vprojection(self.data, n[0])))
        else:
            return self.reflect(tovec2(direction))

    def projection(self, direction):
        """Projection vector for the given direction."""
        cdef dvec2 *n
        if isdirection(direction):
            n = &(<Vec2> direction).data
            return vec2(dot(self.data, n[0]) * n[0])
        else:
            return vec2(vprojection(self.data, tovec2(direction).data))

    def clamp(self, double min, double max):
        """Set length between minimum and maximum values."""

        cdef double len = length(self.data)
        if min <= len <= max:
            return self
        elif len < min:
            return vec2(self.data * (min / len))
        elif len > max:
            return vec2(self.data * (max / len))
        else:
            raise ValueError('cannot re-scale a zero-length vector.')

    def with_length(self, double size):
        """Return a copy with the given length."""

        cdef double len = length(self.data)
        if len != 0:
            return vec2(self.data * (size / len))
        else:
            raise ValueError('cannot re-scale a zero-length vector.')

    def lerp(self, other, double t=0.5):
        """
        Linear interpolation o points.

        Maps the range [0, 1] to a path that goes from self to other.

            t = 0  ==> self
            t = 1  ==> other
            others ==> a proportional combination of self and other
        """

        cdef dvec2 vec
        set_vec(other, &vec)
        return vec2(mix(self.data, vec, t))

    def slerp(self, other, double t=0.5):
        """
        Spherical interpolation o points.

        Like lerp, but uses polar coordinates.

            t = 0  ==> self
            t = 1  ==> other
            others ==> a proportional combination of self and other

        Notes:
            If self and other have the same length, slerp preserves the length
            and interpolates the angle.

        See also:
            https://en.wikipedia.org/wiki/Slerp
        """
        cdef dvec2 u = self.data, v
        set_vec(other, &v)

        # Return lerp on special cases
        ru = length(u)
        rv = length(v)
        if ru == 0 or rv == 0:
            return vec2(mix(u, v, t))

        # Compute slerp on the happy path (here, we are always using radians)
        sin_omega = cross(u, v) / (ru * rv)
        omega = atan2(sin_omega, dot(u, v))
        a = rsin((1 - t) * omega) / sin_omega
        b = rsin(t * omega) / sin_omega
        return vec2(u * a + v * b)

    def midpoint(self, other):
        """Midpoint between two vectors."""

        cdef dvec2 vec
        set_vec(other, &vec)
        return vec2(mix(self.data, vec, 0.5))

    #
    # Queries
    #
    def is_null(self):
        """Return True if vector is null."""

        return self.data.x == 0.0 and self.data.y == 0.0

    def is_almost_null(self, tol=1e-6):
        """Return True if vector is null within given tolerance."""

        return fabs(length(self.data)) < tol

    def is_normalized(self, tol=1e-6):
        """Return True if is normalized within the given tolerance."""

        return fabs(length(self.data) - 1) < tol

    def is_almost_equal(self, other, double tol=1e-6):
        """Return True if two vectors are equal under the given tolerance."""

        return length(self.data - tovec2(other).data) <= tol

    #
    # Norm
    #
    def norm(self):
        return length(self.data)

    def norm_l1(self):
        return fabs(self.data.x) + fabs(self.data.y)

    def norm_linf(self):
        return max(fabs(self.data.x), fabs(self.data.y))

    def norm_sqr(self):
        return length2(self.data)

    def normalize(self):
        return vec2(normalize(self.data))


# ------------------------------------------------------------------------------
# DIRECTION 2D CLASS
# ------------------------------------------------------------------------------

@cython.final(True)
cdef class Direction2(Vec2):
    """
    Normalized 2D vectors.
    """

    @property
    def length(self):
        return 1

    @property
    def length_sqr(self):
        return 1

    def __cinit__(self, double x, double y):
        cdef double norm = sqrt(x * x + y * y)
        if norm == 0:
            raise ValueError('cannot initialize direction from zero-length coordinates')
        self.data = dvec2(x / norm, y / norm)

    def slerp(self, other, double t=0.5):
        if not isdirection(other):
            return super().slerp(other, t)

        cdef dvec2 u = self.data, v
        set_vec(other, &v)

        sin_omega = cross(u, v)
        omega = atan2(sin_omega, dot(u, v))
        a = rsin((1 - t) * omega) / sin_omega
        b = rsin(t * omega) / sin_omega
        return unsafe_direction2(u * a + v * b)

    def is_unity(self, tol=1e-6):
        return True

    def normalize(self, tol=1e-6):
        return self

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors
cdef inline Vec2 newvec2(double x, double y):
    return vec2(dvec2(x, y))

cdef inline Vec2 newdirection2(double x, double y):
    return direction2(dvec2(x, y))

cdef inline Vec2 vec2(dvec2 v):
    cdef PyObject* new = _PyObject_New(Vector2Type)
    (<Vec2> new).data = v
    return <Vec2> new

cdef inline Direction2 direction2(dvec2 v):
    cdef PyObject* new = _PyObject_New(Direction2Type)
    (<Direction2 > new).data = normalize(v)
    return <Direction2 > new

cdef inline Direction2 unsafe_direction2(dvec2 v):
    cdef PyObject* new = _PyObject_New(Direction2Type)
    (<Direction2 > new).data = v
    return <Direction2 > new


cdef inline Vec2 tovec2(u):
    cdef double x, y

    if isvec(u):
        return <Vec2> u
    elif PyObject_TypeCheck(u, TupleType) and len(<tuple> u) == 2:
        with cython.boundscheck(False):
            x = (<tuple> u)[0]
            y = (<tuple> u)[1]
        return newvec2(x, y)
    else:
        raise TypeError

# Queries
cdef inline bint isvec(u):
    return PyObject_TypeCheck(u, Vector2Type)

cdef inline bint isdirection(u):
    return PyObject_TypeCheck(u, Direction2Type)

cdef inline bint istuple(u):
    return PyObject_TypeCheck(u, TupleType)

# Data accessors
cdef inline double x(Vec2 u):
    return u.data.x

cdef inline double y(Vec2 u):
    return u.data.y

cdef inline void set_xy(u, double *x, double *y):
    if isvec(u):
        x[0] = (<Vec2> u).x
        y[0] = (<Vec2> u).y
    elif istuple(u) and len(<tuple> u) != 2:
        with cython.boundscheck(False):
            x[0] = (<tuple> u)[0]
            y[0] = (<tuple> u)[1]
    else:
        raise ValueError('Requires a Vec2 or a 2-tuple')

cdef inline void set_vec(u, dvec2 *vec):
    if isvec(u):
        vec[0] = (<Vec2> u).data
    elif istuple(u) and len(<tuple> u) != 2:
        with cython.boundscheck(False):
            vec[0].x = (<tuple> u)[0]
            vec[0].x = (<tuple> u)[1]
    else:
        raise ValueError('Requires a Vec2 or a 2-tuple')

# Vector functions
cdef inline double length2(dvec2 v):
    return v.x * v.x + v.y * v.y

cdef inline double norm_l1(dvec2 v):
    return fabs(v.x) + fabs(v.y)

cdef inline double cross(dvec2 u, dvec2 v):
    return u.x * v.y - u.y * v.x

cdef inline double angle(dvec2 u, dvec2 v):
    cdef double cos_t = dot(u, v), \
                sin_t = cross(u, v)
    return atan2(sin_t, cos_t)

cdef inline double projection(dvec2 v, dvec2 n):
    return dot(v, n) / length(n)

cdef inline dvec2 vprojection(dvec2 v, dvec2 n):
    return (dot(v, n) / length2(n)) * n

cdef inline dvec2 polar(dvec2 v):
    return dvec2(length(v), atan2(v.y, v.x))

# Constants
cdef Vec2 origin = Vec2(0, 0)
cdef PyTypeObject *Vector2Type = <PyTypeObject*> Vec2
cdef PyTypeObject *Direction2Type = <PyTypeObject*> Direction2
cdef PyTypeObject *TupleType = <PyTypeObject*> tuple