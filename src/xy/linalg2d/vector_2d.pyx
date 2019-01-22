# distutils: language = c++
# cython: cdivision   = True

include "../include/trigonometric.pxi"
include "../include/text.pxi"

# noinspection PyUnresolvedReferences
cimport cython
from libc.math cimport sqrt, atan2, fabs


# ------------------------------------------------------------------------------
# VECTOR 2D CLASS
# ------------------------------------------------------------------------------

# noinspection PyMethodParameters
@cython.freelist(8)
cdef class Vec2:
    """
    Simple 2D Vector
    """
    ndim = 1
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
        return direction(self.data)

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
        return flatvec(x, y)

    @classmethod
    def from_polar(cls, float radius, float angle):
        """Create vector from polar coordinates."""
        return flatvec(radius * cos(angle), radius * sin(angle))

    @classmethod
    def from_rpolar(cls, float radius, float angle):
        """
        Create vector from polar coordinates.

        Uses radians instead of degrees to define angle.
        """
        return flatvec(radius * rcos(angle), radius * rsin(angle))

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
        return flatvec(self.data.x if x is None else x,
                       self.data.y if y is None else y)

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
        return vec(-self.data)

    def __pos__(self):
        return self

    def __abs__(self):
        return length(self.data)

    def __eq__(u, v):
        if isvec(u) and isvec(v):
            return (<Vec2> u).data == (<Vec2> v).data
        try:
            return tovec(u).data == tovec(v).data
        except TypeError:
            return NotImplemented

    def __add__(u, v):
        if isvec(u) and isvec(v):
            return vec((<Vec2> u).data + (<Vec2> v).data)
        try:
            return vec(tovec(u).data + tovec(v).data)
        except TypeError:
            return NotImplemented

    def __sub__(u, v):
        if isvec(u) and isvec(v):
            return vec((<Vec2> u).data - (<Vec2> v).data)
        try:
            return vec(tovec(u).data - tovec(v).data)
        except TypeError:
            return NotImplemented

    def __mul__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec((<Vec2> u).data * m)
            else:
                m = u
                return vec((<Vec2> v).data * m)
        except TypeError:
            return NotImplemented

    def __truediv__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec((<Vec2> u).data / m)
            else:
                return NotImplemented
        except TypeError:
            return NotImplemented

    def __floordiv__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec(floor((<Vec2> u).data / m))
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
        return flatvec(self.x * cos - self.y * sin, self.x * sin + self.y * cos)

    def rotate_around(self, double angle, axis):
        """Rotate vector around given axis."""

        cdef double cos, sin, x0 = 0, y0 = 0, dx, dy
        set_xy(axis, &x0, &y0)
        set_cs(angle, &cos, &sin)
        dx = self.x - x0
        dy = self.y - y0
        return flatvec(x0 + dx * cos - dy * sin, y0 + dx * sin + dy * cos)

    def rotate_rad(self, double angle):
        """Rotate vector around origin (angle in radians)."""

        cdef double cos = rcos(angle), sin = rsin(angle)
        return flatvec(self.x * cos - self.y * sin, self.x * sin + self.y * cos)

    def rotate_around_rad(self, double angle, axis):
        """Rotate vector around origin (angle in radians)."""

        cdef double cos = rcos(angle), sin = rsin(angle), x0 = 0, y0 = 0, dx, dy
        set_xy(axis, &x0, &y0)
        dx = self.x - x0
        dy = self.y - y0
        return flatvec(x0 + dx * cos - dy * sin, y0 + dx * sin + dy * cos)

    def cross(self, other):
        """
        The z component of the cross product between two bi-dimensional
        smallvectors.
        """
        return cross(self.data, tovec(other).data)

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
        return flatvec(-self.data.y, self.data.x) if ccw else flatvec(self.data.y, -self.data.x)

    #
    # Generic vector functions
    #
    def dot(self, other):
        """Dot product with another vector."""

        return dot(self.data, tovec(other).data)

    def distance_to(self, other):
        """Compute distance from vector."""

        return distance(self.data, tovec(other).data)

    def angle_to(self, other):
        """Return the angle to other vector."""

        return angle(self.data, tovec(other).data) * dg

    def angle_to_rad(self, other):
        """Return the angle to other vector (measured in radians)."""

        return angle(self.data, tovec(other).data)

    def reflect(self, direction):
        """Reflection of vector around given normal direction."""

        cdef dvec2 *n

        if isdirection(direction):
            n = &(<Vec2> direction).data
            return vec(self.data - 2 * (self.data - dot(n[0], self.data) * n[0]))
        elif isvec(direction):
            n = &(<Vec2> direction).data
            return vec(self.data - 2 * (self.data - vprojection(self.data, n[0])))
        else:
            return self.reflect(tovec(direction))

    def projection(self, direction):
        """Projection vector for the given direction."""

        cdef dvec2 *n

        if isdirection(direction):
            n = &(<Vec2> direction).data
            return vec(dot(self.data, n[0]) * n[0])
        else:
            return vec(vprojection(self.data, tovec(direction).data))

    # noinspection PyShadowingBuiltins
    def clamp(self, double min, double max):
        """Set length between minimum and maximum values."""

        cdef double norm = length(self.data)

        if min <= norm <= max:
            return self
        elif norm < min:
            return vec(self.data * (min / norm))
        elif norm > max:
            return vec(self.data * (max / norm))
        else:
            raise ValueError('cannot re-scale a zero-length vector.')

    def with_length(self, double size):
        """Return a copy with the given length."""

        cdef double norm = length(self.data)
        if norm != 0:
            return vec(self.data * (size / norm))
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

        cdef dvec2 u = self.data, v
        set_vec(other, &v)
        return vec(mix(u, v, t))

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
            return vec(mix(u, v, t))

        # Compute slerp on the happy path (here, we are always using radians)
        sin_omega = cross(u, v) / (ru * rv)
        omega = atan2(sin_omega, dot(u, v))
        a = rsin((1 - t) * omega) / sin_omega
        b = rsin(t * omega) / sin_omega
        return vec(u * a + v * b)

    def midpoint(self, other):
        """Midpoint between two vectors."""

        cdef dvec2 u = self.data, v
        set_vec(other, &v)
        return vec(mix(u, v, 0.5))

    #
    # Queries
    #
    def is_null(self):
        """Return True if vector is null."""

        return isnull(self.data)

    def is_almost_null(self, tol=1e-6):
        """Return True if vector is null within given tolerance."""

        return fabs(length(self.data)) < tol

    def is_normalized(self, tol=1e-6):
        """Return True if is normalized within the given tolerance."""

        return fabs(length(self.data) - 1) < tol

    def is_almost_equal(self, other, double tol=1e-6):
        """Return True if two vectors are equal under the given tolerance."""

        return length(self.data - tovec(other).data) <= tol

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
        return vec(normalize(self.data))


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
        return unsafe_direction(u * a + v * b)

    def is_unity(self, tol=1e-6):
        return True

    def normalize(self, tol=1e-6):
        return self


# Constants --------------------------------------------------------------------
cdef Vec2 origin = Vec2(0, 0)
Vector2Type = <PyTypeObject*> Vec2
Direction2Type = <PyTypeObject*> Direction2
TupleType = <PyTypeObject*> tuple
