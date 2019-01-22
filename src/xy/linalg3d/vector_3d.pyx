# distutils: language = c++
# cython: cdivision   = True

include "../include/trigonometric.pxi"
include "../include/text.pxi"
cimport cython
from libc.math cimport sqrt, atan2, fabs
cimport xy.linalg2d.vector_2d as V2


# ------------------------------------------------------------------------------
# VECTOR 2D CLASS
# ------------------------------------------------------------------------------

# noinspection PyMethodParameters
@cython.freelist(8)
cdef class Vec3:
    """
    Simple 2D Vector
    """
    ndim = 1
    size = 3
    dtype = float
    shape = (3,)

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
    def direction(self):
        return direction(self.data)

    @property
    def x(self):
        return self.data.x

    @property
    def y(self):
        return self.data.y

    @property
    def z(self):
        return self.data.z

    @property
    def xy(self):
        return V2.flatvec(self.data.x, self.data.y)

    @property
    def yx(self):
        return V2.flatvec(self.data.y, self.data.x)

    @property
    def xz(self):
        return V2.flatvec(self.data.x, self.data.z)

    @property
    def zx(self):
        return V2.flatvec(self.data.z, self.data.x)

    @property
    def yz(self):
        return V2.flatvec(self.data.y, self.data.z)

    @property
    def zy(self):
        return V2.flatvec(self.data.z, self.data.y)


    #
    # Constructors
    #
    @classmethod
    def from_flat(cls, data):
        """Create vector from 3-sequence."""
        x, y, z = data
        return flatvec(x, y, z)

    def __cinit__(self, double x, double y, double z):
        self.data = dvec3(x, y, z)

    #
    # Generic API
    #
    def copy(self, x=None, y=None, z=None):
        """
        Creates a new copy of vector overriding any of the x, y, z components.
        """
        return flatvec(self.data.x if x is None else x,
                       self.data.y if y is None else y,
                       self.data.z if z is None else z)


    #
    # Magic methods
    #
    def __getstate__(self):
        return self.x, self.y, self.z

    def __setstate__(self, st):
        cdef double x, y, z
        x, y, z = st
        self.data = dvec3(x, y, z)

    def __repr__(self):
        return f'{type(self).__name__}({fmt(self.x)}, {fmt(self.y)}, {fmt(self.z)})'

    def __len__(self):
        return 3

    def __iter__(self):
        yield self.data.x
        yield self.data.y
        yield self.data.z

    def __getitem__(self, idx):
        cdef int i
        try:
            i = idx
        except TypeError:
            if isinstance(idx, slice):
                return (self.data.x, self.data.y, self.data.z)[idx]
        else:
            if i == 0 or i == -3:
                return x(self)
            elif i == 1 or i == -2:
                return y(self)
            elif i == 2 or i == -1:
                return z(self)
        raise IndexError(idx)

    def __neg__(self):
        return vec(-self.data)

    def __pos__(self):
        return self

    def __abs__(self):
        return length(self.data)

    def __eq__(u, v):
        if isvec(u) and isvec(v):
            return (<Vec3> u).data == (<Vec3> v).data
        try:
            return tovec(u).data == tovec(v).data
        except TypeError:
            return NotImplemented

    def __add__(u, v):
        if isvec(u) and isvec(v):
            return vec((<Vec3> u).data + (<Vec3> v).data)
        try:
            return vec(tovec(u).data + tovec(v).data)
        except TypeError:
            return NotImplemented

    def __sub__(u, v):
        if isvec(u) and isvec(v):
            return vec((<Vec3> u).data - (<Vec3> v).data)
        try:
            return vec(tovec(u).data - tovec(v).data)
        except TypeError:
            return NotImplemented

    def __mul__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec((<Vec3> u).data * m)
            else:
                m = u
                return vec((<Vec3> v).data * m)
        except TypeError:
            return NotImplemented

    def __truediv__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec((<Vec3> u).data / m)
            else:
                return NotImplemented
        except TypeError:
            return NotImplemented

    def __floordiv__(u, v):
        cdef double m
        try:
            if isvec(u):
                m = v
                return vec(floor((<Vec3> u).data / m))
            else:
                return NotImplemented
        except TypeError:
            return NotImplemented

#     #
#     # Exclusive 3D functions
#     #
#     def cross(self, other):
#         """
#         The z component of the cross product between two bi-dimensional
#         smallvectors.
#         """
#         return vec(cross(self.data, tovec(other).data))
#
#     #
#     # Generic vector functions
#     #
#     def dot(self, other):
#         """Dot product with another vector."""
#
#         return dot(self.data, tovec(other).data)
#
#     def distance_to(self, other):
#         """Compute distance from vector."""
#
#         return distance(self.data, tovec(other).data)
#
#     def angle_to(self, other):
#         """Return the angle to other vector."""
#
#         return angle(self.data, tovec(other).data) * dg
#
#     def angle_to_rad(self, other):
#         """Return the angle to other vector (measured in radians)."""
#
#         return angle(self.data, tovec(other).data)
#
#     def reflect(self, direction):
#         """Reflection of vector around given normal direction."""
#
#         cdef vec *n
#         if isdirection(direction):
#             n = &(<Vec3> direction).data
#             return vec(self.data - 2 * (self.data - dot(n[0], self.data) * n[0]))
#         elif isvec(direction):
#             n = &(<Vec3> direction).data
#             return vec(self.data - 2 * (self.data - vprojection(self.data, n[0])))
#         else:
#             return self.reflect(tovec(direction))
#
#     def projection(self, direction):
#         """Projection vector for the given direction."""
#         cdef vec *n
#         if isdirection(direction):
#             n = &(<Vec3> direction).data
#             return vec(dot(self.data, n[0]) * n[0])
#         else:
#             return vec(vprojection(self.data, tovec(direction).data))
#
#     def clamp(self, double min, double max):
#         """Set length between minimum and maximum values."""
#
#         cdef double len = length(self.data)
#         if min <= len <= max:
#             return self
#         elif len < min:
#             return vec(self.data * (min / len))
#         elif len > max:
#             return vec(self.data * (max / len))
#         else:
#             raise ValueError('cannot re-scale a zero-length vector.')
#
#     def with_length(self, double size):
#         """Return a copy with the given length."""
#
#         cdef double len = length(self.data)
#         if len != 0:
#             return vec(self.data * (size / len))
#         else:
#             raise ValueError('cannot re-scale a zero-length vector.')
#
#     def lerp(self, other, double t=0.5):
#         """
#         Linear interpolation o points.
#
#         Maps the range [0, 1] to a path that goes from self to other.
#
#             t = 0  ==> self
#             t = 1  ==> other
#             others ==> a proportional combination of self and other
#         """
#
#         cdef vec vec
#         set_vec(other, &vec)
#         return vec(mix(self.data, vec, t))
#
#     def slerp(self, other, double t=0.5):
#         """
#         Spherical interpolation o points.
#
#         Like lerp, but uses polar coordinates.
#
#             t = 0  ==> self
#             t = 1  ==> other
#             others ==> a proportional combination of self and other
#
#         Notes:
#             If self and other have the same length, slerp preserves the length
#             and interpolates the angle.
#
#         See also:
#             https://en.wikipedia.org/wiki/Slerp
#         """
#         cdef vec u = self.data, v
#         set_vec(other, &v)
#
#         # Return lerp on special cases
#         ru = length(u)
#         rv = length(v)
#         if ru == 0 or rv == 0:
#             return vec(mix(u, v, t))
#
#         # Compute slerp on the happy path (here, we are always using radians)
#         sin_omega = cross(u, v) / (ru * rv)
#         omega = atan2(sin_omega, dot(u, v))
#         a = rsin((1 - t) * omega) / sin_omega
#         b = rsin(t * omega) / sin_omega
#         return vec(u * a + v * b)
#
#     def midpoint(self, other):
#         """Midpoint between two vectors."""
#
#         cdef vec vec
#         set_vec(other, &vec)
#         return vec(mix(self.data, vec, 0.5))
#
#     #
#     # Queries
#     #
#     def is_null(self):
#         """Return True if vector is null."""
#
#         return isnull(self.data)
#
#     def is_almost_null(self, tol=1e-6):
#         """Return True if vector is null within given tolerance."""
#
#         return fabs(length(self.data)) < tol
#
#     def is_normalized(self, tol=1e-6):
#         """Return True if is normalized within the given tolerance."""
#
#         return fabs(length(self.data) - 1) < tol
#
#     def is_almost_equal(self, other, double tol=1e-6):
#         """Return True if two vectors are equal under the given tolerance."""
#
#         return length(self.data - tovec(other).data) <= tol
#
#     #
#     # Norm
#     #
#     def norm(self):
#         return length(self.data)
#
#     def norm_l1(self):
#         return fabs(self.data.x) + fabs(self.data.y)
#
#     def norm_linf(self):
#         return max(fabs(self.data.x), fabs(self.data.y))
#
#     def norm_sqr(self):
#         return length2(self.data)
#
#     def normalize(self):
#         return vec(normalize(self.data))


# ------------------------------------------------------------------------------
# DIRECTION 2D CLASS
# ------------------------------------------------------------------------------

@cython.final(True)
cdef class Direction2(Vec3):
    """
    Normalized 2D vectors.
    """

    # @property
    # def length(self):
    #     return 1
    #
    # @property
    # def length_sqr(self):
    #     return 1
    #
    # def __cinit__(self, double x, double y):
    #     cdef double norm = sqrt(x * x + y * y)
    #     if norm == 0:
    #         raise ValueError('cannot initialize direction from zero-length coordinates')
    #     self.data = vec(x / norm, y / norm)
    #
    # def slerp(self, other, double t=0.5):
    #     if not isdirection(other):
    #         return super().slerp(other, t)
    #
    #     cdef vec u = self.data, v
    #     set_vec(other, &v)
    #
    #     sin_omega = cross(u, v)
    #     omega = atan2(sin_omega, dot(u, v))
    #     a = rsin((1 - t) * omega) / sin_omega
    #     b = rsin(t * omega) / sin_omega
    #     return unsafe_direction(u * a + v * b)
    #
    # def is_unity(self, tol=1e-6):
    #     return True
    #
    # def normalize(self, tol=1e-6):
    #     return self
