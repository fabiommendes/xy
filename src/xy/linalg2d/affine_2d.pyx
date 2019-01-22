# distutils: language = c++
# cython: cdivision   = True

include "../include/trigonometric.pxi"
include "../include/text.pxi"
include "../include/matrix.pxi"
cimport cython
from libc.math cimport atan2, fabs
from cpython.object cimport PyObject_TypeCheck
cimport vector_2d as V


# ------------------------------------------------------------------------------
# AFFINITY TRANSFORM CLASS
# ------------------------------------------------------------------------------

# noinspection PyMethodParameters
@cython.freelist(8)
cdef class Affine2:
    """
    A 2D affinity transformation.

    Affinity transformations are linear transformations expanded with a
    displacement.
    """
    ndim = 2
    size = 3
    dtype = float
    shape = (3, 2)

    #
    # Properties and accessors
    #
    @property
    def displacement(self):
        return V.vec(self.cvector)

    @property
    def linear(self):
        return M.mat(self.cmatrix)

    @property
    def augmented_matrix(self):
        cdef double a, b, c, d, x, y
        M.set_abcd_(self.cmatrix, &a, &b, &c, &d)
        V.set_xy_(self.cvector, &x, &y)
        return _mat_funcs[3]((a, b, x), (c, d, y), (0., 0., 1.))

    @property
    def rectangular_matrix(self):
        cdef double a, b, c, d, x, y
        M.set_abcd_(self.cmatrix, &a, &b, &c, &d)
        V.set_xy_(self.cvector, &x, &y)
        return _mat_funcs[0]([[a, b, x], [c, d, y]])

    #
    # Constructors
    #
    @classmethod
    def from_flat(cls, data):
        """Create vector from a 6-sequence."""
        a, b, c, d, x, y = data
        return cls.from_coords(a, b, c, d, x, y)

    @classmethod
    def from_rows(cls, rows):
        row1, row2 = rows
        a, b, x = row1
        c, d, y = row2
        return cls.from_coords(a, b, c, d, x, y)

    @classmethod
    def from_coords(cls, double a, double b, double c, double d, double x, double y):
        return affine(M.dmat2(a, c, b, d), V.dvec2(x, y))

    def __cinit__(self, m, v=(0, 0)):
        V.set_vec(v, &self.cvector)
        M.set_mat(m, &self.cmatrix)

    #
    # Generic API
    #
    def copy(self, matrix=None, vector=None):
        """
        Creates a new copy overriding the matrix or vector part of
        transformation.
        """
        raise NotImplementedError

    #
    # Magic methods
    #
    def __getstate__(self):
        return *self.displacement, *self.linear

    def __setstate__(self, tuple st):
        V.set_vec(st[:2], &self.cvector)
        M.set_mat(st[2:], &self.cmatrix)

    def __repr__(self):
        return f'{type(self).__name__}({self.linear}, {self.displacement})'

    def __str__(self):
        return affine_pretty(self.linear, self.displacement)

    def __neg__(self):
        return affine(self.cmatrix, -self.cvector)

    def __pos__(self):
        return self

    def __eq__(u, v):
        if isaffine(u):
            return ((<Affine2> u).cmatrix == (<Affine2> v).cmatrix and
                    (<Affine2> u).cvector == (<Affine2> v).cvector)
        else:
            return NotImplemented

    def __call__(self, v):
        cdef V.dvec2 vec
        V.set_vec(v, &vec)
        return V.vec(self.cmatrix * vec + self.cvector)

    def __add__(u, v):
        cdef V.dvec2 vec
        cdef Affine2 aff

        if isaffine(u):
            V.set_vec(v, &vec)
            aff = <Affine2> u
        elif isaffine(v):
            V.set_vec(u, &vec)
            aff = <Affine2> v
        else:
            return NotImplemented
        return affine(aff.cmatrix, aff.cvector + vec)

    def __sub__(u, v):
        cdef V.dvec2 vec

        if isaffine(u):
            V.set_vec(v, &vec)
            return affine((<Affine2> u).cmatrix, (<Affine2> u).cvector - vec)
        elif isaffine(v):
            V.set_vec(u, &vec)
            return affine((<Affine2> v).cmatrix, vec - (<Affine2> v).cvector)
        else:
            return NotImplemented

    def __mul__(u, v):
        cdef double m
        cdef M.dmat2 mat

        if isaffine(u):
            if M.ismat(v):
                M.set_mat(v, &mat)
                return affine((<Affine2> u).cmatrix * mat, (<Affine2> u).cvector)
            else:
                try:
                    m = v
                except TypeError:
                    return NotImplemented
                else:
                    return affine((<Affine2> u).cmatrix * m, (<Affine2> u).cvector)

        return NotImplemented

    #
    # Transform elements
    #
    def transform_vector(self, u):
        cdef V.dvec2 vec
        V.set_vec(u, &vec)
        return V.vec(self.cmatrix * vec + self.cvector)

    def transform_point(self, double x, double y):
        cdef V.dvec2 vec = V.dvec2(x, y)
        return V.vec(self.cmatrix * vec + self.cvector)

    #
    # Modify affine transform
    #
    def translate(self, u):
        """Translate by given vector."""

        cdef V.dvec2 vec
        V.set_vec(u, &vec)
        return affine(self.cmatrix, self.cvector + vec)

    def translate_xy(self, double x, double y):
        """Translate by given x, y coordinates."""

        return affine(self.cmatrix, self.cvector + V.dvec2(x, y))

    def rotate(self, double angle):
        """Rotate around origin."""

        raise NotImplementedError

    def rotate_around(self, double angle, axis):
        """Rotate around given axis."""

        raise NotImplementedError

    def rotate_rad(self, double angle):
        """Rotate around origin (angle in radians)."""

        raise NotImplementedError

    def rotate_around_rad(self, double angle, axis):
        """Rotate around origin (angle in radians)."""

        raise NotImplementedError

    def perpendicular(self, ccw=True):
        """
        Return rotate 90 degrees to the left (or right if ccw=False).
        """
        raise NotImplementedError

    #
    # Queries
    #
    def is_null(self):
        """Return True if transform have only null constants."""

        return V.isnull(self.cvector) and M.isnull(self.cmatrix)

    def is_almost_null(self, tol=1e-6):
        """Return True if transform is null within given tolerance."""

        raise NotImplementedError

    def is_similarity(self, tol=1e-6):
        """Return True if can be interpreted as a similarity transform within
        the given tolerance."""

        raise NotImplementedError

    def is_linear(self, tol=1e-6):
        """Return True if can be interpreted as a linear transformation within
        the given tolerance."""

        return fabs(V.length(self.cvector) - 1) < tol

    def is_almost_equal(self, other, double tol=1e-6):
        """Return True if two transforms are equal within the given tolerance."""

        raise NotImplementedError


# ------------------------------------------------------------------------------
# Similarity transform
# ------------------------------------------------------------------------------

cdef class Similarity2(Affine2):
    """
    Similarity transformation in 2 dimensions.

    Represents shape preserving transformations like rotations, changes of
    scale, translations and reflections.
    """

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
