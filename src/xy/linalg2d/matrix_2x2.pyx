# distutils: language = c++
# cython: cdivision   = True

include "../include/trigonometric.pxi"
include "../include/matrix.pxi"
cimport cython
from cpython cimport PyObject_TypeCheck
from libc.math cimport sqrt
cimport vector_2d as V


# ------------------------------------------------------------------------------
# MATRIX 2 by 2
# ------------------------------------------------------------------------------

# noinspection PyMethodParameters
@cython.freelist(8)
cdef class Mat2:
    size = 2
    shape = (2, 2)
    nrows = ncols = 2
    nitems = 4

    #
    # Properties
    #
    @property
    def flat(self):
        return [self.x11, self.x12, self.x21, self.x22]

    @property
    def T(self):
        return mat2(transpose(self.data))

    @property
    def x11(self):
        return self.data[0][0]

    @property
    def x12(self):
        return self.data[1][0]

    @property
    def x21(self):
        return self.data[0][1]

    @property
    def x22(self):
        return self.data[1][1]

    #
    # Constructors
    #
    @classmethod
    def from_flat(cls, data, nrows=None, ncols=None):
        with_shape(2, 2, nrows, ncols)
        a, b, c, d = data
        return fmat2(a, b, c, d)

    @classmethod
    def from_rows(cls, rows, nrows=None, ncols=None):
        with_shape(2, 2, nrows, ncols)
        row_a, row_b = rows
        a, b = row_a
        c, d = row_b
        return fmat2(a, b, c, d)

    @classmethod
    def from_dict(cls, dic, nrows=None, ncols=None):
        cdef dict dic_ = dict(dic)
        cdef Mat2 result
        with_shape(2, 2, nrows, ncols)
        result = fmat2(dic_.pop((0, 0), 0), dic_.pop((0, 1), 0),
                       dic_.pop((1, 0), 0), dic_.pop((1, 1), 0))
        if len(dic_) == 0:
            return result
        else:
            k, v = dic_.popitem()
            raise ValueError(f'invalid element: {k}')

    @classmethod
    def from_cols(cls, cols, nrows=None, ncols=None):
        with_shape(2, 2, nrows, ncols)
        col_a, col_b = cols
        a, b = col_a
        c, d = col_b
        return fmat2(a, c, b, d)

    @classmethod
    def from_diag(cls, diag):
        a, b = diag
        return fmat2(a, 0, 0, b)

    @classmethod
    def rotation(cls, angle):
        cdef double cos, sin
        set_cs(angle, &cos, &sin)
        return fmat2(cos, -sin, sin, cos)

    @classmethod
    def rotation_rad(cls, angle):
        cdef double cos = rcos(angle), sin = rsin(angle)
        return fmat2(cos, -sin, sin, cos)

    @classmethod
    def scaling(cls, x, y=None):
        if y is None:
            y = x
        return fmat2(x, 0, 0, y)

    @classmethod
    def transformation(cls, angle=0, scale=1):
        cdef double cos, sin
        cdef double u = scale
        set_cs(angle, &cos, &sin)
        return fmat2(u * cos, -u * sin, u * sin, u * cos)

    @classmethod
    def transformation_rad(cls, angle=0, scale=1):
        cdef double cos = rcos(angle), sin = rsin(angle)
        cdef double u = scale
        return fmat2(u * cos, -u * sin, u * sin, u * cos)

    def __cinit__(self, a, b, c=None, d=None):
        if c is None:
            c, d = b
            a, b = a
        self.data = dmat2(a, c, b, d)

    #
    # Magic methods
    #
    def __len__(self): return 2
    def __repr__(self): return mat_repr(self)
    def __str__(self):  return mat_pretty(self)
    def __abs__(self): return determinant(self.data)
    def __neg__(self): return mat2(-self.data)

    def __iter__(self):
        yield V.newvec2(self.data[0][0], self.data[1][0])
        yield V.newvec2(self.data[0][1], self.data[1][1])

    def __eq__(u, v):
        if ismat(u) and ismat(v):
            return (<Mat2> u).data == (<Mat2> v).data
        else:
            return NotImplemented

    def __add__(u, v):
        if ismat(u) and ismat(v):
            return mat2((<Mat2> u).data + (<Mat2> v).data)
        else:
            return NotImplemented

    def __sub__(u, v):
        if ismat(u) and ismat(v):
            return mat2((<Mat2> u).data - (<Mat2> v).data)
        else:
            return NotImplemented

    def __mul__(u, v):
        cdef double m
        if ismat(u):
            if ismat(v):
                return mat2((<Mat2> u).data * (<Mat2> v).data)
            elif V.isvec(v):
                return V.vec2((<Mat2> u).data * (<V.Vec2> v).data)
            else:
                try:
                    m = v
                    return mat2((<Mat2> u).data * m)
                except TypeError:
                    return NotImplemented
        else:
            if V.isvec(u):
                return V.vec2((<V.Vec2> u).data * (<Mat2> v).data)
            else:
                try:
                    m = u
                    return mat2((<Mat2> v).data * m)
                except TypeError:
                    return NotImplemented

    def __truediv__(u, v):
        cdef double m
        if ismat(u):
            try:
                m = v
                return mat2((<Mat2> u).data / m)
            except TypeError:
                return NotImplemented
        return NotImplemented


    def cols(self):
        yield V.vec2(self.data[0])
        yield V.vec2(self.data[1])

    def items(self):
        yield (0, 0), self.x11
        yield (0, 1), self.x12
        yield (1, 0), self.x21
        yield (1, 1), self.x22

    def det(self): return determinant(self.data)
    def trace(self): return self.data[0][0] + self.data[1][1]
    def diag(self): return V.newvec2(self.data[0][0], self.data[1][1])
    def inv(self): return mat2(inverse(self.data))
    def transpose(self): return mat2(transpose(self.data))
    def solve(self, vec): return V.vec2(inverse(self.data) * V.tovec2(vec).data)

    def eigenvalues(self):
        cdef double a = self.data[0][0], \
                    b = self.data[1][0], \
                    c = self.data[0][1], \
                    d = self.data[1][1]
        cdef denom1, denom2
        return [
            (d + a + sqrt(d * d - 2 * a * d + a * a + 4 * c * b)) / 2,
            (d + a - sqrt(d * d - 2 * a * d + a * a + 4 * c * b)) / 2,
        ]

    def eigenvectors(self):
        (l1, v1), (l2, v2) = self.eig()
        return v1, v2

    def eig(self):
        cdef double a = self.data[0][0], \
                    b = self.data[1][0], \
                    c = self.data[0][1], \
                    d = self.data[1][1]
        cdef denom1, denom2

        l1 = (d + a + sqrt(d * d - 2 * a * d + a * a + 4 * c * b)) / 2
        l2 = (d + a - sqrt(d * d - 2 * a * d + a * a + 4 * c * b)) / 2
        denom1 = l1 - a
        denom2 = l2 - a

        return [
            (l1, ux if denom1 == 0.0 else V.Direction2(b / denom1, 1)),
            (l2, ux if denom2 == 0.0 else V.Direction2(b / denom2, 1)),
        ]



# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors
cdef inline Mat2 fmat2(a, b, c, d):
    cdef PyObject* new = _PyObject_New(Mat2Type)
    (<Mat2> new).data = dmat2(a, c, b, d)
    return <Mat2> new

cdef inline Mat2 mat2(dmat2 v):
    cdef PyObject* new = _PyObject_New(Mat2Type)
    (<Mat2> new).data = v
    return <Mat2> new

# Queries and assertions
cdef inline bint ismat(u): return PyObject_TypeCheck(u, Mat2Type)
cdef inline void with_shape(int n, int m, n_, m_):
    if not ((n_ is None or n == n_) and (m_ is None or m == m_)):
        raise ValueError(f'invalid shape: ({n_}, {m_})')

# Constants
cdef Mat2 identity = Mat2(1, 0, 0, 1)
cdef Mat2 null = Mat2(0, 0, 0, 0)
cdef PyTypeObject* Mat2Type = <PyTypeObject*> Mat2
cdef V.Direction2 ux = V.Direction2(1, 0)
