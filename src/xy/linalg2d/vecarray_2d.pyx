# distutils: language = c++
# cython: cdivision   = True

cimport vector_2d as V
cimport matrix_2x2 as M
from cpython cimport PyObject_TypeCheck, PyMem_Free, PyMem_Malloc


# ------------------------------------------------------------------------------
# VECARRAY CLASS
# ------------------------------------------------------------------------------

# noinspection PyMethodParameters
cdef class VecArray2:
    """
    An array of bi-dimensional vectors.
    """

    def __cinit__(self, data):
        self.size = len(data)
        self.data = <V.dvec2*> PyMem_Malloc(self.size * sizeof(V.dvec2))

    def __dealloc__(self):
        PyMem_Free(self.data)

    #
    # Magic methods
    #
    def __len__(self): return self.size
    def __str__(self): return repr(self)
    def __bool__(self): return self.size != 0

    def __repr__(self):
        tname = type(self).__name__
        data = ', '.join([str(tuple(v)) for v in self])
        return '%s([%s])' % (tname, data)

    def __iter__(self):
        for i in range(self.size):
            yield V.vec2(self.data[i])

    def __getitem__(self, i):
        cdef int idx
        if isinstance(i, int):
            idx = i
            if idx >= self.size or idx < -self.size:
                raise IndexError(i)
            elif idx < 0:
                return V.vec2(self.data[self.size - i])
            else:
                return V.vec2(self.data[i])
        raise TypeError

    def __add__(u, v):
        try:
            return (<VecArray2> u).translate(v) if isvecarray(u) else (<VecArray2> v).translate(u)
        except TypeError:
            return NotImplementedError

    def __sub__(u, v):
        try:
            return (<VecArray2> u).translate(-v) if isvecarray(u) else (<VecArray2> v).translate(-u)
        except TypeError:
            return NotImplementedError

    def __mul__(u, v):
        try:
            return (<VecArray2> u).scale(v) if isvecarray(u) else (<VecArray2> v).scale(u)
        except TypeError:
            return NotImplementedError

    def __truediv__(u, v):
        try:
            return (<VecArray2> u).scale(1/v) if isvecarray(u) else (<VecArray2> v).scale(1/u)
        except TypeError:
            return NotImplementedError

    #
    # Properties
    #
    def norm(self):
        """Return Array with the norm of each vector."""

    def norm_sqr(self):
        """Return Array with the squared norm of each vector."""

    #
    # Transformations
    #
    def normalize(self):
        """Return a VecArray of normalized vectors."""
        return VecArray2([u.normalize() for u in self])

    def linear_transform(self, mat):
        """Apply matrix in each element of the array."""
        return linear(self, mat)

    def rotate(self, double angle):
        """Rotate all vectors by the given angle. (angle in degrees)"""
        return linear(self, M.Mat2.rotation(angle))

    def rotate_rad(self, double angle):
        """Rotate all vectors by the given angle. (angle in degrees)"""
        return linear(self, M.Mat2.rotation_rad(angle))

    def rotate_around(self, double angle, axis):
        """Rotate all vectors by the given angle. (angle in degrees)"""
        raise NotImplementedError

    def rotate_around_rad(self, double angle, axis):
        """Rotate all vectors by the given angle. (angle in degrees)"""
        raise NotImplementedError

    def translate(self, delta):
        """Translate all vectors by the given factor."""
        cdef VecArray2 new = array(self.size)
        cdef V.dvec2 u
        V.set_vec(delta, &u)

        for i in range(self.size):
            new.data[i] = self.data[i] + u
        return new

    def scale(self, double delta):
        """Re-scale all vectors by the given factor."""
        cdef VecArray2 new = array(self.size)
        for i in range(self.size):
            new.data[i] = self.data[i] * delta
        return new

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors
cdef inline VecArray2 array(int size):
    cdef PyObject* new = _PyObject_New(VecArray2Type)
    (<VecArray2> new).size = size
    (<VecArray2> new).data = <V.dvec2*> PyMem_Malloc(size * sizeof(V.dvec2))
    return <VecArray2> new

# Queries and assertions
cdef inline bint isvecarray(u): return PyObject_TypeCheck(u, VecArray2Type)
cdef inline bint issized(u, int n): return PyObject_TypeCheck(u, VecArray2Type) and (<VecArray2> u).size == n

# Transformations
cdef inline VecArray2 linear(VecArray2 orig, M.Mat2 mat):
    cdef VecArray2 new = array(orig.size)
    for i in range(orig.size):
        new.data[i] = mat.data * orig.data[i]
    return new

# Constants
cdef PyTypeObject* VecArray2Type = <PyTypeObject*> VecArray2

