cimport vector_2d as V
cimport matrix_2x2 as M
from cpython cimport PyObject, PyTypeObject


cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


# ------------------------------------------------------------------------------
# VECARRAY CLASS
# ------------------------------------------------------------------------------

cdef class VecArray2:
    cdef int size
    cdef V.dvec2 *data

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

cdef VecArray2 array(int size)
cdef bint isvecarray(u)
cdef bint issized(u, int n)
cdef VecArray2 linear(VecArray2 orig, M.Mat2 mat)
