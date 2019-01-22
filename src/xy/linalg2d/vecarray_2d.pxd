cimport vector_2d as V
cimport matrix_2x2 as M
from cpython.mem cimport  PyMem_Free, PyMem_Malloc
from cpython.object cimport PyObject, PyTypeObject, PyObject_TypeCheck


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

# Constructors -----------------------------------------------------------------
cdef inline VecArray2 array(int size):
    cdef PyObject* new = _PyObject_New(VecArray2Type)
    (<VecArray2> new).size = size
    (<VecArray2> new).data = <V.dvec2*> PyMem_Malloc(size * sizeof(V.dvec2))
    return <VecArray2> new


# Queries and assertions -------------------------------------------------------
cdef inline bint isvecarray(u):
    return PyObject_TypeCheck(u, VecArray2Type)

cdef inline bint issized(u, int n):
    return PyObject_TypeCheck(u, VecArray2Type) and (<VecArray2> u).size == n


# Transformations --------------------------------------------------------------
cdef inline VecArray2 linear(VecArray2 orig, M.Mat2 mat):
    cdef VecArray2 new = array(orig.size)
    for i in range(orig.size):
        new.data[i] = mat.data * orig.data[i]
    return new


# Constants --------------------------------------------------------------------
cdef PyTypeObject* VecArray2Type = <PyTypeObject*> VecArray2
