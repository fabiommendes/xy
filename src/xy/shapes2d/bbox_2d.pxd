from cpython cimport PyObject, PyTypeObject


cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


# ------------------------------------------------------------------------------
# RECTANGULAR BOUNDING BOX 2 by 2
# ------------------------------------------------------------------------------

cdef class BBox2:
    cdef readonly double xmin, xmax, ymin, ymax


# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

cdef BBox2 bbox2(double a, double b, double c, double d)
cdef bint isbbox(u)
cdef bint eq(BBox2 a, BBox2 b)
cdef check_empty(dict d)
