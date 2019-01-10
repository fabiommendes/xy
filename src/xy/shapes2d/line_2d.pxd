from cpython cimport PyObject, PyTypeObject
cimport xy.linalg2d.vector_2d as V


cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


# ------------------------------------------------------------------------------
# CLASSES
# ------------------------------------------------------------------------------

cdef class Segment2:
    cdef V.dvec2 _start, _end
    cdef readonly double radius


cdef class Line2:
    cdef readonly double a, b, c, radius


cdef class Ray2:
    cdef V.dvec2 _start, _direction
    cdef readonly double radius


# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

cdef Segment2 segment(V.dvec2 u, V.dvec2 v, double radius=*)
cdef Line2 line(double a, double b, double c, double radius=*)
cdef bint issegment(u)
cdef check_empty(dict d)