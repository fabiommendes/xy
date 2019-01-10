from cpython cimport PyObject, PyTypeObject
cimport vector_2d as V


cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


cdef extern from "glm/ext/matrix_double2x2.hpp" namespace "glm":
    cdef cppclass dmat2:
        dmat2()
        dmat2(double, double, double, double)
        dmat2(V.dvec2, V.dvec2)
        dmat2 operator +(dmat2, dmat2)
        dmat2 operator -(dmat2, dmat2)
        dmat2 operator *(dmat2, double)
        dmat2 operator *(double, dmat2)
        V.dvec2 operator *(dmat2, V.dvec2)
        V.dvec2  operator *(V.dvec2, dmat2)
        dmat2 operator *(dmat2, dmat2)
        dmat2 operator /(dmat2, double)
        dmat2 operator -()
        bint operator ==(dmat2, dmat2)
        V.dvec2 operator [](int)
    dmat2 transpose(const dmat2)
    dmat2 inverse(const dmat2)
    double determinant(const dmat2)


# ------------------------------------------------------------------------------
# MATRIX 2 by 2
# ------------------------------------------------------------------------------

cdef class Mat2:
    cdef dmat2 data


# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

cdef Mat2 fmat2(a, b, c, d)
cdef Mat2 mat2(dmat2 v)
cdef bint ismat(u)