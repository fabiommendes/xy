from cpython.object cimport PyObject, PyTypeObject


cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


cdef extern from "glm/ext/vector_double2.hpp" namespace "glm":
    cdef cppclass dvec2:
        double x
        double y
        dvec2()
        dvec2(double, double)
        dvec2 operator +(dvec2, dvec2)
        dvec2 operator -(dvec2, dvec2)
        dvec2 operator *(dvec2, double)
        dvec2 operator *(double, dvec2)
        dvec2 operator /(dvec2, double)
        dvec2 operator -()
        bint operator ==(dvec2, dvec2)
        double operator [](int)


cdef extern from "glm/geometric.hpp" namespace "glm":
    double length(const dvec2)
    double dot(const dvec2, const dvec2)
    double distance(const dvec2, const dvec2)
    dvec2 faceforward(const dvec2, const dvec2, const dvec2)
    dvec2 normalize(const dvec2)
    dvec2 reflect(const dvec2, const dvec2)
    dvec2 refract(const dvec2, const dvec2, double)


# ------------------------------------------------------------------------------
# VECTOR CLASSES
# ------------------------------------------------------------------------------

cdef class Vec2:
    cdef dvec2 data


cdef class Direction2(Vec2):
    pass

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

cdef Vec2 vec2(dvec2)
cdef Vec2 direction2(dvec2)
cdef Vec2 newvec2(double, double)
cdef Vec2 tovec2(u)
cdef bint isvec(u)
cdef void set_xy(u, double *x, double *y)
cdef void set_vec(u, dvec2 *vec)
cdef double length2(dvec2)
