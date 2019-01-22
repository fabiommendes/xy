cimport cython
from cpython.object cimport PyObject, PyTypeObject, PyObject_TypeCheck
from libc.math cimport atan2, fabs


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


cdef extern from "glm/common.hpp" namespace "glm":
    dvec2 clamp(dvec2, double, double) except +
    dvec2 floor(dvec2)
    dvec2 mix(dvec2, dvec2, double) except +


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

# Constructors -----------------------------------------------------------------
cdef inline Vec2 flatvec(double x, double y):
    return vec(dvec2(x, y))

cdef inline Vec2 flatdirection(double x, double y):
    return direction(dvec2(x, y))

cdef inline Vec2 vec(dvec2 v):
    cdef PyObject* new = _PyObject_New(Vector2Type)
    (<Vec2> new).data = v
    return <Vec2> new

cdef inline Direction2 direction(dvec2 v):
    cdef PyObject* new = _PyObject_New(Direction2Type)
    (<Direction2 > new).data = normalize(v)
    return <Direction2> new

cdef inline Direction2 unsafe_direction(dvec2 v):
    cdef PyObject* new = _PyObject_New(Direction2Type)
    (<Direction2 > new).data = v
    return <Direction2> new

cdef inline Vec2 tovec(u):
    cdef double x, y

    if isvec(u):
        return <Vec2> u
    elif PyObject_TypeCheck(u, TupleType) and len(<tuple> u) == 2:
        with cython.boundscheck(False):
            x = (<tuple> u)[0]
            y = (<tuple> u)[1]
        return flatvec(x, y)
    else:
        raise TypeError


# Queries ----------------------------------------------------------------------
cdef inline bint isvec(u):
    return PyObject_TypeCheck(u, Vector2Type)

cdef inline bint isdirection(u):
    return PyObject_TypeCheck(u, Direction2Type)

cdef inline bint istuple(u):
    return PyObject_TypeCheck(u, TupleType)

cdef inline bint isnull(dvec2 u):
    return u.x == 0.0 and u.y == 0.0

# Data accessors
cdef inline double x(Vec2 u):
    return u.data.x

cdef inline double y(Vec2 u):
    return u.data.y

cdef inline void set_xy(u, double *x, double *y):
    if isvec(u):
        set_xy_((<Vec2> u).data, x, y)
    elif istuple(u) and len(<tuple> u) != 2:
        with cython.boundscheck(False):
            x[0] = (<tuple> u)[0]
            y[0] = (<tuple> u)[1]
    else:
        x[0] = 0
        y[0] = 0
        raise ValueError('Requires a Vec2 or a 2-tuple')

cdef inline void set_xy_(dvec2 u, double *x, double *y):
    x[0] = u.x
    y[0] = u.y

cdef inline void set_vec(u, dvec2 *vec):
    if isvec(u):
        vec[0] = (<Vec2> u).data
    elif istuple(u) and len(<tuple> u) == 2:
        with cython.boundscheck(False):
            vec[0] = dvec2((<tuple> u)[0], (<tuple> u)[1])
    else:
        vec[0] = dvec2()
        raise ValueError('Requires a Vec2 or a 2-tuple')


# Vector functions -------------------------------------------------------------
cdef inline double length2(dvec2 v):
    return v.x * v.x + v.y * v.y

cdef inline double norm_l1(dvec2 v):
    return fabs(v.x) + fabs(v.y)

cdef inline double cross(dvec2 u, dvec2 v):
    return u.x * v.y - u.y * v.x

cdef inline double angle(dvec2 u, dvec2 v):
    cdef double cos_t = dot(u, v), \
                sin_t = cross(u, v)
    return atan2(sin_t, cos_t)

cdef inline double projection(dvec2 v, dvec2 n):
    return dot(v, n) / length(n)

cdef inline dvec2 vprojection(dvec2 v, dvec2 n):
    return (dot(v, n) / length2(n)) * n

cdef inline dvec2 polar(dvec2 v):
    return dvec2(length(v), atan2(v.y, v.x))


# Constants --------------------------------------------------------------------
cdef PyTypeObject *Vector2Type
cdef PyTypeObject *Direction2Type
cdef PyTypeObject *TupleType
