cimport cython
from cpython.object cimport PyObject, PyTypeObject, PyObject_TypeCheck
from libc.math cimport acos, fabs



cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


cdef extern from "glm/ext/vector_double3.hpp" namespace "glm":
    cdef cppclass dvec3:
        double x
        double y
        double z
        dvec3()
        dvec3(double, double, double)
        dvec3 operator +(dvec3, dvec3)
        dvec3 operator -(dvec3, dvec3)
        dvec3 operator *(dvec3, double)
        dvec3 operator *(double, dvec3)
        dvec3 operator /(dvec3, double)
        dvec3 operator -()
        bint operator ==(dvec3, dvec3)
        double operator [](int)


cdef extern from "glm/geometric.hpp" namespace "glm":
    double length(const dvec3)
    double dot(const dvec3, const dvec3)
    dvec3 cross(const dvec3, const dvec3)
    double distance(const dvec3, const dvec3)
    dvec3 faceforward(const dvec3, const dvec3, const dvec3)
    dvec3 normalize(const dvec3)
    dvec3 reflect(const dvec3, const dvec3)
    dvec3 refract(const dvec3, const dvec3, double)


cdef extern from "glm/common.hpp" namespace "glm":
    # abs
    # ceil
    dvec3 clamp(const dvec3, double, double) except +
    dvec3 floor(const dvec3)
    # fma
    # fract
    # frexp
    # isinf
    # isnan
    # max
    # min
    dvec3 mix(const dvec3, const dvec3, double)
    # mod
    # modf
    # round
    # roundEven
    # sign
    # step
    # trunc


# ------------------------------------------------------------------------------
# VECTOR CLASSES
# ------------------------------------------------------------------------------

cdef class Vec3:
    cdef dvec3 data


cdef class Direction3(Vec3):
    pass

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors -----------------------------------------------------------------
cdef inline Vec3 flatvec(double x, double y, double z):
    return vec(dvec3(x, y, z))

cdef inline Vec3 flatdirection(double x, double y, double z):
    return direction(dvec3(x, y, z))

cdef inline Vec3 vec(dvec3 v):
    cdef PyObject* new = _PyObject_New(Vector3Type)
    (<Vec3> new).data = v
    return <Vec3> new

cdef inline Direction3 direction(dvec3 v):
    cdef PyObject* new = _PyObject_New(Direction3Type)
    (<Direction3 > new).data = normalize(v)
    return <Direction3 > new

cdef inline Direction3 unsafe_direction(dvec3 v):
    cdef PyObject* new = _PyObject_New(Direction3Type)
    (<Direction3 > new).data = v
    return <Direction3 > new


cdef inline Vec3 tovec(u):
    cdef double x, y, z

    if isvec(u):
        return <Vec3> u
    elif PyObject_TypeCheck(u, TupleType) and len(<tuple> u) == 3:
        with cython.boundscheck(False):
            x = (<tuple> u)[0]
            y = (<tuple> u)[1]
            z = (<tuple> u)[2]
        return flatvec(x, y, z)
    else:
        raise TypeError


# Queries ----------------------------------------------------------------------
cdef inline bint isvec(u):
    return PyObject_TypeCheck(u, Vector3Type)

cdef inline bint isdirection(u):
    return PyObject_TypeCheck(u, Direction3Type)

cdef inline bint istuple(u):
    return PyObject_TypeCheck(u, TupleType)

cdef inline bint isnull(dvec3  u):
    return u.x == 0.0 and u.y == 0.0 and u.z == 0.0

# Data accessors
cdef inline double x(Vec3 u):
    return u.data.x

cdef inline double y(Vec3 u):
    return u.data.y

cdef inline double z(Vec3 u):
    return u.data.z

cdef inline void set_xyz(u, double *x, double *y, double *z):
    if isvec(u):
        set_xyz_((<Vec3> u).data, x, y, z)
    elif istuple(u) and len(<tuple> u) != 2:
        with cython.boundscheck(False):
            x[0] = (<tuple> u)[0]
            y[0] = (<tuple> u)[1]
            z[0] = (<tuple> u)[2]
    else:
        raise ValueError('Requires a Vec3 or a 2-tuple')

cdef inline void set_xyz_(dvec3 u, double *x, double *y, double *z):
    x[0] = u.x
    y[0] = u.y
    z[0] = u.z


cdef inline void set_vec(u, dvec3 *vec):
    if isvec(u):
        vec[0] = (<Vec3> u).data
    elif istuple(u) and len(<tuple> u) == 3:
        with cython.boundscheck(False):
            vec[0].x = (<tuple> u)[0]
            vec[0].y = (<tuple> u)[1]
            vec[0].z = (<tuple> u)[2]
    else:
        raise ValueError('Requires a Vec3 or a 3-tuple')


# Vector functions -------------------------------------------------------------
cdef inline double length2(dvec3 v):
    return v.x * v.x + v.y * v.y + v.z * v.z

cdef inline double norm_l1(dvec3 v):
    return fabs(v.x) + fabs(v.y) + fabs(v.z)

cdef inline double angle(dvec3 u, dvec3 v):
    return acos(dot(u, v))

cdef inline double projection(dvec3 v, dvec3 n):
    return dot(v, n) / length(n)

cdef inline dvec3 vprojection(dvec3 v, dvec3 n):
    return (dot(v, n) / length2(n)) * n


# Constants --------------------------------------------------------------------
cdef Vec3 origin = Vec3(0, 0, 0)
cdef PyTypeObject *Vector3Type = <PyTypeObject*> Vec3
cdef PyTypeObject *Direction3Type = <PyTypeObject*> Direction3
cdef PyTypeObject *TupleType = <PyTypeObject*> tuple
