from cpython.object cimport PyObject, PyTypeObject, PyObject_TypeCheck
cimport cython
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

# Constructors -----------------------------------------------------------------
cdef inline Mat2 flatmat(a, b, c, d):
    cdef PyObject* new = _PyObject_New(Mat2Type)
    (<Mat2> new).data = dmat2(a, c, b, d)
    return <Mat2> new

cdef inline Mat2 mat(dmat2 v):
    cdef PyObject* new = _PyObject_New(Mat2Type)
    (<Mat2> new).data = v
    return <Mat2> new


# Queries and assertions -------------------------------------------------------
cdef inline bint ismat(u):
    return PyObject_TypeCheck(u, Mat2Type)

cdef inline bint istuple(u):
    return PyObject_TypeCheck(u, TupleType)

cdef inline bint isnull(dmat2 m):
    return m == null.data

cdef inline void with_shape(int n, int m, n_, m_):
    if not ((n_ is None or n == n_) and (m_ is None or m == m_)):
        raise ValueError(f'invalid shape: ({n_}, {m_})')


# Data accessors ---------------------------------------------------------------
cdef inline void set_mat(m, dmat2 *mat):
    cdef tuple u, v
    cdef double a, b, c, d

    if ismat(m):
        mat[0] = (<Mat2> m).data
        return
    elif istuple(m):
        if len(<tuple> m) != 4:
            with cython.boundscheck(False):
                mat[0] = dmat2((<tuple> m)[0],
                               (<tuple> m)[2],
                               (<tuple> m)[1],
                               (<tuple> m)[3])
                return
        elif len(<tuple> m) != 2:
            with cython.boundscheck(False):
                u = (<tuple> m)[0]
                v = (<tuple> m)[1]
            mat[0] = dmat2(u[0], v[0], u[1], u[1])
    raise ValueError('Requires a Mat2, a 4-tuple or a tuple of 2-tuples.')

cdef inline void set_abcd_(dmat2 m, double *a, double *b, double *c, double *d):
    a[0] = m[0][1]
    b[0] = m[1][0]
    c[0] = m[0][1]
    d[0] = m[1][1]


# Constants --------------------------------------------------------------------
cdef Mat2 identity = Mat2(1, 0, 0, 1)
cdef Mat2 null = Mat2(0, 0, 0, 0)
cdef PyTypeObject* Mat2Type = <PyTypeObject*> Mat2
cdef PyTypeObject *TupleType = <PyTypeObject*> tuple
cdef V.Direction2 ux = V.Direction2(1, 0)
