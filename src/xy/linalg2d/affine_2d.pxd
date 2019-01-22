from cpython.object cimport PyObject, PyTypeObject, PyObject_TypeCheck
cimport vector_2d as V
cimport matrix_2x2 as M


cdef extern from "Python.h":
    PyObject* _PyObject_New(PyTypeObject*)


# ------------------------------------------------------------------------------
# TRANSFORMATION CLASSES
# ------------------------------------------------------------------------------

cdef class Affine2:
    cdef M.dmat2 cmatrix
    cdef V.dvec2 cvector


cdef class Similarity2(Affine2):
    pass


# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors -----------------------------------------------------------------
cdef inline Affine2 affine(M.dmat2 mat, V.dvec2 v):
    cdef PyObject* new = _PyObject_New(Affine2Type)
    (<Affine2> new).cmatrix = mat
    (<Affine2> new).cvector = v
    return <Affine2> new

cdef inline Similarity2 similarity(M.dmat2 mat, V.dvec2 v):
    cdef PyObject* new = _PyObject_New(Similarity2Type)
    (<Similarity2> new).cvector = v
    (<Similarity2> new).cmatrix = mat
    return <Similarity2> new


# Queries ----------------------------------------------------------------------
cdef inline bint issimilarity(u):
    return PyObject_TypeCheck(u, Affine2Type)

cdef inline bint isaffine(u):
    return PyObject_TypeCheck(u, Similarity2Type)


# Data accessors ---------------------------------------------------------------
cdef inline V.Vec2 vector(Affine2 t):
    return V.vec(t.cvector)

cdef inline M.Mat2 matrix(Affine2 t):
    return M.mat(t.cmatrix)

cdef inline void set_matvec(u, M.dmat2 *m, V.dvec2 *v):
    if isaffine(u):
        m[0] = (<Affine2> u).cmatrix
        v[0] = (<Affine2> u).cvector
    else:
        raise ValueError('Requires a affinity transform.')


# Transformation functions -----------------------------------------------------
cdef inline V.dvec2 transform_vector(Affine2 t, V.dvec2 v):
    return t.cmatrix * v + t.cvector


# Constants --------------------------------------------------------------------
cdef Affine2 null = Affine2()
cdef PyTypeObject *Affine2Type = <PyTypeObject*> Affine2
cdef PyTypeObject *Similarity2Type = <PyTypeObject*> Similarity2
cdef PyTypeObject *TupleType = <PyTypeObject*> tuple
_mat_funcs = [None, M.Mat2, None, None]