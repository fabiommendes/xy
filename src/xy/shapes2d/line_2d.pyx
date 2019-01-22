# distutils: language = c++
# cython: cdivision = True

cimport cython
from cpython cimport PyObject_TypeCheck


# ------------------------------------------------------------------------------
# LINE SEGMENT
# ------------------------------------------------------------------------------

@cython.freelist(8)
cdef class Segment2:
    """
    A finite directed 2D line segment.
    """

    @property
    def x1(self):
        return self._start.x

    @property
    def x2(self):
        return self._end.x

    @property
    def y1(self):
        return self._start.y

    @property
    def y2(self):
        return self._end.y

    @property
    def start(self):
        return V.vec(self._start)

    @property
    def end(self):
        return V.vec(self._end)

    @property
    def displacement(self):
        return V.vec(self._end - self._start)

    @property
    def direction(self):
        return V.direction(self._end - self._start)

    @property
    def midpoint(self):
        return V.vec(self._end * 0.5 + self._start * 0.5)

    @property
    def bisection(self):
        raise NotImplementedError

    @property
    def length(self):
        return V.length(self._end - self._start)

    @property
    def length_sqr(self):
        return V.length2(self._end - self._start)

    #
    # Constructors
    #
    def __cinit__(self, start, end, radius=0.0):
        V.set_vec(start, &self._start)
        V.set_vec(end, &self._end)
        self.radius = radius

    #
    # Geometric transformations
    #
    def rotate(self, double angle):
        """Rotate segment by angle (in degrees)"""
        raise NotImplementedError

    def rotate_rad(self, double angle):
        """Rotate segment by angle (in radians)"""
        raise NotImplementedError

    def rotate_around(self, double angle, axis):
        """Rotate segment by angle (in degrees)"""
        raise NotImplementedError

    def rotate_around_rad(self, double angle, axis):
        """Rotate segment by angle (in radians)"""
        raise NotImplementedError

    def translate(self, delta):
        """Displace segment by given factor"""
        cdef V.dvec2 vec
        V.set_vec(delta, &vec)
        return segment(self._start + vec, self._end + vec, self.radius)

    def scale(self, double factor):
        """Scale around origin"""
        raise NotImplementedError

    def scale_around(self, double factor, axis):
        """Scale around origin"""
        raise NotImplementedError

    #
    # Geometric properties
    #
    def distance_to(self, double x, double y):
        """Return distance from line to given point without considering
        radius."""
        raise NotImplementedError


# ------------------------------------------------------------------------------
# INFINITE LINE
# ------------------------------------------------------------------------------

@cython.freelist(8)
cdef class Line2:
    """
    Infinite line with given slope and intercept Point.

    It obeys the equation a * x + b * y = c
    """

    @property
    def slope(self):
        return self.b / self.a

    @property
    def intercept_x(self):
        return self.c / self.a

    @property
    def intercept_y(self):
        return self.c / self.b

    #
    # Constructors
    #
    def __cinit__(self, double a, double b, double c = 0, radius=0):
        self.a = a
        self.b = b
        self.c = c
        self.radius = radius

    #
    # Magic methods
    #
    def __getstate__(self):
        return self.a, self.b, self.c, self.radius

    def __setstate__(self, st):
        a, b, c, radius = st
        self.__cinit__(a, b, c, radius)

    #
    # Geometric transformations
    #
    def rotate(self, double angle):
        """Rotate line by angle (in degrees)"""
        raise NotImplementedError

    def rotate_around(self, double angle, axis):
        """Rotate line by angle (in degrees)"""
        raise NotImplementedError

    def translate(self, delta):
        """Translate line by given displacement vector"""
        cdef double x, y
        V.set_xy(delta, &x, &y)
        return line(self.a, self.b, self.c + self.a * x + self.b * y)

    def scale(self, double factor):
        """Rescale line (only affects radius)"""
        return line(self.a, self.b, self.c, self.radius * factor)

    def scale_around(self, double factor, axis):
        """Rescale line (only affects radius)"""
        raise NotImplemented

    def flip_x(self):
        return line(-self.a, self.b, self.c, self.radius)

    def flip_y(self):
        return line(self.a, -self.b, self.c, self.radius)


    #
    # Geometric properties
    #
    def distance_to(self, pt):
        """Return distance from line to given point without considering
        radius."""
        raise NotImplementedError

    #
    # API
    #
    def get_x(self, double y):
        """Return the x-coordinate corresponding to a given y value."""
        return (self.c - self.b * y) / self.a

    def get_y(self, double x):
        """Return the y-coordinate corresponding to a given y value."""
        return (self.c - self.a * x) / self.b

    # Queries
    def is_left_point(self, pt):
        """
        Return True if point is on the left region of the line.
        """
        raise NotImplementedError

    def is_right_point(self, pt):
        """
        Return True if point is on the right region of the line.
        """
        raise NotImplementedError

    def is_inner_point(self, pt):
        """
        Return True if point is on the inner region of the line.
        """
        return self.distance_to(pt) <= self.radius


# ------------------------------------------------------------------------------
# DIRECTED RAY (SEMI-INFINITE LINE)
# ------------------------------------------------------------------------------

@cython.freelist(8)
cdef class Ray2:
    """
    A directed semi-infinite line.
    """


# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors
cdef inline Segment2 segment(V.dvec2 u, V.dvec2 v, double radius=0):
    cdef PyObject* new = _PyObject_New(Segment2Type)
    (<Segment2> new)._start = u
    (<Segment2> new)._end = v
    (<Segment2> new).radius = radius
    return <Segment2> new

cdef inline Line2 line(double a, double b, double c, double radius=0):
    cdef PyObject* new = _PyObject_New(Line2Type)
    (<Line2> new).a = a
    (<Line2> new).b = b
    (<Line2> new).c = c
    (<Line2> new).radius = radius
    return <Line2> new

# Queries and assertions
cdef inline bint issegment(u): return PyObject_TypeCheck(u, Segment2Type)
cdef inline bint isline(u): return PyObject_TypeCheck(u, Line2Type)
cdef inline bint isray(u): return PyObject_TypeCheck(u, Ray2Type)

# Arithmetic and comparison

# Generic
cdef inline check_empty(dict d):
    if not len(d) == 0:
        k, _ = d.popitem()
        raise TypeError(f'invalid argument: {k}')

# Constants
cdef PyTypeObject* Segment2Type = <PyTypeObject*> Segment2
cdef PyTypeObject* Line2Type = <PyTypeObject*> Line2
cdef PyTypeObject* Ray2Type = <PyTypeObject*> Ray2

oo = float('inf')
cdef V.Vec2 origin = V.Vec2(0, 0)
