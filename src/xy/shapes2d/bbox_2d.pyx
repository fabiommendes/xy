# distutils: language = c++
# cython: cdivision = True

include "../include/text.pxi"
cimport cython
cimport xy.linalg2d.vector_2d as V
from libc.math cimport sqrt
from cpython cimport PyObject_TypeCheck


# ------------------------------------------------------------------------------
# RECTANGULAR BOUNDING BOX 2 by 2
# ------------------------------------------------------------------------------

@cython.freelist(8)
cdef class BBox2:
    """
    A rectangular axis aligned bounding box.

    Attributes:
        xmin, xmax, ymin, ymax (float):
            AABB limits in x and y directions.
        vertices:
            Sequence of vertices that form the AABB. Starts with the 1st
            quadrant and runs counter clockwise.

    Example:
        Create an BBox by specifying its limits in the x and y directions

        >>> a = BBox2(0, 50, 0, 100)
        >>> a.pos  # center point
        Vec2(25, 50)
    """

    @property
    def pos(self):
        return V.flatvec((self.xmin + self.xmax) / 2, (self.ymin + self.ymax) / 2)

    @property
    def bbox(self):
        return self

    @property
    def vertices(self):
        return (
            V.flatvec(self.xmin, self.ymin), V.flatvec(self.xmax, self.ymin),
            V.flatvec(self.xmax, self.ymax), V.flatvec(self.xmin, self.ymax)
        )

    @property
    def vertice_ne(self):
        return V.flatvec(self.xmax, self.ymax)

    @property
    def vertice_nw(self):
        return V.flatvec(self.xmin, self.ymax)

    @property
    def vertice_se(self):
        return V.flatvec(self.xmax, self.ymin)

    @property
    def vertice_sw(self):
        return V.flatvec(self.xmin, self.ymin)

    @property
    def outer_radius(self):
        dx = self.xmax - self.xmin
        dy = self.ymax - self.ymin
        return sqrt(dx * dx + dy * dy) / 2

    @property
    def shadow_x(self):
        return self.xmin, self.xmax

    @property
    def shadow_y(self):
        return self.ymin, self.ymax

    @property
    def width_x(self):
        return self.xmax - self.xmin

    @property
    def width_y(self):
        return self.ymax - self.ymin

    @property
    def coords(self):
        yield self.xmin
        yield self.xmax
        yield self.ymin
        yield self.ymax

    #
    # Constructors
    #
    @classmethod
    def from_coords(cls, xmin, xmax, ymin, ymax):
        """
        Creates a new AABB from xmin, xmax, ymin, ymax coordinates.

        This constructor do not accept alternative signatures and is slightly
        faster than the default one.
        """
        if xmin > xmax:
            raise ValueError('xmax < xmin')
        if ymin > ymax:
            raise ValueError('ymax < ymin')
        return bbox2(xmin, xmax, ymin, ymax)

    @classmethod
    def from_points(self, points):
        """Construct the smallest bounding box that can fit the given sequence
        of points"""
        cdef double xmin, xmax, ymin, ymax, x, y
        it = iter(points)
        x, y = next(it)
        xmin = xmax = x
        ymin = ymax = y

        for pt in it:
            try:
                V.set_xy(pt, &x, &y)
            except TypeError:
                x, y = pt
            if x < xmin: xmin = x
            if x > xmax: xmax = x
            if y < ymin: ymin = y
            if y > ymax: ymax = y

        return bbox2(xmin, xmax, ymin, ymax)

    def __cinit__(self, *args, **kwargs):
        cdef double xmin, xmax, ymin, ymax, x, y, width = 0, height = 0

        if args and kwargs:
            raise TypeError('cannot set positional and keyword arguments simultaneously')
        elif len(args) == 1:
            xmin, xmax, ymin, ymax = args[0]
        elif args:
            xmin, xmax, ymin, ymax = args
        elif 'pos' in kwargs:
            x, y = kwargs.pop('pos')
            if 'shape' in kwargs:
                width, height = kwargs.pop('shape')
            elif 'width' in kwargs and 'height' in kwargs:
                width = kwargs.pop('width')
                height = kwargs.pop('height')
            xmin = x - width / 2
            xmax = x + width / 2
            ymin = y + height / 2
            ymax = y + height / 2
        elif 'shape' in kwargs:
            width, height = kwargs.pop('shape')
            xmin, xmax = -width / 2, width / 2
            ymin, ymax = -height / 2, height / 2
        elif 'width' in kwargs and 'height' in kwargs:
            width = kwargs.pop('width')
            height = kwargs.pop('height')
            xmin, xmax = -width / 2, width / 2
            ymin, ymax = -height / 2, height / 2
        else:
            raise TypeError('invalid arguments')

        # Final consistency checks
        check_empty(kwargs)
        if xmin > xmax:
            raise ValueError('xmax < xmin')
        if ymin > ymax:
            raise ValueError('ymax < ymin')

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    #
    # Magic methods
    #
    def __repr__(self):
        coords = tuple(fmt(x) for x in self.coords)
        return f'{type(self).__name__}(%g, %g, %g, %g)' % coords

    def __eq__(self, other):
        return eq(<BBox2> self, <BBox2> other) if isbbox(other) else NotImplemented

    #
    # Geometric properties
    #
    def area(self):
        return (self.xmax - self.xmin) * (self.ymax - self.ymin)

    def gyration_radius(self):
        cdef double a = self.xmax - self.xmin, b = self.ymax - self.ymin
        return sqrt((a * a + b * b) / 12)

    def gyration_sqr(self):
        cdef double a = self.xmax - self.xmin, b = self.ymax - self.ymin
        return (a * a + b * b) / 12

    #
    # Geometric transformations
    #
    def translate(self, vec):
        cdef double x, y
        V.set_xy(vec, &x, &y)
        return bbox2(self.xmin + x, self.xmax + x, self.ymin + y, self.ymax + y)

    def scale(self, double scale, point=None):
        return bbox2(self.xmin * scale, self.xmax * scale,
                     self.ymin * scale, self.ymax * scale)

    def scale_around(self, double scale, point):
        cdef double x, y, u
        V.set_xy(point, &x, &y)
        u = 1 - scale
        return bbox2(self.xmin * scale + u * x, self.xmax * scale + u * x,
                     self.ymin * scale + u * y, self.ymax * scale + u * y)

    def flip_x(self): return self
    def flip_y(self): return self
    def flip_xy(self):
        """Flip width and height of bounding box."""
        return bbox2(self.ymin, self.ymax, self.xmin, self.xmax)

    #
    # Queries
    #
    def contain_xy(self, double x, double y):
        return self.xmin <= x <= self.xmax and self.ymin <= y <= self.ymax

    def contain_point(self, point):
        cdef double x, y
        V.set_xy(point, &x, &y)
        return self.contains_xy(x, y)

    def contain_bbox(self, BBox2 other):
        return self.xmin <= other.xmin and self.ymin <= other.ymin \
               and self.xmax >= other.xmax and self.ymax >= other.ymax

    def contain_circle(self, other):
        # return x <= other.xmin <= x + xm and y <= other.ymin <= y + ym
        raise NotImplementedError


# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

# Constructors
cdef inline BBox2 bbox2(double a, double b, double c, double d):
    cdef PyObject* new = _PyObject_New(BBox2Type)
    (<BBox2> new).xmin = a
    (<BBox2> new).xmax = b
    (<BBox2> new).ymin = c
    (<BBox2> new).ymax = d
    return <BBox2> new

# Queries and assertions
cdef inline bint isbbox(u): return PyObject_TypeCheck(u, BBox2Type )

# Arithmetic and comparison
cdef inline bint eq(BBox2 a, BBox2 b): return a.xmin == b.xmin and a.xmax == b.xmax and a.ymin == b.ymin and a.ymax == b.ymax

# Generic
cdef inline check_empty(dict d):
    if not len(d) == 0:
        k, _ = d.popitem()
        raise TypeError(f'invalid argument: {k}')

# Constants
cdef BBox2 unity = bbox2(0, 1, 0, 1)
cdef PyTypeObject* BBox2Type = <PyTypeObject*> BBox2
