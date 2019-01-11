from functools import singledispatch
from ..linalg2d import Vec2
from ..shapes2d import Segment2, Line2, Segment2, BBox2
from .style import as_style
import hyperpython as hp


class Shape:
    """
    A shape
    """
    def __init__(self, shape, style=None, transform=None):
        self.shape = shape
        self.affine = transform
        self.style = as_style(style)

    def move(self, vec):
        self.shape = self.shape.translate(vec)

    def rotate(self, vec):
        self.shape = self.shape.translate(vec)

    def to_json(self):
        return {'shape': self.shape.to_json(), 'style': self.to_json()}

    def to_svg(self):
        pass

    def to_svg_string(self, embed=True):
        pass


def circle(radius, pos=Vec2(0, 0), **kwargs):
    return Shape(...)


@singledispatch
def svg(x):
    raise TypeError(f'invalid argument type: {type(x).__name__}')


@hp_args.register(Segment2)
def _(x):
    raise NotImplementedError


@hp_args.register(Circle2)
def _(pt):
    return 'circle', {'cx': pt.x, 'cy': pt.y, 'r': pt.radius}
