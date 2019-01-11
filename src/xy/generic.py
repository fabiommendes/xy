from generic import generic
from numbers import Number

__all__ = ['simeq', 'assert_simeq']


@generic
def simeq(x, y, tol=1e-6):
    """Test approximate equality between both arguments."""

    try:
        func = x.is_almost_equal
    except AttributeError:
        try:
            func = y.is_almost_equal
        except AttributeError:
            raise TypeError(
                f'cannot check approximate equality of '
                f'{x.__class__.__name__} and {y.__class__.__name__}')
        else:
            return func(x, tol=tol)
    else:
        return func(y, tol=tol)


@simeq.register(Number, Number)
def _(x, y, tol=1e-6):
    return abs(x - y) < tol


@simeq.register(list, list)
@simeq.register(tuple, tuple)
def _(x, y, tol=1e-6):
    eq = simeq
    return len(x) == len(y) and all(eq(xi, yi, tol=tol) for xi, yi in zip(x, y))


def assert_simeq(u, v, tol=1e-6):
    """Assert approximate equality between both arguments."""

    assert simeq(u, v, tol=tol), f'{u} is different from {v}'
