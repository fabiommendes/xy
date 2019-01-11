from math import pi, sqrt

import pytest

from xy import simeq


class VectorInterface:
    """
    Generic tests for the vector API
    """

    base_cls = None
    size = property(lambda self: self.base_cls.size)

    def test_equality_against_tuples(self, u):
        assert tuple(u) == u
        assert u == tuple(u)

    def test_vec_has_no_dict(self, u):
        with pytest.raises(AttributeError):
            print('dict:', u.__dict__)

    def test_from_flat_data(self):
        args = range(self.size)
        assert self.base_cls(*args) == self.base_cls.from_flat(args)

    def test_clamp_to_value(self, unity):
        assert simeq(unity.with_length(2), 2 * unity)
        assert simeq(unity.with_length(0.5), 0.5 * unity)

    def test_clamp_interval(self, unity):
        assert unity.clamp(0.5, 2) == unity

    def test_clamp_missing_interval(self, unity):
        assert simeq(unity.clamp(2, 3), 2 * unity)
        assert simeq(unity.clamp(0.1, 0.5), 0.5 * unity)

    def test_lerp(self, u, v):
        assert simeq(u.lerp(v), v.lerp(u))
        assert simeq(u.midpoint(v), u.lerp(v))
        assert simeq(u.lerp(v, 0), v)
        assert simeq(u.lerp(v, 1), u)

    def test_middle(self, unity, null):
        assert simeq(unity.midpoint(null), null.midpoint(unity))
        assert simeq(unity.midpoint(null), unity / 2)

    def test_distance(self, unity, null):
        assert simeq(unity.distance_to(unity), 0)
        assert simeq(unity.distance_to(null), 1)
        assert simeq(unity.distance_to(-unity), 2)

    def test_angle(self, unity):
        assert simeq(unity.angle_to(unity), 0)
        assert simeq(unity.angle_to(-unity), 180)
        assert simeq(unity.angle_to_rad(-unity), pi)

    def test_vector_norm_defaults_to_euclidean(self):
        vec = self.base_cls(*(1 for _ in range(self.size)))
        assert simeq(vec.norm(), sqrt(self.size))
        assert simeq(abs(vec), sqrt(self.size))

    def test_l1_norm(self):
        u = self.base_cls(*(1 for _ in range(self.size)))
        assert abs(u.norm_l1() - self.size) < 1e-6

    def test_floordiv(self):
        u = self.base_cls(*(3 for _ in range(self.size))) // 2
        assert list(u) == [1] * self.size, u


class VectorInvalidOperations:
    """
    Tests if invalid vector operations raise the correct errors
    """

    def test_invalid_scalar_operations(self, u):
        with pytest.raises(TypeError):
            print(u + 1)
        with pytest.raises(TypeError):
            print(u - 1)

    def test_invalid_mul_tuple(self, u):
        with pytest.raises(TypeError):
            print(u * (1, 2))

    def test_invalid_mul_vec(self, u):
        with pytest.raises(TypeError):
            print(u * u)

    def test_invalid_div_tuple(self, u):
        with pytest.raises(TypeError):
            print(u / (1, 2))
        with pytest.raises(TypeError):
            # noinspection PyUnresolvedReferences
            print((1, 2) / u)

    def test_invalid_div_vec(self, u):
        with pytest.raises(TypeError):
            print(u / u)

    def test_invalid_div_scalar(self, u):
        with pytest.raises(TypeError):
            print(1 / u)

    def test_vec_almost_equal(self, u, v):
        v = u + v / 1e9
        w = u + 0 * v
        assert u.is_almost_equal(v)
        assert v.is_almost_equal(u)
        assert u.is_almost_equal(w)
        assert w.is_almost_equal(u)