from math import pi, sqrt

from tests import interfaces as i
from xy import Vec2, simeq


class TestVector2DInterfaces(i.LinearFixtures,
                             i.SequenceInterface,
                             i.VectorInterface,
                             i.VectorInvalidOperations,
                             i.ElementwiseAddition,
                             i.ScalarMultiplication,
                             i.Normed):
    base_cls = Vec2
    norm_list = ['l1', 'linf']

    def test_vector_repr(self):
        assert repr(Vec2(1, 2)) == 'Vec2(1, 2)'
        assert repr(Vec2(0.5, 0.5)) == 'Vec2(0.5, 0.5)'

    def test_polar_coordinates(self):
        assert Vec2.from_polar(1, 0) == Vec2(1, 0)

    def test_rotations(self):
        v = Vec2(1, 0)
        assert v.rotate(90) == Vec2(0, 1)
        assert v.rotate_around(90, Vec2(1, 0)) == v
        assert simeq(v.rotate_rad(pi / 2), Vec2(0, 1))
        assert simeq(v.rotate_around_rad(pi / 2, Vec2(1, 0)), v)

    def test_cross(self):
        assert Vec2(1, 0).cross(Vec2(0, 1)) == 1

    def test_polar(self):
        r, t = Vec2(1, 1).polar()
        assert simeq(r, sqrt(2))
        assert simeq(t, 45)

        r, t = Vec2(1, 1).polar_rad()
        assert simeq(r, sqrt(2))
        assert simeq(t, pi / 4)

    def test_perp(self):
        assert Vec2(1, 0).perpendicular() == Vec2(0, 1)
        assert Vec2(1, 0).perpendicular(ccw=False) == Vec2(0, -1)

    def test_rotated_is_new(self, u):
        assert u.rotate(0.0) is not u
        assert u.rotate(0.0) == u

    def test_rotated_keeps_norm(self, u):
        for t in range(5):
            Z1 = u.norm()
            Z2 = u.rotate(360 * t / 6).norm()
            assert simeq(Z1, Z2)

    def test_triangular_identity_2D(self):
        self.assert_triangular_identity(Vec2(1, 2), Vec2(3, 4))
        self.assert_triangular_identity(Vec2(1, 1), Vec2(1, 1))
        self.assert_triangular_identity(Vec2(1, 2), Vec2(0, 0))


class TestVec2Examples:
    def test_class_parameters(self):
        vec2 = Vec2
        assert vec2.shape == (2,)
        assert vec2.dtype == float
        assert vec2.__name__ == 'Vec2'

    def test_correct_type_promotion_on_vec_creation(self):
        assert isinstance(Vec2(1.0, 2.0), Vec2)
        assert isinstance(Vec2(1, 2.0), Vec2)
        assert isinstance(Vec2(1.0, 2), Vec2)

    def test_vec_equality(self):
        assert Vec2(1, 2) == Vec2(1, 2)
        assert Vec2(1, 2) == Vec2(1.0, 2.0)

    def test_vec_equality_with_tuples(self):
        assert Vec2(1, 2) == (1, 2)
        assert Vec2(1, 2) == (1.0, 2.0)

    def test_reverse_vec_equality_with_tuples(self):
        assert (1, 2) == Vec2(1, 2)
        assert (1.0, 2.0) == Vec2(1, 2)

    def test_tuple_promotion_on_arithmetic_operations(self):
        u = Vec2(1, 2)
        v = (0.0, 0.0)
        assert u + v == u
        assert u - v == u
        assert isinstance(u + v, Vec2)
        assert isinstance(u - v, Vec2)


class TestRegressions:
    def test_clamp(sef):
        u = Vec2(3, 4)
        assert u.clamp(1, 10) == u
        assert u.clamp(2, 4) == u.normalize() * 4
        assert u.with_length(10) == 2 * u
