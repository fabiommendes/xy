from math import sqrt

import pytest


class LinearFixtures:
    """
    Automatic fixtures for items items that can be created from a list of
    scalar data.
    """
    base_cls = None

    @property
    def size(self):
        return self.base_cls.size

    @pytest.fixture
    def unity(self):
        return self.base_cls(*(1 / sqrt(self.size) for _ in range(self.size)))

    @pytest.fixture
    def null(self):
        return self.base_cls(*(0.0 for _ in range(self.size)))

    @pytest.fixture
    def u(self):
        return self.base_cls(*(x + 1 for x in range(self.size )))

    @pytest.fixture
    def v(self):
        return self.base_cls(*(x + 1 for x in reversed(range(self.size))))

    @pytest.fixture
    def e1(self):
        args = [0 for _ in range(self.size)]
        args[0] = 1
        return self.base_cls(*args)

    @pytest.fixture
    def e2(self):
        args = [0 for _ in range(self.size)]
        args[1] = 1
        return self.base_cls(*args)

    @pytest.fixture
    def e3(self):
        args = [0 for _ in range(self.size)]
        args[2] = 1
        return self.base_cls(*args)

    @pytest.fixture
    def e4(self):
        args = [0 for _ in range(self.size)]
        args[3] = 1
        return self.base_cls(*args)


class SequenceInterface:
    """
    For objects that behave like sequences.
    """
    base_cls = None

    def test_has_basic_sequence_interface(self, u):
        assert len(u) == self.base_cls.size

    def test_equality_with_non_tuple_sequences(self, u):
        assert u != '12'
        assert u is not None
        assert u != list(u) + [0]

    def test_can_get_item(self, u):
        assert u[0] is not None
        assert u[-1] is not None

    def test_can_slice(self, u):
        assert u[:] == list(u)

    def test_index_error_for_invalid_indexes(self, u):
        n = len(u)
        with pytest.raises(IndexError):
            print('element:', u[n])
        with pytest.raises(IndexError):
            print('element:', u[-(n + 1)])

    def test_copy_is_equal_to_itself(self, u):
        assert u == u.copy()