import abc


class Normed(metaclass=abc.ABCMeta):
    """
    Base class for all objects that have a notion of norm
    """

    @abc.abstractmethod
    def __abs__(self):
        return self.norm()

    @abc.abstractmethod
    def dot(self, other):
        """Dot (inner) product with 'other' object."""

        raise NotImplementedError

    @abc.abstractmethod
    def norm(self):
        """Returns the norm of an object."""

        tname = type(self).__name__
        raise TypeError('%s object does not define a norm' % tname)

    @abc.abstractmethod
    def norm_sqr(self):
        """
        Returns the squared norm of an object using the desired metric.

        If object does not define norm_sqr, it tries to compute the norm using
        obj.norm() and squaring it.
        """

        value = self.norm()
        return value * value

    @abc.abstractmethod
    def normalize(self):
        """Normalizes according to the default norm."""

        try:
            return self / self.norm()
        except ZeroDivisionError:
            raise ValueError('null element cannot be normalize')

    @abc.abstractmethod
    def is_normalized(self, tol=1e-6):
        """
        Return True if the norm equals one within the given tolerance.
        """

        return abs(self.norm() - 1) <= tol

    @abc.abstractmethod
    def is_null(self, tol=1e-6):
        """Return true if object has a null norm."""

        return abs(self.norm() - tol) <= tol
