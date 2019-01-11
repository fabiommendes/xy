from . import linalg2d


class Property:
    def __init__(self, getter, slot=None, setter=None):
        self.slot = slot
        self.getter = getter
        self.setter = setter

    def __get__(self, obj, tt):
        return self.getter(obj, tt)

    def __set__(self, obj, value):
        return self.setter(obj, value)

    def __set_name__(self, name):
        pass

    def setter(self, setter):
        return type(self)(self.getter, setter=setter, slot=self.slot)


class TypedProperty(Property):
    dtype = None
    coerce = None

    def __set__(self, obj, value):
        if not isinstance(value, self.dtype):
            value = self.coerce(value)
        super().__set__(obj, value)


#
# Create specialized properties
#
def typed_property(dtype, transform=None):
    transform = transform or dtype

    def decorator(func):
        prop = TypedProperty(func)
        prop.dtype = dtype
        prop.transform = transform
        return prop

    return decorator


vec2_property = typed_property(linalg2d.Vec2)
direction2_property = typed_property(linalg2d.Direction2)
mat2_property = typed_property(linalg2d.Mat2)
