class Style(dict):
    """
    A mapping between style properties and values.
    """

    def __getattr__(self, attr):
        if attr.startswith('_'):
            raise AttributeError(attr)
        try:
            return self[attr]
        except KeyError:
            raise AttributeError

    def to_json(self):
        return self

    def to_svg_attrs(self):
        return {k: str(v) for k, v in self.items()}


def as_style(st):
    return st if isinstance(st, Style) else Style(st)
