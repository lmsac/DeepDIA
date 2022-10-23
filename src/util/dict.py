def get_value_from_dict(d, path, default=None):
    if not isinstance(path, list):
        path = [path]

    value = None
    for name in path:
        if isinstance(d, dict):
            value = d.get(name, default)
            d = value
        else:
            return None
    return value


def set_value_to_dict(d, path, value):
    if not isinstance(path, list):
        path = [path]

    dd = d
    for i, name in enumerate(path):
        if i == len(path) - 1:
            dd[name] = value
        else:
            x = dd.get(name, None)
            if not isinstance(x, dict):
                x = {}
                dd[name] = x
            dd = x
    return d

