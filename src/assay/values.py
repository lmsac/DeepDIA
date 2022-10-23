import numpy as np
from collections import OrderedDict

from util import get_value_from_dict, set_value_to_dict


def get_assay_values(assay, params, func_args=None):
    d = OrderedDict()

    for p in params:
        value = None

        if not isinstance(p, dict):
            d[p] = get_value_from_dict(assay, p)
            continue

        name = p.get('name')
        
        func = p.get('function', None)
        if callable(func):
            if func_args is None:
                value = func(assay)
            elif isinstance(func_args, dict):
                value = func(assay, **func_args)
            elif isinstance(func_args, list):
                value = func(assay, *func_args)
            else:
                value = func(assay, func_args)
            d[name] = value
            continue
        
        path = p.get('path', name)
        if path is not None:
            value = get_value_from_dict(assay, path)
            if value is not None:
                convert = p.get('convert', None)
                if callable(convert):
                    value = convert(value)            
            else:
                value = p.get('default', None)
            d[name] = value
        
    return d


def set_assay_values(assay, params, func_args=None, **values):
    for p in params:
        if not isinstance(p, dict):
            if p in values:
                set_value_to_dict(assay, p, values[p])
            continue

        name = p.get('name')
        if name in values:
            value = values.get(name)
        else:
            if 'default' in p:
                value = p.get('default')
            else:
                continue        

        func = p.get('function', None)
        if callable(func):
            func(assay, value, **func_args)
            
        path = p.get('path', name)
        if path is not None:
            convert = p.get('convert', None)
            if callable(convert):
                if p.get('list', False):
                    value = list(map(convert, value))
                else:
                    value = convert(value)

            set_value_to_dict(assay, path, value)
            continue

        

    return assay


def assign_assays_values(assays, keys, params, data):
    key_columns = []
    key_params = []
    for k in keys:
        name = k.get('name') if isinstance(k, dict) else k
        if name in data.columns:
            key_columns.append(name)
            key_params.append(k)

    data = data.set_index(key_columns)

    for assay in assays:
        assay_id = tuple(get_assay_values(assay, key_params).values())

        values = data.loc[assay_id].to_dict()
        for name, value in list(values.items()):
            if isinstance(value, np.generic):
                 values[name] = value.item()

        set_assay_values(assay, params, **values)

    return assays


