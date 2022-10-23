import os
import re
import itertools
import pickle
import json
import numpy as np

def list_files(path='.', pattern=None, recursive=False, include_dirs=False):
    if recursive:
        return list(itertools.chain.from_iterable(
            [
                os.path.join(t[0], x) 
                for x in (t[1] + t[2] if include_dirs else t[2])
                if pattern is None or \
                re.search(pattern, x) is not None
            ]
            for t in os.walk(path)
        ))
    else:
        return [
             t
             for t in os.listdir(path)
             if pattern is None or \
             re.search(pattern, t) is not None
        ]


# https://github.com/mpld3/mpld3/issues/434#issuecomment-340255689
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
            np.int16, np.int32, np.int64, np.uint8,
            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
            np.float64)):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)		


def save_json(data, file, **kwargs):
    with open(file, 'w') as f:
        json.dump(data, f, cls=NumpyEncoder, **kwargs)        

def load_json(file, **kwargs):
    with open(file, 'r') as f:
        return json.load(f, **kwargs)
    

def save_pickle(data, file, **kwargs):
    with open(file, 'wb') as f:
        pickle.dump(data, f, **kwargs)
        
def load_pickle(file, **kwargs):
    with open(file, 'rb') as f:
        return pickle.load(f, **kwargs)
    
    