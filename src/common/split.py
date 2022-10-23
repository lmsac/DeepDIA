import numpy as np


def split_train_validate(x, y, validate_percent=.33, seed=None):
    if seed is None:
        random = np.random
    else:
        random = np.random.RandomState(seed=seed)

    length = len(x)
    indexs = random.permutation(length)
    train_end = int((1 - validate_percent) * length)
    train_indexs = indexs[:train_end]
    validate_indexs = indexs[train_end:]
    x_train = x[train_indexs]
    y_train = y[train_indexs]
    x_validate = x[validate_indexs]
    y_validate = y[validate_indexs]
    return x_train, y_train, \
        x_validate, y_validate, \
        train_indexs, validate_indexs

