import numpy as np
from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger

from .options import PeptideMS2Options
from .preprocessing import DataConverter
from .modeling import build_model, load_model


def split_train_validate(x, y, validate_percent=.33, seed=None):
    length = len(x)
    np.random.seed(seed)
    indexs = np.random.permutation(length)
    train_end = int((1 - validate_percent) * length)
    train_indexs = indexs[:train_end]
    validate_indexs = indexs[train_end:]
    x_train = x[train_indexs]
    y_train = y[train_indexs]
    x_validate = x[validate_indexs]
    y_validate = y[validate_indexs]
    return x_train, y_train, x_validate, y_validate, train_indexs, validate_indexs


class PeptideMS2Trainer:

    def __init__(self, options=PeptideMS2Options.default(), save_path='epoch_{epoch:03d}.hdf5', save_best_only=True, log_path='training.log'):
        self.options = options
        self.converter = DataConverter(self.options)
        self.model = None
        self.save_path = save_path
        self.save_best_only = save_best_only
        self.log_path = log_path
    

    def load_model(self, model_path, **kwargs):
        self.model = load_model(model_path, **kwargs)
    
    def save_model(self, path):
        self.model.save(path)
    
    
    def train_with_tensor(self, x_train, y_train, x_validate, y_validate, epochs=100, patience=15):
        csvlogger = CSVLogger(self.log_path)
        earlystopper = EarlyStopping(monitor='val_loss', patience=patience, verbose=1)
        if self.save_path is not None:
            checkpointer = ModelCheckpoint(filepath=self.save_path, verbose=1, save_best_only=self.save_best_only)
            callbacks = [checkpointer, csvlogger, earlystopper]
        else:
            callbacks = [csvlogger, earlystopper]

        if self.model is None:
            self.model = build_model(options=self.options)

        history = self.model.fit(
            x_train, y_train, epochs=epochs,
            validation_data=(x_validate, y_validate),
            callbacks=callbacks
        )

        return {
            'history': history.history
        }


    def train(self, data, epochs=100, patience=15, validate_percent=.33, seed=0):
        x, y = self.converter.data_to_tensor(data)
        x_train, y_train, x_validate, y_validate, train_indexs, validate_indexs = split_train_validate(x, y, validate_percent=validate_percent, seed=seed)
        split = {
            'validate_percent': validate_percent,
            'seed': seed,
            'train': train_indexs.tolist(), 
            'validate': validate_indexs.tolist()
        }

        result = self.train_with_tensor(
            x_train, y_train, x_validate, y_validate, 
            epochs=epochs, patience=patience
        )

        result['split'] = split
        return result