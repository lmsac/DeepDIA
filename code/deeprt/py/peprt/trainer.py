from . import models

import numpy as np
import pandas as pd
from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from keras.models import load_model

def data_to_tensors(data, rt_min=-50, rt_max=150):
    peptides = data[["sequence"]].values.flatten()
    rt = data[["irt"]].values    
    x = models.peptide_to_tensor(peptides)
    y = models.normalize(rt, min=rt_min, max=rt_max)
    return x, y


def split_train_validate(x, y, validate_percent=.33, seed=None):
    length = len(x)
    np.random.seed(seed)
    indexs = np.random.permutation(length)    
    train_end = int((1 - validate_percent) * length)
    train_indexs = indexs[:train_end]
    validate_indexs = indexs[train_end:]
    x_train = x[train_indexs]
    y_train = y[train_indexs]
    x_validate =x[validate_indexs]
    y_validate =y[validate_indexs]
    return x_train, y_train, x_validate, y_validate, train_indexs, validate_indexs


class PeptideRTTrainer:

    def __init__(self, model_path=None, model=None, rt_min=-50, rt_max=150, save_path="bestmodel.hdf5", save_best_only=True, log_path='training.log'):
        if model_path is not None:
            model = load_model(model_path)
        elif model is None:
            model = models.build_model()        
        self.model = model
        self.rt_min = rt_min
        self.rt_max = rt_max
        self.save_path = save_path
        self.save_best_only = save_best_only
        self.log_path = log_path
    
    def get_model(self):
        return self.model

    def save_model(self, path):
        self.model.save(path)
    
    def train(self, data, epochs=100, patience=15, validate_percent=.33, seed=0):
        x, y = data_to_tensors(data, rt_min=self.rt_min, rt_max=self.rt_max)
        x_train, y_train, x_validate, y_validate, train_indexs, validate_indexs = split_train_validate(x, y, validate_percent=0.33, seed=0)
        split = {
            'validate_percent': validate_percent,
            'seed': seed,
            'train': train_indexs.tolist(), 
            'validate': validate_indexs.tolist()                
        }

        csvlogger = CSVLogger(self.log_path)
        earlystopper = EarlyStopping(monitor='val_loss', patience=patience, verbose=1)
        if self.save_path is not None:
            checkpointer = ModelCheckpoint(filepath=self.save_path, verbose=1, save_best_only=self.save_best_only)
            callbacks = [checkpointer, csvlogger, earlystopper]
        else:
            callbacks = [csvlogger, earlystopper]                    

        history = self.model.fit(
            x_train, y_train, epochs=epochs,
            validation_data=(x_validate, y_validate),
            callbacks=callbacks)

        return {
            'split': split,
            'history': history.history
        }

