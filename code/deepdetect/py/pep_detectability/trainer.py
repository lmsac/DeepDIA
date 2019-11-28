from . import models

import numpy as np
import pandas as pd
from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from keras.models import load_model

def data_to_tensors(data):    
    sequences = models.get_peptide_sequence(data)
    x = models.peptide_to_tensor(sequences)
    y = data[["detectability"]].values  
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


class PeptideDetectabilityTrainer:

    def __init__(self, model_path=None, model=None, save_path="bestmodel.hdf5", save_best_only=True, log_path='training.log'):
        if model_path is not None:
            model = load_model(model_path)
        elif model is None:
            model = models.build_model()        
        self.model = model
        self.save_path = save_path
        self.save_best_only = save_best_only
        self.log_path = log_path
    
    def get_model(self):
        return self.model

    def save_model(self, path):
        self.model.save(path)
        
    
    def train(self, data1, data2, epochs=100, patience=15, validate_percent=.33, seed=0):
        x1, y1 = data_to_tensors(data1)
        x1_train, y1_train, x1_validate, y1_validate, \
            train_indexs1, validate_indexs1 = \
            split_train_validate(x1, y1, validate_percent=0.33, seed=0)

        x2, y2 = data_to_tensors(data2)
        x2_train, y2_train, x2_validate, y2_validate, \
            train_indexs2, validate_indexs2 = \
            split_train_validate(x2, y2, validate_percent=0.33, seed=0)
        split = {
            'validate_percent': validate_percent,
            'seed': seed,
            'train': {                
                'positive': train_indexs1.tolist(), 
                'negative': train_indexs2.tolist()
            },
            'validate': {                
                'positive': validate_indexs1.tolist(), 
                'negative': validate_indexs2.tolist()
            }            
        }
        x_train = np.concatenate((x1_train, x2_train))
        y_train = np.concatenate((y1_train, y2_train))
        x_validate = np.concatenate((x1_validate, x2_validate))
        y_validate = np.concatenate((y1_validate, y2_validate))

        result = self.train_with_tensor(
            x_train, y_train, x_validate, y_validate, 
            epochs=epochs, patience=patience
        )
        result['split'] = split
        return result


    def train_with_tensor(self, x_train, y_train, x_validate, y_validate, epochs=100, patience=15):
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
            'history': history.history
        }

