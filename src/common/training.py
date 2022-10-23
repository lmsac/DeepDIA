from keras.callbacks import ModelCheckpoint, CSVLogger, \
    EarlyStopping, ReduceLROnPlateau
from keras.models import load_model
import keras.backend as K

from .split import split_train_validate


class TrainerBase:
    def __init__(self, options, converter,
                 save_path='epoch_{epoch:03d}.hdf5', 
                 save_best_only=True,
                 log_path='training.log'):
        self.options = options
        self.converter = converter
        self.model = None
        self.save_path = save_path
        self.save_best_only = save_best_only
        self.log_path = log_path


    def use_reduced_lr(self, factor=0.1,
                       patience=5,
                       verbose=1,
                       **kwargs):
        self.reduce_lr = dict(
            factor=factor,
            patience=patience,
            verbose=verbose,
            **kwargs
        )


    def build_model(self):
        raise NotImplementedError()


    def load_model(self, model_path, **kwargs):
        self.model = load_model(model_path, **kwargs)
        return self.model


    def save_model(self, path):
        self.model.save(path)


    def train_with_tensor(self, x_train, y_train, 
                          x_validate, y_validate,
                          epochs=100, patience=15, 
                          lr=None):
        callbacks = []
        
        if self.log_path is not None:
            callbacks.append(CSVLogger(self.log_path))

        if self.save_path is not None:
            callbacks.append(ModelCheckpoint(
                filepath=self.save_path, verbose=1,
                save_best_only=self.save_best_only
            ))

        if patience is not None:
            callbacks.append(EarlyStopping(
                monitor='val_loss',
                patience=patience, verbose=1
            ))
        
        reduce_lr = getattr(self, 'reduce_lr')
        if reduce_lr:
            callbacks.append(ReduceLROnPlateau(
                **reduce_lr
            ))

        if self.model is None:
            self.build_model()

        if lr is not None:
            K.set_value(
               self. model.optimizer.lr, 
               lr
            )

        history = self.model.fit(
            x_train, y_train, epochs=epochs,
            validation_data=(x_validate, y_validate),
            callbacks=callbacks
        )

        return {
            'history': history.history
        }


    def train(self, data, 
              validate_percent=.33, seed=None,
              **kwargs):
        x, y = self.converter.data_to_tensor(data)
        x_train, y_train, x_validate, y_validate, \
            train_indexs, validate_indexs = \
            split_train_validate(
                x, y,
                validate_percent=validate_percent,
                seed=seed
            )
        split = {
            'validate_percent': validate_percent,
            'seed': seed,
            'train': train_indexs.tolist(),
            'validate': validate_indexs.tolist()
        }

        result = self.train_with_tensor(
            x_train, y_train, x_validate, y_validate,
             **kwargs
        )

        result['split'] = split
        return result

