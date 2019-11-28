import keras.backend as K
from keras.layers.convolutional import Conv1D
from keras.layers.core import Dense, Dropout, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model as keras_load_model
from keras.models import Sequential

from .options import PeptideMS2Options


def cosine_similarity(y_true, y_pred):        
    length = K.int_shape(y_pred)[1]
    y_true = K.batch_flatten(y_true)
    y_pred = K.batch_flatten(y_pred)
    y_true = K.l2_normalize(y_true, axis=-1)
    y_pred = K.l2_normalize(y_pred, axis=-1)
    cos = K.sum(y_true * y_pred, axis=-1, keepdims=True)
    result = -K.repeat_elements(cos, rep=length, axis=1)
    return result


def build_model(options, metrics=[cosine_similarity]):
    model = Sequential()
    model.add(
        Conv1D(
            filters=64,
            kernel_size=2,
            activation="relu",
            input_shape=(options.max_sequence_length, options.amino_acid_size())))
    model.add(Masking(mask_value=0.))
    model.add(Bidirectional(LSTM(128, return_sequences=True)))
    model.add(Dropout(0.5))
    model.add(TimeDistributed(Dense(options.intensity_size(), activation='relu')))
    model.compile(
        loss="mean_squared_error",
        optimizer="adam",
        metrics=metrics)
    return model

def load_model(file, custom_metrics={'cosine_similarity': cosine_similarity}):
    model = keras_load_model(file, custom_objects=custom_metrics)
    return model

def build_model_from_weights(options, weights_path, metrics=[cosine_similarity]):
    model = build_model(options=options, metrics=metrics)
    model.load_weights(weights_path)
    return model