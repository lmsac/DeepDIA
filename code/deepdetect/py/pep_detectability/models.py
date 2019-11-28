import keras.backend as K
import numpy as np
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.core import Dense, Dropout, Masking, Flatten
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional
from keras.models import Sequential
from keras.regularizers import l2
from keras.preprocessing.sequence import pad_sequences
import pandas as pd

aa_size = 20 + 2
max_aa_length = 50
terminal_length = 7
max_sequence_length = max_aa_length + terminal_length * 2

def get_peptide_sequence(data):
    return data['nTerminal'].str.slice(start=-terminal_length) \
        .str.cat(data['sequence'], sep='.') \
        .str.cat(data['cTerminal'].str.slice(stop=terminal_length), sep='.') \
        .values

def peptide_to_tensor(sequences):
    def aa_to_vector(aa):
        vec = np.zeros(aa_size, dtype=int)
        vec['_.ARNDCEQGHILKMFPSTWYV'.index(aa)] = 1
        return vec

    def seq_to_tensor(seq):
        return [aa_to_vector(aa) for aa in seq]

    return pad_sequences(
        [seq_to_tensor(seq) for seq in sequences],
        maxlen=max_sequence_length,
        padding='post')

def build_model():
    model = Sequential()
    model.add(Conv1D(
        filters=64,
        kernel_size=5,
        activation="relu",
        input_shape=(max_sequence_length, aa_size),
        kernel_regularizer=l2(l=1e-4),
        bias_regularizer=l2(l=1e-4)
    ))            
    model.add(MaxPooling1D(pool_size=2, strides=2))
    model.add(Bidirectional(LSTM(
        128, return_sequences=True,
        activation="tanh",
        recurrent_regularizer=l2(l=1e-5), 
        kernel_regularizer=l2(l=1e-5),
        bias_regularizer=l2(l=1e-5)
    )))
    model.add(Dropout(0.5))
    model.add(Flatten())
    model.add(Dense(
        64, activation='relu', 
        kernel_regularizer=l2(l=0.001),
        bias_regularizer=l2(l=0.001)
    ))
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='relu'))
    model.compile(
        loss="mean_squared_error",
        optimizer="adam",
        metrics=["mean_absolute_error"])
    return model

