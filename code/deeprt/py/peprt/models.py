import keras.backend as K
import numpy as np
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.core import Dense, Dropout, Masking, Flatten
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional
from keras.models import Sequential
from keras.preprocessing.sequence import pad_sequences

aa_size = 20
max_sequence_length = 50


def peptide_to_tensor(sequences):
    def aa_to_vector(aa):
        vec = np.zeros(aa_size, dtype=int)
        vec['ARNDCEQGHILKMFPSTWYV'.index(aa)] = 1
        return vec

    def seq_to_tensor(seq):
        return [aa_to_vector(aa) for aa in seq]

    return pad_sequences(
        [seq_to_tensor(seq) for seq in sequences],
        maxlen=max_sequence_length,
        padding='post')


def normalize(x, min, max):
    return (x - min) / (max - min)

def denormalize(x, min, max):
    return x * (max - min) + min


def build_model():
    model = Sequential()
    model.add(
        Conv1D(
            filters=64,
            kernel_size=5,
            activation="relu",
            input_shape=(max_sequence_length, aa_size)))            
    model.add(MaxPooling1D(pool_size=2, strides=2))
    model.add(Bidirectional(LSTM(128, return_sequences=True)))
    model.add(Dropout(0.5))
    model.add(Flatten())
    model.add(Dense(512, activation='relu'))
    model.add(Dense(256, activation='relu'))
    model.add(Dense(1, activation='relu'))
    model.compile(
        loss="mean_absolute_error",
        optimizer="adam",
        metrics=["mean_absolute_percentage_error"])
    return model


def predict(model, sequences, rt_min=-50, rt_max=150):
    x = peptide_to_tensor(sequences)
    y = model.predict(x)
    return denormalize(y, min=rt_min, max=rt_max)