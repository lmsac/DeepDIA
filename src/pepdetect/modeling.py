from keras.models import Sequential
from keras.regularizers import l2

try:
    from keras.layers import Conv1D, MaxPooling1D, \
        Dense, Dropout, Flatten, LSTM, Bidirectional
except ImportError:
    from keras.layers.convolutional import Conv1D, MaxPooling1D
    from keras.layers.core import Dense, Dropout, Flatten
    from keras.layers.recurrent import LSTM
    from keras.layers.wrappers import Bidirectional


def build_model(options, metrics=["mean_absolute_error"]):
    model = Sequential()
    model.add(Conv1D(
        filters=64,
        kernel_size=5,
        activation="relu",
        input_shape=(options.max_sequence_length, options.amino_acid_size()),
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
        metrics=metrics
    )
    return model

