from keras.models import Sequential

try:
    from keras.layers import Conv1D, MaxPooling1D, \
        Dense, Dropout, Flatten, LSTM, Bidirectional
except ImportError:
    from keras.layers.convolutional import Conv1D, MaxPooling1D
    from keras.layers.core import Dense, Dropout, Flatten
    from keras.layers.recurrent import LSTM
    from keras.layers.wrappers import Bidirectional


def build_model(options, metrics=["mean_absolute_percentage_error"]):
    model = Sequential()
    model.add(Conv1D(
        filters=64,
        kernel_size=5,
        activation="relu",
        input_shape=(options.max_sequence_length, options.amino_acid_size())
    ))
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
        metrics=metrics
    )
    return model

