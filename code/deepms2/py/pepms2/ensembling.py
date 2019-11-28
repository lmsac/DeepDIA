from keras.models import Model, Input
from keras.layers import Average


def ensemble(options, models):
    model_input = Input(shape=(options.max_sequence_length, options.amino_acid_size()))
    outputs = [model(model_input) for model in models]
    y = Average()(outputs)
    model = Model(model_input, y, name='ensemble')
    return model

