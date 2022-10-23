import numpy as np
from common.preprocessing import PeptideDataConverter

try:
    from keras.utils import pad_sequences
except ImportError:
    from keras.preprocessing.sequence import pad_sequences


def normalize_by_max(x):
    x = np.clip(x, a_min=0, a_max=None)
    xmax = np.max(x)
    if xmax > 0:
        return np.divide(x, xmax)
    else:
        return x


class PeptideMS2DataConverter(PeptideDataConverter):
    def __init__(self, options):
        super(PeptideMS2DataConverter, self).__init__(options=options)

        self.transform_input = lambda x: np.log2(x + 1)
        self.transform_output = lambda x: np.power(2, x) - 1
        self.normalize = normalize_by_max


    def ions_to_tensor(self, ions):
        def column_stack_ions(x):
            arrays = [
                (x[frag[0]] if not frag[1] else x[frag[0]][::-1]) \
                if frag[0] in x else [0]
                for frag in self.options.fragments
            ]
            maxlen = max(map(len, arrays))
            arrays = [
                x if len(x) == maxlen else x + [0] * (maxlen - len(x))
                for x in arrays
            ]
            return np.column_stack(arrays)

        return pad_sequences(
            [
                self.transform_input(self.normalize(column_stack_ions(x)))
                for x in ions
            ],
            maxlen=self.options.max_sequence_length - 1,
            padding='post',
            dtype='float'
        )


    def data_to_tensor(self, data):
        sequences = [spec['peptide'] for spec in data]
        modifications = [
            spec['modification'] if 'modification' in spec else None
            for spec in data
        ]
        ions = [spec['ions'] for spec in data]
        x = self.peptides_to_tensor(sequences, modifications)
        y = self.ions_to_tensor(ions)
        return x, y


    def tensor_to_ions(self, tensor, sequence_lengths):
        def decompose_ions(y):
            return {
                frag[0]: y[:, i].tolist() \
                if not frag[1] else y[:, i].tolist()[::-1]
                for i, frag in enumerate(self.options.fragments)
            }

        return [
            decompose_ions(
                self.transform_output(self.normalize(y)) \
                    [:(sequence_lengths[i] - 1)]
            )
            for i, y in enumerate(tensor)
        ]

