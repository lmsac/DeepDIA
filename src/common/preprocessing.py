import re
import numpy as np

try:
    from keras.utils import pad_sequences
except ImportError:
    from keras.preprocessing.sequence import pad_sequences


class PeptideDataConverter:
    def __init__(self, options):
        self.options = options


    def amino_acid_to_vector(self, amino_acid, modification=None):
        vec = np.zeros(self.options.amino_acid_size(), dtype=int)
        if modification == None:
            idx = self.options.amino_acids.find(amino_acid)
        else:
            idx = 0
            find = False
            for k, v in self.options.modifications.items():
                if k == amino_acid:
                    for i, m in enumerate(v):
                        if m == modification:
                            idx += i
                            find = True
                            break
                    break
                else:
                    idx += len(v)
            if find:
                idx += len(self.options.amino_acids)
            else:
                idx = self.options.amino_acids.find(amino_acid)
        if idx < 0:
            raise ValueError('invalid amino acid: ' + str(amino_acid))
        vec[idx] = 1
        return vec


    def peptide_to_array(self, sequence, modification=None):
        sequence = sequence.upper()
        if modification == None or modification == '' or str(modification) == 'nan':
            return [self.amino_acid_to_vector(aa) for aa in sequence]

        mod_list = [None] * len(sequence)
        for mod in modification.split(';'):
            mat = re.search('[A-Z]([0-9]+)\\((.*)\\)', mod)
            if mat is not None:
                pos = int(mat.group(1)) - 1
                name = mat.group(2)
                mod_list[pos] = name
            #else:
            #    raise ValueError(modification)

        return [
            self.amino_acid_to_vector(aa, mod)
            for aa, mod in zip(sequence, mod_list)
        ]


    def peptides_to_tensor(self, sequences, modifications=None):
        if modifications is None:
            modifications = [None] * len(sequences)
        return pad_sequences(
            [
                self.peptide_to_array(seq, mods)
                for seq, mods in zip(sequences, modifications)
            ],
            maxlen=self.options.max_sequence_length,
            padding='post'
        )

