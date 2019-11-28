import re
import numpy as np
from keras.preprocessing.sequence import pad_sequences


def normalize_by_max(x):
    x = np.clip(x, a_min=0, a_max=None)
    xmax = np.max(x)
    if xmax > 0:
        return np.divide(x, xmax)
    else:
        return x


class DataConverter:

    def __init__(self, options):
        self.options = options
        self.transform_input = lambda x: np.log2(x + 1)
        self.transform_output = lambda x: np.power(2, x) - 1
        self.normalize = normalize_by_max

    def amino_acid_to_vector(self, amino_acid, modification=None):
        vec = np.zeros(self.options.amino_acid_size(), dtype=int)
        if modification == None:
            idx = self.options.amino_acids.index(amino_acid)
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
                idx = self.options.amino_acids.index(amino_acid)
        vec[idx] = 1
        return vec

    
    def peptide_to_array(self, sequence, modification=None):
        sequence = sequence.upper()
        if modification == None or modification == '':
            return [self.amino_acid_to_vector(aa) for aa in sequence]
        
        mod_list = [None] * len(sequence)
        for mod in modification.split(';'):
            pos = int(re.search('[A-Z]([0-9]+)\\(', mod).group(1)) - 1
            name = re.search('\\((.*)\\)', mod).group(1)
            mod_list[pos] = name
        
        return [self.amino_acid_to_vector(aa, mod) for aa, mod in zip(sequence, mod_list)]


    def peptides_to_tensor(self, sequences, modifications=None):
        if modifications is None:
            modifications = [None] * len(sequences)
        return pad_sequences(
            [self.peptide_to_array(seq, mods) for seq, mods in zip(sequences, modifications)],
            maxlen=self.options.max_sequence_length,
            padding='post'
        )


    def ions_to_tensor(self, ions):
        def column_stack_ions(x):   
            return np.column_stack(tuple(
                x[frag[0]] if not frag[1] else x[frag[0]][::-1]
                for frag in self.options.fragments
            ))

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
        modifications = [spec['modification'] if 'modification' in spec else None for spec in data]
        ions = [spec['ions'] for spec in data]
        x = self.peptides_to_tensor(sequences, modifications)
        y = self.ions_to_tensor(ions)
        return x, y


    def tensor_to_ions(self, tensor, sequence_lengths):
        def decompose_ions(y):                
            return {
                frag[0]: y[:, i].tolist() if not frag[1] else y[:, i].tolist()[::-1]
                for i, frag in enumerate(self.options.fragments)
            }

        return [
            decompose_ions(self.transform_output(self.normalize(y))[:(sequence_lengths[i] - 1)])
            for i, y in enumerate(tensor)
        ]