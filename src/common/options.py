class PeptideOptions:
    def __init__(self, max_sequence_length=50):
        self.max_sequence_length = max_sequence_length

        self.amino_acids = 'ARNDCEQGHILKMFPSTWYV'
        self.modifications = {}


    def amino_acid_size(self):
        return len(self.amino_acids) + sum([
            len(mods) for mods in self.modifications.values()
        ])
    
    