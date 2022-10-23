from common.options import PeptideOptions


class PeptideDetectabilityOptions(PeptideOptions):
    def __init__(self, max_peptide_length=50, terminal_length=7):
        super(PeptideDetectabilityOptions, self).__init__(
            max_sequence_length=max_peptide_length + terminal_length * 2
        )

        self.amino_acids = '_.' + self.amino_acids
        self.max_aa_length = max_peptide_length
        self.terminal_length = terminal_length


    @staticmethod
    def default():
        return PeptideDetectabilityOptions()

