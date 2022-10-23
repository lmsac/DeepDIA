from common.preprocessing import PeptideDataConverter


class PeptideDetectabilityDataConverter(PeptideDataConverter):
    def __init__(self, options):
        super(PeptideDetectabilityDataConverter, self) \
            .__init__(options=options)


    def sequences_to_tensor(self, data):
        sequences = data['nTerminal'] \
            .str.slice(start=-self.options.terminal_length) \
            .str.cat(data['sequence'], sep='.') \
            .str.cat(data['cTerminal'] \
                     .str.slice(stop=self.options.terminal_length), sep='.') \
            .values

        return self.peptides_to_tensor(sequences)


    def data_to_tensor(self, data):
        x = self.sequences_to_tensor(data)
        y = data[['detectability']].values
        return x, y


    def tensor_to_detectability(self, tensor):
        return tensor.flatten()

