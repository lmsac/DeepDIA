from .options import PeptideMS2Options
from .preprocessing import DataConverter
from .modeling import load_model


def predict(model, sequences, modifications=None):
    x = peptide_to_tensor(sequences, modifications)
    y = model.predict(x)
    return tensor_to_ions(y, [len(seq) for seq in sequences])

class PeptideMS2Predictor:

    def __init__(self, options=PeptideMS2Options.default(), model_path=None, model=None):
        self.options = options
        self.converter = DataConverter(self.options)
        if model_path is not None:
            self.model = load_model(model_path)
        else:
            self.model = model
    
    
    def load_model(self, model_path, **kwargs):
        self.model = load_model(model_path, **kwargs)


    def predict(self, sequences, modifications=None):
        if modifications is None:
            modifications = [None] * len(sequences)
        x = self.converter.peptides_to_tensor(sequences, modifications)
        y = self.model.predict(x)
        pred = self.converter.tensor_to_ions(y, [len(seq) for seq in sequences])
        return [
            {
                'peptide': seq,
                'modification': mod,
                'ions': ions
            } for seq, ions, mod in zip(sequences, pred, modifications)
        ]