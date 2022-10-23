from .options import PeptideMS2Options
from .preprocessing import PeptideMS2DataConverter
from .modeling import build_model # , load_model as _load_model


class PeptideMS2Predictor:
    def __init__(self, options=PeptideMS2Options.default(),
                 model_path=None, model=None):
        self.options = options
        self.converter = PeptideMS2DataConverter(self.options)
        if model_path is not None:
            self.load_model(model_path)
        else:
            self.model = model


    def load_model(self, model_path, **kwargs):
        # self.model = _load_model(model_path, **kwargs)
        model = build_model(self.options, **kwargs)
        model.load_weights(model_path)
        self.model = model


    def predict(self, sequences, modifications=None):
        x = self.converter.peptides_to_tensor(sequences, modifications)
        y = self.model.predict(x)
        pred = self.converter.tensor_to_ions(
            y, [len(seq) for seq in sequences]
        )

        if modifications is not None:
            result = [
                {
                    'peptide': seq,
                    'modification': mod,
                    'ions': ions
                } for seq, ions, mod in zip(sequences, pred, modifications)
            ]
        else:
            result = [
                {
                    'peptide': seq,
                    'modification': None,
                    'ions': ions
                } for seq, ions in zip(sequences, pred)
            ]
        return result

