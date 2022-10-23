from .options import PeptideDetectabilityOptions
from .preprocessing import PeptideDetectabilityDataConverter
from .modeling import build_model


class PeptideDetectabilityPredictor:

    def __init__(self, options=PeptideDetectabilityOptions.default(),
                 model_path=None, model=None):
        self.options = options
        self.converter = PeptideDetectabilityDataConverter(self.options)
        if model_path is not None:
            self.load_model(model_path)
        else:
            self.model = model


    def load_model(self, model_path, **kwargs):
        model = build_model(self.options, **kwargs)
        model.load_weights(model_path)
        self.model = model


    def predict(self, data):
        x = self.converter.sequences_to_tensor(data)
        y = self.model.predict(x)
        pred = self.converter.tensor_to_detectability(y)

        result = data[['sequence', 'nTerminal', 'cTerminal']]
        result = result.assign(detectability=pred)
        return result

