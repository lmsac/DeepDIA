from . import models

from keras.models import load_model
import pandas as pd

class PeptideDetectabilityPredictor:

    def __init__(self, model_path=None, model=None):
        if model_path is not None:
            model = load_model(model_path)
        elif model is None:
            model = models.build_model()        
        self.model = model

    def predict(self, data):
        sequences = models.get_peptide_sequence(data)
        x = models.peptide_to_tensor(sequences)
        pred = self.model.predict(x)
        return pd.DataFrame.from_items([
            ('sequence', sequences),
            ('detectability', pred.flatten())
        ])

    def load_weights(self, path=None):
        self.model.load_weights(path)