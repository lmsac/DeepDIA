from . import models

import pandas as pd

class PeptideRTPredictor:

    def __init__(self, model_path=None, model=None, rt_min=-50, rt_max=150):
        if model_path is not None:
            model = models.load_model(model_path)
        elif model is None:
            model = models.build_model()        
        self.model = model
        self.rt_min = rt_min
        self.rt_max = rt_max

    def predict(self, sequences):
        pred = models.predict(self.model, sequences, rt_min=self.rt_min, rt_max=self.rt_max)
        return pd.DataFrame.from_items([
            ('sequence', sequences),
            ('irt', pred.flatten())
        ])

    def load_weights(self, path=None):
        self.model.load_weights(path)