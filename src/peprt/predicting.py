from .options import PeptideRTOptions
from .preprocessing import PeptideRTDataConverter
from .modeling import build_model

import pandas as pd
from collections import OrderedDict


class PeptideRTPredictor:
    def __init__(self, options=PeptideRTOptions.default(),
                 model_path=None, model=None, value_name='irt'):
        self.options = options
        self.value_name = value_name
        self.converter = PeptideRTDataConverter(
            self.options, 
            value_name=value_name
        )
        if model_path is not None:
            self.load_model(model_path)
        else:
            self.model = model


    def load_model(self, model_path, **kwargs):
        model = build_model(self.options, **kwargs)
        model.load_weights(model_path)
        self.model = model


    def predict(self, sequences, modifications=None):
        x = self.converter.peptides_to_tensor(sequences, modifications)
        y = self.model.predict(x)
        pred = self.converter.tensor_to_rt(y)

        return pd.DataFrame.from_dict(OrderedDict([
            ('sequence', sequences),
            ('modification', modifications if modifications is not None \
                 else [None] * len(sequences)),
            (self.value_name, pred)
        ]))

