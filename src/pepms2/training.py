from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger

from .options import PeptideMS2Options
from .preprocessing import PeptideMS2DataConverter
from .modeling import build_model, load_model
from common import TrainerBase


class PeptideMS2Trainer(TrainerBase):
    def __init__(self, options=PeptideMS2Options.default(),
                 **kwargs):
        super(PeptideMS2Trainer, self).__init__(
            options=options,
            converter=PeptideMS2DataConverter(options),
            **kwargs
        )


    def load_model(self, model_path, **kwargs):
        self.model = load_model(model_path, **kwargs)
        return self.model


    def build_model(self):
        self.model = build_model(self.options)
        return self.model

