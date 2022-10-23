from .options import PeptideRTOptions
from .preprocessing import PeptideRTDataConverter
from .modeling import build_model
from common import TrainerBase


class PeptideRTTrainer(TrainerBase):
    def __init__(self, options=PeptideRTOptions.default(),
                 value_name='irt',
                 **kwargs):
        super(PeptideRTTrainer, self).__init__(
            options=options,
            converter=PeptideRTDataConverter(options, value_name=value_name),
            **kwargs
        )


    def build_model(self):
        self.model = build_model(self.options)
        return self.model

