import numpy as np

from .options import PeptideDetectabilityOptions
from .preprocessing import PeptideDetectabilityDataConverter
from .modeling import build_model
from common import TrainerBase
from common.split import split_train_validate


class PeptideDetectabilityTrainer(TrainerBase):
    def __init__(self, options=PeptideDetectabilityOptions.default(),
                 **kwargs):
        super(PeptideDetectabilityTrainer, self).__init__(
            options=options,
            converter=PeptideDetectabilityDataConverter(options),
            **kwargs
        )


    def build_model(self):
        self.model = build_model(self.options)
        return self.model


    def train(self, data1, data2,
              validate_percent=.33, seed=None,
              **kwargs):
        x1, y1 = self.converter.data_to_tensor(data1)
        x1_train, y1_train, x1_validate, y1_validate, \
            train_indexs1, validate_indexs1 = \
            split_train_validate(
                x1, y1,
                validate_percent=validate_percent,
                seed=seed
            )

        x2, y2 = self.converter.data_to_tensor(data2)
        x2_train, y2_train, x2_validate, y2_validate, \
            train_indexs2, validate_indexs2 = \
            split_train_validate(
                x2, y2,
                validate_percent=validate_percent,
                seed=seed
            )
        split = {
            'validate_percent': validate_percent,
            'seed': seed,
            'train': {
                'positive': train_indexs1.tolist(),
                'negative': train_indexs2.tolist()
            },
            'validate': {
                'positive': validate_indexs1.tolist(),
                'negative': validate_indexs2.tolist()
            }
        }
        x_train = np.concatenate((x1_train, x2_train))
        y_train = np.concatenate((y1_train, y2_train))
        x_validate = np.concatenate((x1_validate, x2_validate))
        y_validate = np.concatenate((y1_validate, y2_validate))

        result = self.train_with_tensor(
            x_train, y_train, x_validate, y_validate,
            **kwargs
        )
        result['split'] = split
        return result

