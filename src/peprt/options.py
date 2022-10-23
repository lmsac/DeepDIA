from common.options import PeptideOptions


class PeptideRTOptions(PeptideOptions):
    def __init__(self, max_sequence_length=50, min_value=-50, max_value=150):
        super(PeptideRTOptions, self).__init__(
            max_sequence_length=max_sequence_length
        )

        self.min_value = min_value
        self.max_value = max_value


    @staticmethod
    def default():
        return PeptideRTOptions()

