from common.options import PeptideOptions


class PeptideMS2Options(PeptideOptions):
    def __init__(self, max_sequence_length=50):
        super(PeptideMS2Options, self).__init__(
            max_sequence_length=max_sequence_length
        )

        self.fragments = [
            ('b1', False), ('b2', False),
            ('bn1', False), ('bn2', False),
            ('bo1', False), ('bo2', False),
            ('y1', True), ('y2', True),
            ('yn1', True), ('yn2', True),
            ('yo1', True), ('yo2', True)
        ]


    def intensity_size(self):
        return len(self.fragments)


    @staticmethod
    def default():
        return PeptideMS2Options()

