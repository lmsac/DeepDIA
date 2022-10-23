from common.preprocessing import PeptideDataConverter


def normalize_by_min_max(x, min_value, max_value):
    return (x - min_value) / (max_value - min_value)

def denormalize_by_min_max(x, min_value, max_value):
    return x * (max_value - min_value) + min_value


class PeptideRTDataConverter(PeptideDataConverter):
    def __init__(self, options, value_name='irt'):
        super(PeptideRTDataConverter, self).__init__(options=options)

        self.normalize = normalize_by_min_max
        self.denormalize = denormalize_by_min_max
        self.value_name = value_name
    

    def data_to_tensor(self, data):
        sequences = data['sequence'].values
        modifications = data.get('modification', None)
        rt = data[[self.value_name]].values
        x = self.peptides_to_tensor(sequences, modifications)
        y = self.normalize(
            rt, min_value=self.options.min_value, max_value=self.options.max_value
        )
        return x, y


    def tensor_to_rt(self, tensor):
        return self.denormalize(
            tensor,
            min_value=self.options.min_value,
            max_value=self.options.max_value
        ).flatten()

