from peprt import PeptideRTOptions
from peprt import PeptideRTTrainer
from peprt import PeptideRTPredictor


def ion_mobility_options(max_sequence_length=50, 
                    min_value=0.5, max_value=1.5):
    return PeptideRTOptions(
        max_sequence_length=max_sequence_length,
        min_value=min_value, max_value=max_value
    ) 
    

def ion_mobility_trainer(options=ion_mobility_options(), 
                         **kwargs):
    return PeptideRTTrainer(
        options=options, value_name='ionMobility',
        **kwargs
    )

def ion_mobility_predictor(options=ion_mobility_options(), 
                           **kwargs):
    return PeptideRTPredictor(
        options=options, value_name='ionMobility',
        **kwargs
    )


__all__ = [
    'ion_mobility_options', 
    'ion_mobility_trainer', 
    'ion_mobility_predictor'
]


