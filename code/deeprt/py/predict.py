# coding: utf-8
from __future__ import print_function

import pandas as pd

from peprt import PeptideRTPredictor
from peprt.models import max_sequence_length

import itertools
import json
import os
import re

def readlines(file):
    with open(file) as f:
        content = f.readlines()
    return [x.strip() for x in content]


predict_dir = '.'

model_path = [ 
    os.path.join(predict_dir, 'models', f) 
    for f in os.listdir(os.path.join(predict_dir, 'models'))
    if re.match(r'^epoch_[0-9]+\.hdf5$', f) is not None
][-1]

predictor = PeptideRTPredictor()
predictor.load_weights(model_path)

filenames = [
    os.path.join(predict_dir, f) 
    for f in os.listdir(predict_dir) 
    if f.endswith('.peptide.csv')
    # if f.endswith('.irt.csv')
    # if f.endswith('.sequence.txt')
]

for file in filenames:
    sequences = pd.read_csv(file)['sequence']
    # sequences = readlines(file)

    prediction = predictor.predict(sequences)
    
    with open(re.sub(r'\.peptide\.csv$', '.prediction.irt.csv', file), 'w') as f:
    # with open(re.sub(r'\.irt\.csv$', '.prediction.irt.csv', file), 'w') as f:
    # with open(re.sub(r'\.sequence\.txt$', '.prediction.irt.csv', file), 'w') as f:
        prediction.to_csv(f, index=False)

        