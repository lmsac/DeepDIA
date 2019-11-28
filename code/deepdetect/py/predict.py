# coding: utf-8
from __future__ import print_function

import pandas as pd

from pep_detectability import PeptideDetectabilityPredictor

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

predictor = PeptideDetectabilityPredictor()
predictor.load_weights(model_path)

filenames = [
    os.path.join(predict_dir, f) 
    for f in os.listdir(predict_dir) 
    if f.endswith('.peptide.csv')
    # if f.endswith('.detectability.csv')
]

for file in filenames:
    data = pd.read_csv(file)

    prediction = predictor.predict(data)
    
    with open(re.sub(r'\.peptide\.csv$', '.prediction.detectability.csv', file), 'w') as f:
    # with open(re.sub(r'\.detectability\.csv$', '.prediction.detectability.csv', file), 'w') as f:
        prediction.to_csv(f, index=False)

        