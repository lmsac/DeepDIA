# coding: utf-8
from __future__ import print_function

import pandas as pd

from pep_detectability import PeptideDetectabilityPredictor
from pep_detectability.models import max_sequence_length, aa_size

from keras.models import Model, Input
from keras.layers import Average
from keras.models import load_model

import itertools
import json
import os
import re

def ensemble(models):
    model_input = Input(shape=(max_sequence_length, aa_size))
    outputs = [model(model_input) for model in models]
    y = Average()(outputs)
    model = Model(model_input, y, name='ensemble')
    return model

def readlines(file):
    with open(file) as f:
        content = f.readlines()
    return [x.strip() for x in content]


predict_dir = '.'

model_paths = [
    [ 
        os.path.join(predict_dir, train_dir, 'models', f) 
        for f in os.listdir(os.path.join(predict_dir, train_dir, 'models'))
        if re.match(r'^epoch_[0-9]+\.hdf5$', f) is not None
    ][-1]
    for train_dir in os.listdir(os.path.join(predict_dir))
    if re.match(r'^training_', train_dir) is not None and os.path.isdir(os.path.join(predict_dir, train_dir))
]

models = [
    load_model(model_path)
    for model_path in model_paths
]

ensembled_model = ensemble(models=models)
# ensembled_model.save(
#     'ensembled_' + '_'.join([
#         re.search(r'epoch_([0-9]+)\.hdf5$', model_path).group(1)
#         for model_path in model_paths
#     ]) + '.hdf5'
# )

predictor = PeptideDetectabilityPredictor(model=ensembled_model)

filenames = [
    os.path.join(predict_dir, f) 
    for f in os.listdir(predict_dir) 
    # if f.endswith('.peptide.csv')
    if f.endswith('.detectability.csv')
]

for file in filenames:
    data = pd.read_csv(file)

    prediction = predictor.predict(data)
    
    # with open(re.sub(r'\.peptide\.csv$', '.prediction.detectability.csv', file), 'w') as f:
    with open(re.sub(r'\.detectability\.csv$', '.prediction.detectability.csv', file), 'w') as f:
        prediction.to_csv(f, index=False)

        