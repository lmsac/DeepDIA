# coding: utf-8
from __future__ import print_function

import os
import re

from pepms2 import PeptideMS2Options, PeptideMS2Predictor
from pepms2.modeling import build_model_from_weights
from pepms2.utils import save_data_json, load_peptides_csv


options = PeptideMS2Options.default()

predict_dir = '.'

model_path = [ 
    os.path.join(predict_dir, 'models', f) 
    for f in os.listdir(os.path.join(predict_dir, 'models'))
    if re.match(r'^epoch_[0-9]+\.hdf5$', f) is not None
][-1]

predictor = PeptideMS2Predictor(options)
predictor.model = build_model_from_weights(options=options, weights_path=model_path)

data_files = [
    os.path.join(predict_dir, f) 
    for f in os.listdir(predict_dir) 
    if f.endswith('.peptide.csv')
]

for file in data_files:
    peptides = load_peptides_csv(file)
    sequences = peptides['sequence']
    modifications = peptides['modification']
    
    prediction = predictor.predict(sequences, modifications)

    save_data_json(
        data=prediction,
        file=re.sub(r'\.peptide\.csv$', '.prediction.ions.json', file)
    )

        