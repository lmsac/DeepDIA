# coding: utf-8
from __future__ import print_function

import pandas as pd

from pep_detectability import PeptideDetectabilityTrainer

import json
import os


def load_data(file):
    with open(file, 'r') as f:
        data = pd.read_csv(file)
    return data

def save_data_json(data, file):
    with open(file, 'w') as f:
        json.dump(data, f)


working_dir = '.'
filenames_positive = [        
    f for f in os.listdir(working_dir) 
    if f.endswith(".detectability.csv") and \
    not f.endswith("negative.detectability.csv")
]
data_positive = pd.concat([load_data(os.path.join(working_dir, f)) for f in filenames_positive])

filenames_negative = [        
    f for f in os.listdir(working_dir) 
    if f.endswith("negative.detectability.csv")
]
data_negative = pd.concat([load_data(os.path.join(working_dir, f)) for f in filenames_negative])
print('Positive {n}, negative {m}'.format(
    n=len(data_positive), m=len(data_negative)
))

os.makedirs(os.path.join(working_dir, 'models'), exist_ok=True)
trainer = PeptideDetectabilityTrainer(
    save_path=os.path.join(working_dir, 'models', 'epoch_{epoch:03d}.hdf5'),
    log_path=os.path.join(working_dir, 'training.log')
)
result = trainer.train(data_positive, data_negative)
result['files'] = {
    'positive': filenames_positive,
    'negative': filenames_negative
}

trainer.save_model(os.path.join(working_dir, 'models', 'last_epoch.hdf5'))

save_data_json(result, os.path.join(working_dir, 'training.json'))

