# coding: utf-8
from __future__ import print_function

import pandas as pd

from peprt import PeptideRTTrainer
from peprt.models import max_sequence_length

import json
import os


def load_data(file):
    with open(file, 'r') as f:
        data = pd.read_csv(file)
    return data

def load_data_dir(dir):
    filenames = [        
        f for f in os.listdir(dir) 
        if f.endswith(".irt.csv")
    ]
    data = [load_data(os.path.join(dir, f)) for f in filenames]
    data = pd.concat(data)
    return data, filenames


train_dir = '.'
data, data_files = load_data_dir(train_dir)

os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)
trainer = PeptideRTTrainer(
    save_path=os.path.join(train_dir, 'models', 'epoch_{epoch:03d}.hdf5'),
    log_path=os.path.join(train_dir, 'training.log')
)
result = trainer.train(data)
result['files'] = data_files

trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))

with open(os.path.join(train_dir, 'training.json'), 'w') as f:
    json.dump(result, f)

