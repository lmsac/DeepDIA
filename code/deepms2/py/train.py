# coding: utf-8
from __future__ import print_function
import os

from pepms2 import PeptideMS2Options, PeptideMS2Trainer
from pepms2.utils import load_data_dir, save_data_json #, filter_data


options = PeptideMS2Options.default()

train_dir = '.'
data, data_files = load_data_dir(train_dir)
# data, filter_indexs = filter_data(data, max_sequence_length=options.max_sequence_length, threshold=0.1)

os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)

trainer = PeptideMS2Trainer(
    options=options,
    save_path=os.path.join(train_dir, 'models', 'epoch_{epoch:03d}.hdf5'),
    log_path=os.path.join(train_dir, 'training.log')
)
result = trainer.train(data)
result['files'] = data_files
# result['filter'] = filter_indexs

trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))

save_data_json(result, os.path.join(train_dir, 'training.json'))

