# coding: utf-8
from __future__ import print_function
import numpy as np
import os

from .options import PeptideMS2Options
from .training import PeptideMS2Trainer
from .modeling import build_model, cosine_similarity
from .preprocessing import DataConverter
from .utils import load_data_dir, save_data_json


def split_cross_validation(data, k=5, seed=None):
    length = len(data)
    np.random.seed(seed)
    indexes = np.random.permutation(length)
    chunks = np.array_split(indexes, k)
    return [
        (np.delete(indexes, chunks[i]), indexes[chunks[i]])
        for i in range(0, k)
    ]


def train_cross_validation(options, data, k=5, seed=None, 
    round_dir='training_{round}', model_file='epoch_{epoch:03d}.hdf5', save_best_only=True, log_file='training.log',
    epochs=100, patience=15, metrics=[cosine_similarity]):
    
    converter = DataConverter(options=options)
    
    x, y = converter.data_to_tensor(data)
    train_validate_indexes = split_cross_validation(data, k=k, seed=seed)
    
    results = []
    for i in range(0, k):
        print('Round {round}/{k}'.format(round=i, k=k))

        train_dir = round_dir.format(round=i)
        os.makedirs(train_dir, exist_ok=True)
        os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)
        
        train_indexes = train_validate_indexes[0][0]
        validate_indexes = train_validate_indexes[0][1]
        split = {
            'seed': seed,
            'cross_validation_k': k,        
            'cross_validation_nth': i,
            'train': train_indexes.tolist(), 
            'validate': validate_indexes.tolist()
        }

        trainer = PeptideMS2Trainer(
            options=options,
            save_path=os.path.join(train_dir, 'models', model_file),
            log_path=os.path.join(train_dir, log_file)
        )
        trainer.model = build_model(options=options, metrics=metrics)

        result = trainer.train_with_tensor(
            x_train = x[train_indexes], y_train = y[train_indexes], 
            x_validate = x[validate_indexes], y_validate= y[validate_indexes],
            epochs=epochs, patience=patience
        )
        result['split'] = split
        
        trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))

        save_data_json(result, os.path.join(train_dir, 'training.json'))

        results.append(result)

    return {
        'seed': seed,
        'cross_validation_k': k,
        'rounds': results
    }
    
