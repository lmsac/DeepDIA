# coding: utf-8
from __future__ import print_function
import numpy as np
import pandas as pd
import os
import json

from pep_detectability import PeptideDetectabilityTrainer


def load_data(file):
    with open(file, 'r') as f:
        data = pd.read_csv(file)
    return data

def save_data_json(data, file):
    with open(file, 'w') as f:
        json.dump(data, f)

def split_subsets(data, k=5, seed=None):
    length = len(data)
    np.random.seed(seed)
    indexes = np.random.permutation(length)
    chunks = np.array_split(indexes, k)
    return [
        indexes[chunks[i]]
        for i in range(0, k)
    ]

def balance_data(data1, data2, seed=None):
    ratio = data1.shape[0] / data2.shape[0]
    if (ratio >= 1):
        k = round(ratio)
        if (k > 1):
            indexes1 = split_subsets(data1, k, seed=seed)            
        else:
            indexes1 = [np.arange(0, data1.shape[0])]
        indexes2 = [np.arange(0, data2.shape[0])]
        return (indexes1, indexes2)
    else:
        indexes2, indexes1 = balance_data(data2, data1, seed=seed)
        return (indexes1, indexes2)


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

seed = 0
indexes_positive, indexes_negative = balance_data(data_positive, data_negative, seed=seed)
print('Positive {n} subsets, negative {m} subsets'.format(
    n=len(indexes_positive), m=len(indexes_negative)
))

result_summary = {
    'seed': seed,
    'files': {
        'positive': filenames_positive,
        'negative': filenames_negative
    },
    'subsets': {
        'positive': [np.ndarray.tolist(x) for x in indexes_positive],
        'negative': [np.ndarray.tolist(x) for x in indexes_negative]
    }    
}
save_data_json(result_summary, os.path.join(working_dir, 'subsets.json'))

results = []
for i in range(0, len(indexes_positive)):
    for j in range(0, len(indexes_negative)):
        print('Train on positive subset {i} ({n} samples), negative subset {j} ({m} samples)'.format(
            i=i, n=len(indexes_positive[i]), j=j, m=len(indexes_negative[j])
        ))

        train_dir = '{working_dir}/training_pos{i}_neg{j}'.format(working_dir=working_dir, i=i, j=j)
        os.makedirs(train_dir, exist_ok=True)
        os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)
        
        data1 = data_positive.iloc[indexes_positive[i]] 
        data2 = data_negative.iloc[indexes_negative[j]]

        trainer = PeptideDetectabilityTrainer(
            save_path=os.path.join(train_dir, 'models', 'epoch_{epoch:03d}.hdf5'),
            log_path=os.path.join(train_dir, 'training.log')
        )
        result = trainer.train(data1, data2, seed=seed)
        trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))
        
        result['split'] = {
            'seed': seed,
            'positive_subset': i,        
            'negative_subset': j,
            'train': {
                'positive': indexes_positive[i] \
                    [np.array(result['split']['train']['positive'])].tolist(),
                'negative': indexes_negative[j] \
                    [np.array(result['split']['train']['negative'])].tolist()
            },                        
            'validate': {
                'positive': indexes_positive[i] \
                    [np.array(result['split']['validate']['positive'])].tolist(),
                'negative': indexes_negative[j] \
                    [np.array(result['split']['validate']['negative'])].tolist()
            }            
        }       
        result['files'] = filenames_positive + filenames_negative
        trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))
        save_data_json(result, os.path.join(train_dir, 'training.json'))        
        results.append(result)

result_summary['rounds'] = results
save_data_json(result_summary, os.path.join(working_dir, 'subsets.json'))

