# coding: utf-8
from __future__ import print_function
import numpy as np
import pandas as pd
import os
import json
import re

from pep_detectability import PeptideDetectabilityTrainer
from pep_detectability import PeptideDetectabilityPredictor


def load_data(file):
    with open(file, 'r') as f:
        data = pd.read_csv(file)
    return data

def save_data_json(data, file):
    with open(file, 'w') as f:
        json.dump(data, f)


def get_initial_negative_subset(data1, data2, ratio=1, seed=None):
    if (data1.shape[0] * ratio >= data2.shape[0]):
        return np.arange(data2.shape[0])
    else:
        np.random.seed(seed)
        indexes = np.random.permutation(data2.shape[0])
        return indexes[:(data1.shape[0] * ratio)]

def evaluate_model(data1, data2, model_path=None, model=None, positive_threshold=0.5):
    predictor = PeptideDetectabilityPredictor(model_path=model_path, model=model)
    prediction1 = predictor.predict(data1)
    prediction2 = predictor.predict(data2)
    tp = np.where(prediction1["detectability"].values >= positive_threshold)[0]
    fn = np.where(prediction1["detectability"].values < positive_threshold)[0]
    fp = np.where(prediction2["detectability"].values >= positive_threshold)[0]
    tn = np.where(prediction2["detectability"].values < positive_threshold)[0]
    tpr = len(tp) / (len(tp) + len(fn))
    fpr = len(fp) / (len(fp) + len(tn))
    precision = len(tp) / (len(tp) + len(fp))
    return tpr, fpr, precision, \
        tp, fn, fp, tn, \
        prediction1[["detectability"]].values, prediction2[["detectability"]].values

def get_hard_negative_subset(data1, data2, model_path=None, model=None, \
    positive_threshold=0.5, min_ratio=None, max_ratio=None):
    tpr, fpr, precision, \
        tp, fn, fp, tn, \
        prediction1, prediction2 = \
        evaluate_model(
            data1, data2, 
            model_path=model_path, model=model, 
            positive_threshold=positive_threshold
        )

    indexes = fp

    if min_ratio is not None and min_ratio > 0:
        min_count = data1.shape[0] * min_ratio
    else:
        min_count = None    
    if max_ratio is not None and max_ratio > 0:
        max_count = data1.shape[0] * max_ratio
    else:
        max_count = None      
    if min_count is not None and \
        max_count is not None and \
        min_count > max_count:
            raise ValueError('min_count > max_count')
    
    if min_count is not None and min_count > len(indexes):
        indexes = np.concatenate(
            (indexes, tn[np.argsort(-prediction2[tn], axis=None)])
        )[:min_count]
    elif max_count < len(indexes):
        indexes = indexes[np.argsort(-prediction2[indexes], axis=None)[:max_count]]
    
    return indexes, tpr, fpr, precision
    

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
max_rounds = 10
positive_threshold = 0.5
fpr_threshold = 0.01

result_summary = {
    'seed': seed,
    'files': {
        'positive': filenames_positive,
        'negative': filenames_negative
    }  
}
save_data_json(result_summary, os.path.join(working_dir, 'training.json'))

results = []
for i in range(0, max_rounds + 1):
    if i == 0:
        model_path=None
        indexes_negative = get_initial_negative_subset(
            data_positive, data_negative, seed=seed
        )        
    else:
        predict_dir = '{working_dir}/training_{i}'.format(working_dir=working_dir, i=i - 1)
        model_path = model_path = [ 
            os.path.join(predict_dir, 'models', f) 
            for f in os.listdir(os.path.join(predict_dir, 'models'))
            if re.match(r'^epoch_[0-9]+\.hdf5$', f) is not None
        ][-1]
        indexes_negative, tpr, fpr, precision = get_hard_negative_subset(
            data_positive, data_negative, 
            model_path=model_path, 
            positive_threshold=positive_threshold,
            min_ratio=1
        )

        print('TPR: {tpr}, FPR: {fpr}, Precision: {precision}'.format(
            tpr=tpr, fpr=fpr, precision=precision
        ))
        if (len(results) > 0):
            results[-1]['evaluate'] = {
                'tpr': tpr,
                'fpr': fpr,
                'precision': precision
            }
        if (fpr <= fpr_threshold):
            print('Early stopping')
            break
        if (i >= max_rounds):
            break
    
    print('Round {i}: train on {n} positive samples, {m} negative samples'.format(
        i=i, n=len(data_positive), m=len(indexes_negative)
    ))

    train_dir = '{working_dir}/training_{i}'.format(working_dir=working_dir, i=i)
    os.makedirs(train_dir, exist_ok=True)
    os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)

    data1 = data_positive
    data2 = data_negative.iloc[indexes_negative]
    
    trainer = PeptideDetectabilityTrainer(
        model_path=model_path,
        save_path=os.path.join(train_dir, 'models', 'epoch_{epoch:03d}.hdf5'),
        log_path=os.path.join(train_dir, 'training.log')
    )
    result = trainer.train(data1, data2, seed=seed)
    trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))
        
    result['split'] = {
        'seed': seed,
        'round': i,
        'train': {
            'positive': result['split']['train']['positive'],
            'negative': indexes_negative \
                [np.array(result['split']['train']['negative'])].tolist()
        },                        
        'validate': {
            'positive': result['split']['validate']['positive'],
            'negative': indexes_negative \
                [np.array(result['split']['validate']['negative'])].tolist()
        }
    }           
    result['files'] = filenames_positive + filenames_negative
    trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))
    save_data_json(result, os.path.join(train_dir, 'training.json'))        
    results.append(result)


result_summary['rounds'] = results
save_data_json(result_summary, os.path.join(working_dir, 'training.json'))
