import itertools
import json
import os
import re
import pandas as pd


def load_data_json(file):
    with open(file, 'r') as f:
        data = json.load(f)
    return data

def load_data_dir(dir, pattern=r'\.ions\.json$'):
    filenames = [
        f for f in os.listdir(dir)
        if re.search(pattern, f) is not None
    ]
    data = [load_data_json(os.path.join(dir, f)) for f in filenames]
    data = list(itertools.chain.from_iterable(data))
    return data, filenames

def filter_data(data, max_sequence_length=50, threshold_field='assigned', threshold=0.1, threshold_upper=True):
    indexs = [
        i for i, d in enumerate(data)
        if len(d['peptide']) <= max_sequence_length
            and (d[threshold_field] >= threshold if threshold_upper else d[threshold_field] <= threshold)
    ]
    filtered_data = [data[i] for i in indexs]
    return filtered_data, indexs

def save_data_json(data, file):
    with open(file, 'w') as f:
        json.dump(data, f)

def readlines(file):
    with open(file) as f:
        content = f.readlines()
    return [x.strip() for x in content]
    
def load_peptides_csv(file, sep=','):
    df = pd.read_csv(file, sep=sep, na_filter=False)
    return {
        "sequence": df['sequence'].values,
        "modification": df['modification'].values 
            if 'modification' in df else None
    }

