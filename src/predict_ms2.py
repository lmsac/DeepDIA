import argparse

parser = argparse.ArgumentParser(
    description='Predict peptide MS2.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input peptide list files'
)
parser.add_argument(
    '--model',
    help='model file'
)
parser.add_argument(
    '--charge', type=int,
    help='precursor charge'
)
parser.add_argument(
    '--out', nargs='+',
    help='output ions files'
)
parser.add_argument(
    '--score', action='store_true', default=False,
    help='calculate similarties between predicted and experimental spectra (input files must contain experimental ions)'
)


args = parser.parse_args()
peptide_files = getattr(args, 'in')
model_file = args.model
charge = args.charge
out_files = args.out
score = args.score

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('peptide_files', None) is None:
    peptide_files = list_files(
        path='.',
        pattern='\\.peptide\\.csv$' if not globals().get('score', False) \
            else '(?<!prediction)\\.ions\\.json$'
    )

if len(peptide_files) == 0:
    raise ValueError('no peptide list files')

# %%
import os

if globals().get('model_file', None) is None:
    model_files = list_files(
        path='models' if os.path.isdir('models') else '.',
        pattern='^epoch_[0-9]+\\.hdf5$',
        recursive=True
    )

    if len(model_files) == 0:
        raise ValueError('no model file')

    model_file = model_files[-1]

# %%
import re

if globals().get('out_files', None) is None:
    out_files = []
    for peptide_file in peptide_files:
        out_file = os.path.splitext(peptide_file)[0]
        if out_file.endswith('.peptide'):
            out_file = out_file[:-len('.peptide')]
        elif out_file.endswith('.ions'):
            out_file = out_file[:-len('.ions')]
            out_file = re.sub('_charge[0-9]+$', '', out_file)
        out_file += ('_charge' + str(charge) if charge is not None else '') + \
            '.prediction.ions.json'
        out_files.append(out_file)

if len(out_files) != len(peptide_files):
    raise ValueError('numbers of peptide list files and output files not match')


# %%
import pandas as pd

from pepms2 import PeptideMS2Predictor, PeptideMS2Options
from util import save_json


options = PeptideMS2Options.default()


logging.info('use model: ' + model_file)

predictor = PeptideMS2Predictor(
    options=options,
    model_path=model_file
)


# %%
import itertools
from assay.similarity import dot_product
from util import load_json

def calculate_similarity(ions1, ions2):
    keys1 = set(ions1['ions'].keys())
    keys2 = set(ions2['ions'].keys())
    keys_common = keys1.intersection(keys2)
    
    intensity1 = list(itertools.chain.from_iterable(
        ions1['ions'][k]
        for k in keys_common
    ))
    intensity2 = list(itertools.chain.from_iterable(
        ions2['ions'][k]
        for k in keys_common
    ))
    intensity1 += list(itertools.chain.from_iterable(
        ions1['ions'][k]
        for k in keys1.difference(keys_common)
    ))
    intensity2 += [0.0] * (len(intensity1) - len(intensity2))
    intensity2 += list(itertools.chain.from_iterable(
        ions2['ions'][k]
        for k in keys2.difference(keys_common)
    ))
    intensity1 += [0.0] * (len(intensity2) - len(intensity1))

    return dot_product(intensity1, intensity2)
    

# %%
for peptide_file, out_file in zip(peptide_files, out_files):
    if not globals().get('score', False) and not peptide_file.endswith('ions.json'):
        logging.info('load peptides: ' + peptide_file)
        
        peptides = pd.read_csv(peptide_file)

        peptides = peptides.loc[
            (peptides['sequence'].str.len() <= options.max_sequence_length) & \
            peptides['sequence'].map(lambda s: \
                all(map(lambda a: a in options.amino_acids, s))) \
        , :]

        logging.info('peptides loaded: {0} valid peptides' \
                    .format(len(peptides)))

        sequences = peptides['sequence']
        modifications = peptides.get('modification', None)

    else:
        logging.info('load peptide ions: ' + peptide_file)

        ions = load_json(peptide_file)

        ions = list(filter(
            lambda d: len(d['peptide']) <= options.max_sequence_length and \
                all(map(lambda a: a in options.amino_acids, d['peptide'])),
            ions
        ))

        logging.info('peptide ions loaded: {0}, {1} valid peptides with charge {2}+' \
                    .format(peptide_file, len(ions), charge))
        
        sequences = [d['peptide'] for d in ions]
        modifications = [d.get('modification', None) for d in ions]

    logging.info('predict peptide MS2: ' + peptide_file)
    
    prediction = predictor.predict(sequences, modifications)

    logging.info('peptide MS2 predicted: {0} spectra' \
                 .format(len(prediction)))
    
    if charge is not None:
        for pred in prediction:
            pred['charge'] = charge

    logging.info('saving peptide MS2: {0}' \
                 .format(out_file))

    save_json(prediction, out_file)

    logging.info('peptide MS2 saved: {0}, {1} spectra' \
        .format(out_file, len(prediction)))

    if globals().get('score', False):
        logging.info('calculate peptide MS2 similarity: ' + out_file)

        scores = pd.DataFrame.from_dict({
            'sequence': sequences,
            'modification' : modifications,
            'charge': [d['charge'] for d in ions],
            'similarity': [
                calculate_similarity(d, pred)
                for d, pred in zip(ions, prediction)
            ]
        })

        out_score_file = os.path.splitext(out_file)[0]
        if out_score_file.endswith('.ions'):
            out_score_file = out_score_file[:-len('.ions')]
        out_score_file += '.ions_score.csv' 

        scores.to_csv(out_score_file, index=False)

        logging.info('scores saved: {0}, {1} peptides' \
            .format(out_score_file, len(scores)))

        logging.info('intensity similarity: median={0}, quantile=({1}, {2})'.format(
            scores['similarity'].median(),
            scores['similarity'].quantile(0.25),
            scores['similarity'].quantile(0.75)
        ))

