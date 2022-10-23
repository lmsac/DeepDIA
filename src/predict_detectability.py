import argparse

parser = argparse.ArgumentParser(
    description='Predict peptide detectability.'
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
    '--out', nargs='+',
    help='output detectability files'
)


args = parser.parse_args()
peptide_files = getattr(args, 'in')
model_file = args.model
out_files = args.out

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
        pattern='\\.peptide\\.csv$'
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
if globals().get('out_files', None) is None:
    out_files = []
    for peptide_file in peptide_files:
        out_file = os.path.splitext(peptide_file)[0]
        if out_file.endswith('.peptide'):
            out_file = out_file[:-len('.peptide')]
        out_file += '.prediction.detectability.csv'
        out_files.append(out_file)

if len(out_files) != len(peptide_files):
    raise ValueError('numbers of peptide list files and output files not match')


# %%
import pandas as pd

from pepdetect import PeptideDetectabilityPredictor, PeptideDetectabilityOptions


options = PeptideDetectabilityOptions.default()


logging.info('use model: ' + model_file)

predictor = PeptideDetectabilityPredictor(
    options=options,
    model_path=model_file
)


# %%
for peptide_file, out_file in zip(peptide_files, out_files):
    logging.info('load peptides: ' + peptide_file)

    peptides = pd.read_csv(peptide_file)

    peptides = peptides.loc[
        (peptides['sequence'].str.len() <= options.max_sequence_length) & \
        peptides['sequence'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) & \
        peptides['nTerminal'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) & \
        peptides['cTerminal'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) \
    , :]

    logging.info('peptides loaded: {0} valid peptides' \
                 .format(len(peptides)))

    logging.info('predict peptide detectability: ' + peptide_file)

    prediction = predictor.predict(peptides)

    logging.info('peptide detectability predicted: {0} peptides' \
                 .format(len(prediction)))

    logging.info('saving peptide detectability: {0}' \
                 .format(out_file))

    prediction.to_csv(out_file, index=False)

    logging.info('peptide detectability saved: {0}, {1} peptides' \
        .format(out_file, len(prediction)))

