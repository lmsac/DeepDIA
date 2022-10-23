import argparse

parser = argparse.ArgumentParser(
    description='Train a model for peptide detectability.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input training data files'
)
parser.add_argument(
    '--model',
    help='pretrained model file'
)
parser.add_argument(
    '--wkdir',
    help='working dir'
)
parser.add_argument(
    '--rounds', type=int, default=5,
    help='number of rounds (default: %(default)s)'
)
parser.add_argument(
    '--positive_threshold', type=float, default=0.5,
    help='score threshold to determine positive data (default: %(default)s)'
)
parser.add_argument(
    '--fpr_threshold', type=float, default=0.05,
    help='FPR threshold for early stopping (default: %(default)s)'
)
parser.add_argument(
    '--epochs', type=int, default=100,
    help='number of epochs in each round (default: %(default)s)'
)
parser.add_argument(
    '--patience', type=int, default=10,
    help='number of patience epochs for early stopping in each round (default: %(default)s)'
)

parser.add_argument(
    '--lr', type=float, default=1e-3,
    help='learning rate (default: %(default)s)'
)
parser.add_argument(
    '--reduce_lr_patience', type=int, default=5,
    help='number of epochs with no improvement after which learning rate will be reduced; ' + \
         'values <= 0 indicate not reducing learning rate (default: %(default)s)'
)
parser.add_argument(
    '--reduce_lr_factor', type=float, default=0.1,
    help='reducing the learning rate by a factor once learning stagnates (default: %(default)s)'
)

parser.add_argument(
    '--validate_percent', type=float, default=0.33,
    help='percent of data for validation in each round (default: %(default)s)'
)
parser.add_argument(
    '--seed', type=int,
    help='run in test mode with fixed seed for data spliting'
)


args = parser.parse_args()
data_files = getattr(args, 'in')
model_file = args.model
working_dir = args.wkdir
seed = args.seed
max_rounds = args.rounds
positive_threshold = args.positive_threshold
fpr_threshold = args.fpr_threshold
epochs = args.epochs
patience = args.patience
lr = args.lr
reduce_lr_factor = args.reduce_lr_factor
reduce_lr_patience = args.reduce_lr_patience
validate_percent = args.validate_percent
seed = args.seed

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
import os
from util import list_files

if globals().get('working_dir', None) is None:
    if globals().get('data_files', None) is None or len(data_files) == 0:
        working_dir = '.'
    else:
        working_dir = os.path.dirname(data_files[0]) or '.'

if globals().get('data_files', None) is None:
    data_files = list_files(
        path=working_dir,
        pattern='(?<!prediction)\\.detectability\\.csv$',
        recursive=True
    )

if len(data_files) == 0:
    raise ValueError('no data files')

# %%
if globals().get('model_file', None) is None:
    model_files = list_files(
        path='.',
        pattern='pretrained.*\\.hdf5$',
        recursive=True
    )

    if len(model_files) > 0:
        model_file = model_files[-1]
    else:
        model_file = None

# %%
from pepdetect import PeptideDetectabilityOptions
from pepdetect.hardnegative import train_hard_negative_iter

options = PeptideDetectabilityOptions.default()

# %%
import pandas as pd
from util import save_json

data_files_positive = [
    f for f in data_files
    if f.find('negative') < 0
]
data_files_negative = [
    f for f in data_files
    if f.find('negative') >= 0
]

if len(data_files_positive) == 0:
    raise ValueError('no positive data files')

if len(data_files_negative) == 0:
    raise ValueError('no negative data files')

data_positive = []
data_negative = []

for data_file, positive in zip(
    data_files_positive + data_files_negative,
    [True] * len(data_files_positive) + [False] * len(data_files_negative)
):
    logging.info('load ' + 'positive' if positive else 'negative' + \
                 ' data: ' + data_file)

    data_ = pd.read_csv(data_file)

    data_ = data_.loc[
        (data_['sequence'].str.len() <= options.max_sequence_length) & \
        data_['sequence'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) & \
        data_['nTerminal'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) & \
        data_['cTerminal'].map(lambda s: \
            all(map(lambda a: a in options.amino_acids, s))) \
    , :]

    if positive:
        data_positive.append(data_)
    else:
        data_negative.append(data_)
    logging.info('{0} data loaded: {1}, {2} valid peptides'.format(
        'positive' if positive else 'negative',
        data_file, len(data_)
    ))

data_positive = pd.concat(data_positive, ignore_index=True)
data_negative = pd.concat(data_negative, ignore_index=True)
logging.info(('data loaded: {0} valid positive peptides and ' + \
              '{1} valid negative peptides totally') \
             .format(len(data_positive), len(data_negative)))

# %%
result_summary = {
    'seed': seed,
    'files': {
        'positive': data_files_positive,
        'negative': data_files_negative
    }
}
save_json(result_summary, os.path.join(working_dir, 'training.json'))

if model_file is not None:
    logging.info('use model: ' + model_file)

results = []

train_iter = train_hard_negative_iter(
    data_positive, data_negative, options=options,
    model_path=model_file,
    working_dir=working_dir,
    max_rounds=max_rounds,
    positive_threshold=positive_threshold,
    fpr_threshold=fpr_threshold,
    seed=seed,
    validate_percent=validate_percent,
    lr=lr,
    reduce_lr_factor=reduce_lr_factor,
    reduce_lr_patience=reduce_lr_patience,
    log=logging.info
)

for i, result in enumerate(train_iter):
    result['files'] = data_files_positive + data_files_negative,

    save_json(result, os.path.join(
        working_dir,
        'training_{i}'.format(i=i),
        'training.json'
    ))

    results.append(result)

result_summary['rounds'] = results
save_json(result_summary, os.path.join(working_dir, 'training.json'))

metric_summary = pd.DataFrame.from_records([
    x['evaluate']
    for x in result_summary['rounds']
])
metric_summary.to_csv(os.path.join(working_dir, 'summary.csv'))
print(metric_summary)

