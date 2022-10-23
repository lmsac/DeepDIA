import argparse

parser = argparse.ArgumentParser(
    description='Train a model for peptide MS2.'
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
    '--epochs', type=int, default=100,
    help='number of epochs (default: %(default)s)'
)
parser.add_argument(
    '--patience', type=int, default=15,
    help='number of patience epochs for early stopping (default: %(default)s)'
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
    help='percent of data for validation (default: %(default)s)'
)
parser.add_argument(
    '--seed', type=int,
    help='run in test mode with fixed seed for data spliting'
)


args = parser.parse_args()
data_files = getattr(args, 'in')
model_file = args.model
train_dir = args.wkdir
seed = args.seed
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

if globals().get('train_dir', None) is None:
    if globals().get('data_files', None) is None or len(data_files) == 0:
        train_dir = '.'
    else:
        train_dir = os.path.dirname(data_files[0]) or '.'

if globals().get('data_files', None) is None:
    data_files = list_files(
        path=train_dir,
        pattern='(?<!prediction)\\.ions\\.json$',
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
from pepms2 import PeptideMS2Options, PeptideMS2Trainer

options = PeptideMS2Options.default()

# %%
from util import load_json, save_json

data = []
for data_file in data_files:
    logging.info('load data: ' + data_file)

    data_ = load_json(data_file)

    data_ = list(filter(
        lambda d: len(d['peptide']) <= options.max_sequence_length and \
            all(map(lambda a: a in options.amino_acids, d['peptide'])),
        data_
    ))
    data.extend(data_)
    logging.info('data loaded: {0}, {1} valid peptides' \
                 .format(data_file, len(data_)))

logging.info('data loaded: {0} valid peptides totally' \
             .format(len(data)))

# %%
os.makedirs(os.path.join(train_dir, 'models'), exist_ok=True)

trainer = PeptideMS2Trainer(
    options=options,
    save_path=os.path.join(train_dir, 'models', 'epoch_{epoch:03d}.hdf5'),
    log_path=os.path.join(train_dir, 'training.log')
)

if model_file is not None:
    logging.info('use model: ' + model_file)

    trainer.load_model(model_path=model_file)

if lr is not None:
    logging.info('use lr: ' + str(lr))

if reduce_lr_patience:
    trainer.use_reduced_lr(
        patience=reduce_lr_patience,
        factor=reduce_lr_factor
    )

result = trainer.train(
    data, epochs=epochs, patience=patience,
    validate_percent=validate_percent, seed=seed,
    lr=lr
)
result['files'] = data_files

trainer.save_model(os.path.join(train_dir, 'models', 'last_epoch.hdf5'))

save_json(result, os.path.join(train_dir, 'training.json'))

