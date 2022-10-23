import argparse

parser = argparse.ArgumentParser(
    description='Concatenate spectral library text files.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input spectral library text files'
)
parser.add_argument(
    '--out',
    help='output spectral library textfile'
)

args = parser.parse_args()
data_files = getattr(args, 'in')
out_file = args.out

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('data_files', None) is None:
    data_files = list_files(
        path='.',
        pattern='\\.library\\.(csv|xls|tsv)$'
    )

if len(data_files) == 0:
    raise ValueError('no spectral library text files')

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(data_files[0])[0]
    if out_file.endswith('.library'):
        out_file = out_file[:-len('.library')]
        out_file += '_' + str(len(data_files))
    out_file += '_concat.library.xls'

# %%
import pandas as pd

data = []
for data_file in data_files:
    logging.info('load data: ' + data_file)

    data_ = pd.read_csv(
        data_file, 
        sep=',' if data_file.endswith('.csv') \
            else '\t'
    )

    data.append(data_)
    logging.info('data loaded: {0}, {1} transitions' \
                 .format(data_file, len(data_)))

# %%
data = pd.concat(data, ignore_index=True)
logging.info('library concatenated: {0} transitions totally' \
             .format(len(data)))

# %%
logging.info('saving library: {0}' \
    .format(out_file))

data.to_csv(out_file, sep='\t', index=False)

logging.info('library saved: {0}, {1} transitions' \
    .format(out_file, len(data)))
