import argparse

parser = argparse.ArgumentParser(
    description='Convert assays to Spectronaut spectral library.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output spectral library file'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('assay_files', None) is None:
    assay_files = list_files(
        path='.',
        pattern='\\.assay\\.pickle$'
    )

if len(assay_files) == 0:
    raise ValueError('no assay files')

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    out_file += '.library.xls'

# %%
from util import load_pickle
from assay.assay2table import AssayToDataFrameConverter
from formatting.spectronaut import Spectronaut_assay_library_columns

# %%
assays = []
for assay_file in assay_files:
    logging.info('loading assays: ' + assay_file)

    assay_data = load_pickle(assay_file)
    assays.extend(assay_data)

    logging.info('assays loaded: {0}, {1} spectra' \
        .format(assay_file, len(assay_data)))

logging.info('assays loaded: {0} spectra totally' \
    .format(len(assays)))

# %%
converter = AssayToDataFrameConverter(
    columns=Spectronaut_assay_library_columns()
)

# %%
logging.info('converting assays to table')

data = converter.assays_to_dataframe(assays)

logging.info('assays converted: {0} transitions' \
    .format(len(data)))

# %%
logging.info('saving table: {0}' \
    .format(out_file))

data.to_csv(
    out_file,
    index=False, 
    sep='\t' if not out_file.endswith('.csv') \
        else ','
)

logging.info('table saved: {0}, {1} transitions' \
    .format(out_file, len(data)))

