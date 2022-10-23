import argparse

parser = argparse.ArgumentParser(
    description='Build assays from SpectroMine report or Spectronaut library.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input SpectroMine report or Spectronaut library'
)
parser.add_argument(
    '--out',
    help='output assay file'
)
parser.add_argument(
    '--type', choices=['SpectroMine', 'Spectronaut'], 
    default='SpectroMine',
    help='input report type (default: %(default)s)'
)

filter_group = parser.add_argument_group('entries filters')
filter_group.add_argument(
    '--precursor_charge', type=int, nargs='+', default=[2,3],
    help='list of allowed charge states of precursor ions (default: %(default)s)'
)
filter_group.add_argument(
    '--min_fragment_number', default=6, type=int,
    help='remove assays with < N fragments (default: %(default)s)'
)
filter_group.add_argument(
    '--min_peptide_length', type=int, default=7,
    help='lower sequence length limit of peptides (default: %(default)s)'
)
filter_group.add_argument(
    '--max_peptide_length', type=int, default=50,
    help='upper sequence length limit of peptides (default: %(default)s)'
)
filter_group.add_argument(
    '--modification_action', choices=[
        None, 'keep_if_any', 'keep_if_exclusive',
        'remove_if_any', 'remove_if_exclusive'
    ], default='remove_if_exclusive',
    help='filter precursors according to modifications (default: %(default)s)'
)
filter_group.add_argument(
    '--modification_list', nargs='+',
    default=['Carbamidomethyl'],
    help='selected modifications (default: %(default)s)'
)

filter_group.add_argument(
	'--fragment_type', type=str, nargs='+', default=['b', 'y'],
	help='list of fragment types (default: %(default)s)'
)
filter_group.add_argument(
	'--fragment_charge', type=int, nargs='+', default=[1, 2],
	help='list of allowed charge states of fragment ions (default: %(default)s)'
)
filter_group.add_argument(
	'--fragment_loss_type', type=str, nargs='+', default=['noloss', 'NH3', 'H2O'],
	help='list of neutral loss types of %sfragment ions (default: %(default)s)'
)

args = parser.parse_args()
report_files = getattr(args, 'in')
out_file = args.out
report_type = args.type

filter_args = vars(args)
filter_args.pop('in')
filter_args.pop('out')
filter_args.pop('type')

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('report_files', None) is None:
    report_files = list_files(
        path='.',
        pattern='Fragment( )?Report\\.(csv|xls)$'
    )

if len(report_files) == 0:
    raise ValueError('no report files')

# %%
import os
import re

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(report_files[0])[0]
    out_file = re.sub(
        '[\\._]?Fragment( )?Report$',
        '', out_file
    )
    if len(report_files) > 1:
        out_file += '_' + str(len(report_files))
    out_file += '.assay.pickle'

# %%
import pandas as pd

from util import save_pickle
from assay import AssayBuilder
from assay.table2assay import DataFrameToAssayConverter

if report_type == 'Spectronaut':
    from formatting.spectronaut import \
        Spectronaut_library_parsing_columns as assay_parsing_columns
else:
    from formatting.spectronaut import \
        SpectroMine_assay_parsing_columns as assay_parsing_columns

# %%
logging.info('loading report(s): ' + '; '.join(report_files))

report = pd.concat(
    (
        pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
        for f in report_files
    ),
    ignore_index=True
)

logging.info('report(s) loaded: {0} rows' \
    .format(len(report)))

# %%
converter = DataFrameToAssayConverter(
    columns=assay_parsing_columns()
)

assays = converter.dataframe_to_assays(
    data=report,
    return_generator=True
)

logging.info(
    'filtering assays using the following parameters: \n' + \
    '\n'.join((
        k + '=' + str(v)
        for k, v in filter_args.items()
        if v is not None
    ))
)
assay_builder = AssayBuilder()

assays = assay_builder.filter_assays(
    assays,
    return_generator=True,
    **filter_args
)

logging.info('converting report to assays')

assays = [
    assay for i, assay in enumerate(assays)
    if (i % 1000 == 0 and print(i, end='\t', flush=True)) or True
]

print('')
logging.info('assays converted: {0} spectra' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

