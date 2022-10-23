import argparse

parser = argparse.ArgumentParser(
    description='Filter assays.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)

parser.add_argument(
    '--swath_windows',
    help='SWATH isolation window file'
)

assay_filter_group = parser.add_argument_group('assay filters')
assay_filter_group.add_argument(
    '--min_precursor_mz', type=float,
    help='lower m/z limit of precursor ions  (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--max_precursor_mz', type=float,
    help='upper m/z limit of precursor ions  (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--precursor_charge', type=int, nargs='+', default=[2,3],
    help='list of allowed charge states of precursor ions (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--min_fragment_number', default=6, type=int,
    help='remove assays with < N fragments (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--min_peptide_length', type=int, default=7,
    help='lower sequence length limit of peptides  (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--max_peptide_length', type=int, default=50,
    help='upper sequence length limit of peptides  (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--modification_action', choices=[
        None, 'keep_if_any', 'keep_if_exclusive',
        'remove_if_any', 'remove_if_exclusive'
    ], default='remove_if_exclusive',
    help='filter precursors according to modifications (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--modification_list', nargs='+',
    default=['Carbamidomethyl'],
    help='selected modifications (default: %(default)s)'
)

frag_filter_group = parser.add_argument_group('fragment filters')
frag_filter_group.add_argument(
	'--max_fragment_number', type=int,
	help='maximal number of fragments (default: %(default)s)'
)
frag_filter_group.add_argument(
	'--fragment_type', type=str, nargs='+', default=['b', 'y'],
	help='list of fragment types (default: %(default)s)'
)
frag_filter_group.add_argument(
	'--min_fragment_amino_acid_number', type=int,
	help='lower limit of amino acid number of fragment ions (default: %(default)s)'
)
frag_filter_group.add_argument(
	'--fragment_charge', type=int, nargs='+', default=[1, 2],
	help='list of allowed charge states of fragment ions (default: %(default)s)'
)
frag_filter_group.add_argument(
	'--fragment_loss_type', type=str, nargs='+', default=['noloss', 'NH3', 'H2O'],
	help='list of neutral loss types of %sfragment ions (default: %(default)s)'
)
frag_filter_group.add_argument(
	'--min_fragment_mz', type=float,
	help='lower m/z limit of fragment ions (default: %(default)s)'
)
frag_filter_group.add_argument(
	'--max_fragment_mz', type=float,
	help='upper m/z limit of fragment ions (default: %(default)s)'
)
frag_filter_group.add_argument(
    '--min_relative_fragment_intensity', type=float,
    help='lower relative intensity limit of fragment ions (default: %(default)s)'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
swath_window_file = args.swath_windows

filter_args = vars(args)
filter_args.pop('in')
filter_args.pop('out')
filter_args.pop('swath_windows')

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
    out_file += '_filtered.assay.pickle'

# %%
from util import save_pickle, load_pickle
from assay import AssayBuilder
import pandas as pd

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
if swath_window_file is not None:
    logging.info('loading SWATH windows: ' + swath_window_file)

    swath_windows = pd.read_csv(swath_window_file, sep='\t')

    logging.info('SWATH windows loaded: {0} windows' \
                 .format(len(swath_windows)))
else:
    swath_windows = None

# %%
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
    swath_windows=swath_windows,
    **filter_args
)

logging.info('assays filtered: {0} spectra remaining' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

