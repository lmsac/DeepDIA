import argparse

parser = argparse.ArgumentParser(
    description='Filter peptides by detectability.'
)
parser.add_argument(
    '--peptide', nargs='+',
    help='input peptide list files'
)
parser.add_argument(
    '--detect', nargs='+',
    help='input peptide detectability files'
)
parser.add_argument(
    '--out', nargs='+',
    help='output peptide list file'
)

parser.add_argument(
    '--min_detectability', type=float, default=0.5,
    help='detectability threshold of peptides  (default: %(default)s)'
)

group_duplicated_group = parser.add_mutually_exclusive_group(required=False)
group_duplicated_group.add_argument(
    '--group_duplicated',
    dest='group_duplicated', action='store_true',
    help='group duplicated peptides (default: %(default)s)'
)
group_duplicated_group.add_argument(
    '--no-group_duplicated',
    dest='group_duplicated', action='store_false',
    help='not group duplicated peptides (default: True)'
)
parser.set_defaults(group_duplicated=True)

args = parser.parse_args()
peptide_files = args.peptide
detectability_files = args.detect
out_files = args.out
min_detectability = args.min_detectability
group_duplicated = args.group_duplicated

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
import os
from util import list_files

if globals().get('peptide_files', None) is None:
    peptide_files = list_files(
        path='.',
        pattern='(?<!detectability.*)\\.peptide\\.csv$',
        recursive=True
    )

if globals().get('detectability_files', None) is None:
    peptide_files = list_files(
        path='.',
        pattern='prediction\\.detectability\\.csv$',
        recursive=True
    )

if len(detectability_files) == 0:
    raise ValueError('no peptide detectability files')

# %%
if globals().get('out_files', None) is None:
    out_files = []
    if peptide_files:
        for peptide_file in peptide_files:
            out_file = os.path.splitext(peptide_file)[0]
            if out_file.endswith('.peptide'):
                out_file = out_file[:-len('.peptide')]
            out_file += '_detectability{0:02.0f}.peptide.csv' \
                .format(min_detectability * 100)
            out_files.append(out_file)
    else:
        for peptide_file in detectability_files:
            out_file = os.path.splitext(peptide_file)[0]
            if out_file.endswith('.detectability'):
                out_file = out_file[:-len('.detectability')]
            if out_file.endswith('.prediction'):
                out_file = out_file[:-len('.prediction')]
            out_file += '_detectability{0:02.0f}.peptide.csv' \
                .format(min_detectability * 100)
            out_files.append(out_file)

if peptide_files and len(out_files) != len(peptide_files):
    raise ValueError('numbers of peptide list files and output files not match')
elif len(out_files) != len(detectability_files):
    raise ValueError('numbers of peptide detectability files and output files not match')

# %%
import pandas as pd
from sequence.digest import group_duplicated_peptides

if peptide_files:
    logging.info('loading peptide detectability: ' + '; '.join(detectability_files))

    detectability = pd.concat(
        (
            pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
            for f in detectability_files
        ),
        ignore_index=True
    )

    logging.info('peptide detectability loaded: {0} entries' \
        .format(len(detectability)))

else:
    detectability = None
    peptide_files = detectability_files


for peptide_file, out_file in zip(peptide_files, out_files):
    logging.info('load peptides: ' + peptide_file)

    peptides = pd.read_csv(peptide_file)

    logging.info('peptides loaded: {0} peptides' \
        .format(len(peptides)))

    if detectability is not None:
        peptides = pd.merge(
            peptides, detectability, 
            on=['sequence', 'nTerminal', 'cTerminal'],
            how='inner', copy=False
        )
    peptides = peptides.loc[
        peptides['detectability'] >= min_detectability
    , :]

    logging.info('peptide list filtered: {0} entries with detectability >= {1}' \
        .format(len(peptides), min_detectability))
    
    if group_duplicated:
        logging.info('grouping duplicated peptides')

        peptides = group_duplicated_peptides(peptides)

        logging.info('duplicated peptides grouped: {0} non-redundant peptides' \
                    .format(len(peptides)))

    logging.info('saving peptide list: {0}' \
        .format(out_file))

    peptides.to_csv(out_file, index=False)

    logging.info('peptide list saved: {0}, {1} entries' \
        .format(out_file, len(peptides)))

