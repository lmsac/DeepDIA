import argparse

parser = argparse.ArgumentParser(
    description='Extract ion mobility values from SpectroMine/Spectronaut report.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input SpectroMine/Spectronaut report'
)
parser.add_argument(
    '--out',
    help='output ion mobility file'
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
    '--min_peptide_length', type=int, default=7,
    help='lower sequence length limit of peptides  (default: %(default)s)'
)
filter_group.add_argument(
    '--max_peptide_length', type=int, default=50,
    help='upper sequence length limit of peptides  (default: %(default)s)'
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
        pattern='PSM( )?Report\\.(csv|xls)$'
    )

if len(report_files) == 0:
    raise ValueError('no report files')

# %%
import os
import re

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(report_files[0])[0]
    out_file = re.sub(
        '[\\._]?PSM( )?Report$',
        '', out_file
    )
    if len(report_files) > 1:
        out_file += '_' + str(len(report_files))
    out_file += '.ionMobility.csv'

# %%
import pandas as pd

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
from formatting.generic import PeptideReportCleaner

if report_type == 'Spectronaut':
    from formatting.spectronaut import \
        Spectronaut_im_report_columns as im_report_columns
else:
    from formatting.spectronaut import \
        SpectroMine_im_report_columns as im_report_columns

cleaner = PeptideReportCleaner(columns=im_report_columns())

logging.info('parsing ion mobility report')

data = cleaner.parse_report(report)

logging.info('ion mobility report parsed: {0} entries' \
    .format(len(data)))

# %%
logging.info('remove duplicated entries')

data = cleaner.remove_duplicates(data)

logging.info('duplicated entries removed: {0} non-redundant entries' \
    .format(len(data)))

# %%
logging.info(
    'filtering entries using the following parameters: \n' + \
    '\n'.join((
        k + '=' + str(v)
        for k, v in filter_args.items()
        if v is not None
    ))
)

data = cleaner.filter_peptides(
    data, **filter_args
)

data = cleaner.finalize(data)

logging.info('entries filtered: {0} entries' \
    .format(len(data)))


# %%
for charge in data['charge'].unique():
    data_charge = data.loc[data['charge'] == charge]

    out_file_charge = os.path.splitext(out_file)[0]
    if out_file_charge.endswith('.ionMobility'):
        out_file_charge = out_file_charge[:-len('.ionMobility')]
    out_file_charge += '_charge' + str(charge) + '.ionMobility.csv'

    logging.info('saving ion mobility report: {0}, charge {1}+' \
        .format(out_file_charge, charge))

    data_charge.to_csv(out_file_charge, index=False)

    logging.info('ion mobility report saved: {0}, charge {1}+, {2} entries' \
        .format(out_file_charge, charge, len(data_charge)))


