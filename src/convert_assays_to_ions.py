import argparse

parser = argparse.ArgumentParser(
    description='Convert assays to ions.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output ions file'
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
    out_file += '.ions.json'

# %%
from util import save_json, load_pickle

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
from formatting.generic import AssayToIonsConverter

converter = AssayToIonsConverter()

logging.info('converting assays to ions')

ions = converter.assays_to_ions(assays)

logging.info('ions converted: {0} spectra' \
    .format(len(ions)))

# %%
for charge in set((x['charge'] for x in ions)):
    ions_charge = [x for x in ions if x['charge'] == charge]

    out_file_charge = os.path.splitext(out_file)[0]
    if out_file_charge.endswith('.ions'):
        out_file_charge = out_file_charge[:-len('.ions')]
    out_file_charge += '_charge' + str(charge) + '.ions.json'

    logging.info('saving ions: {0}, charge {1}+' \
        .format(out_file_charge, charge, len(ions_charge)))

    save_json(ions_charge, out_file_charge)

    logging.info('ions saved: {0}, charge {1}+, {2} spectra' \
        .format(out_file_charge, charge, len(ions_charge)))

