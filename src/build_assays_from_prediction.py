import argparse

parser = argparse.ArgumentParser(
    description='Build assays from prediction files.'
)
parser.add_argument(
    '--peptide', nargs='+',
    help='input peptide list files'
)
parser.add_argument(
    '--ions', nargs='+',
    help='input predicted ions files'
)
parser.add_argument(
    '--rt', nargs='+',
    help='input predicted retention time/iRT files'
)
parser.add_argument(
    '--im', nargs='+',
    help='input predicted ion mobility files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)

args = parser.parse_args()
peptide_files = args.peptide
ions_files = args.ions
rt_files = args.rt
im_files = args.im
out_file = args.out

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
        pattern='\\.peptide\\.csv$',
        recursive=True
    )

if globals().get('ions_files', None) is None:
    ions_files = list_files(
        path='.',
        pattern='prediction\\.ions\\.json$',
        recursive=True
    )

if len(ions_files) == 0:
    raise ValueError('no ions files')

if globals().get('rt_files', None) is None:
    rt_files = list_files(
        path='.',
        pattern='prediction\\.irt\\.csv$',
        recursive=True
    )

if globals().get('im_files', None) is None:
    im_files = list_files(
        path='.',
        pattern='prediction\\.ionMobility\\.csv$',
        recursive=True
    )

# %%
import os

if globals().get('out_file', None) is None:
    if peptide_files:
        out_file = os.path.splitext(peptide_files[0])[0]
        if out_file.endswith('.peptide'):
            out_file = out_file[:-len('.peptide')]
        if len(peptide_files) > 1:
            out_file += '_' + str(len(peptide_files))
    else:
        out_file = os.path.splitext(ions_files[0])[0]
        if out_file.endswith('.ions'):
            out_file = out_file[:-len('.ions')]
        if len(ions_files) > 1:
            out_file += '_' + str(len(ions_files))
    out_file += '.assay.pickle'

# %%
import pandas as pd

if peptide_files:
    logging.info('loading peptide list: ' + '; '.join(peptide_files))

    peptides = pd.concat(
        (
            pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
            for f in peptide_files
        ),
        ignore_index=True
    )

    logging.info('peptide list loaded: {0} entries' \
        .format(len(peptides)))

else:
    logging.info('no peptide list files')

    peptides = None

# %%
from util import load_json

ions = []
for ions_file in ions_files:
    logging.info('load ions: ' + ions_file)

    ions_ = load_json(ions_file)

    ions.extend(ions_)
    logging.info('ions loaded: {0}, {1} entries' \
                 .format(ions_file, len(ions_)))

logging.info('ions loaded: {0} entriess totally' \
             .format(len(ions)))

# %%
if rt_files:
    logging.info('loading retention time/iRT: ' + '; '.join(rt_files))

    rt = pd.concat(
        (
            pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
            for f in rt_files
        ),
        ignore_index=True
    )

    logging.info('retention time/iRT loaded: {0} entries' \
        .format(len(rt)))

else:
    logging.info('no retention time/iRT files')

    rt = None

if im_files:
    logging.info('loading ion mobility: ' + '; '.join(im_files))

    im = pd.concat(
        (
            pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
            for f in im_files
        ),
        ignore_index=True
    )

    logging.info('ion mobility loaded: {0} entries' \
        .format(len(im)))

else:
    logging.info('no ion mobility files')

    im = None

# %%
from formatting.generic import IonsToAssayConverter
from formatting.generic.mod import stringify_modification
from assay.values import assign_assays_values

converter = IonsToAssayConverter()

logging.info('converting ions to assays')

assays = converter.ions_to_assays(ions=ions)

logging.info('assays converted: {0} spectra' \
             .format(len(assays)))

if peptides is not None:
    logging.info('assigning peptide info')
    assign_assays_values(
        assays, data=peptides,
        keys=[{'name': 'sequence', 'path': 'peptideSequence'}],
        params=[{'name': 'protein', 'path': ['metadata', 'protein']}]
    )

if rt is not None:
    logging.info('assigning peptide retention time/iRT')
    assign_assays_values(
        assays, data=rt,
        keys=[
            {'name': 'sequence', 'path': 'peptideSequence'}, 
            {'name': 'modification', 'path': 'modification',
             'convert': stringify_modification}
        ],
        params=['rt', {'name': 'irt', 'path': 'iRT'}]
    )

if im is not None:
    logging.info('assigning peptide ion mobility')
    assign_assays_values(
        assays, data=im,
        keys=[
            {'name': 'sequence', 'path': 'peptideSequence'}, 
            {'name': 'modification', 'path': 'modification',
             'convert': stringify_modification},
            {'name': 'charge', 'path': 'precursorCharge'},
        ],
        params=['ionMobility']
    ) 

# %%
from util import save_pickle

logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

