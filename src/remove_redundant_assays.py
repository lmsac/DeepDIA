import argparse

parser = argparse.ArgumentParser(
    description='Remove redundant assays.'
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
    '--action', choices=['consensus', 'best', 'first'], default='best',
    help='build consensus assays, keep the best replicates, or keep the first occurrence (default: %(default)s)'
)
parser.add_argument(
    '--score', default='score', 
    help='score used for ranking replicates (when "--action best" is set) or used as weights (when "--action consensus" is set) (default: %(default)s)'
)
parser.add_argument(
    '--score_ascending', default=False, action='store_true',
    help='lower scores are better, only valid when "--action best" is set (default: %(default)s)'
)

withinrun_group = parser.add_mutually_exclusive_group(required=False)
withinrun_group.add_argument(
    '--within_run', 
    dest='within_run', action='store_true',
    help='remove redundant assays within each run (default: %(default)s)'
)
withinrun_group.add_argument(
    '--across_run', 
    dest='within_run', action='store_false',
    help='remove redundant assays across all runs (default: True)'
)
parser.set_defaults(within_run=False)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
action = args.action
score = args.score
score_ascending = args.score_ascending
within_run = args.within_run

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
if globals().get('action', None) is None:
    action = 'consensus'    
    
logging.info('use action: ' + str(action))

if action == 'consensus' or action == 'best':
    if globals().get('score', None) is None:
        score = 'score'
    if globals().get('score_ascending', None) is None:
        score_ascending = False

    if action == 'consensus' and score_ascending:
        raise ValueError('ascending score not supported when "--action consensus" is set')

    logging.info('use score: ' + str(score) + \
        (' (ascending)' if score_ascending else ''))
    
# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    if action == 'consensus':
        out_file += '_consensus.assay.pickle'
    else:
        out_file += '_nonredundant.assay.pickle'

# %%
from util import save_pickle, load_pickle
from assay.combine import peptide_group_key

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
group_key = peptide_group_key(
    within_run=within_run
)

# %%
if action == 'consensus':
    from assay.consensus import ConsensusAssayCombiner
    combiner = ConsensusAssayCombiner(
        group_key=group_key, 
        replicate_weight=score
    )
elif action == 'best':
    from assay.combine import BestReplicateAssayCombiner
    combiner = BestReplicateAssayCombiner(
        group_key=group_key,
        score=score
    )
else:
    from assay.combine import AssayCombiner
    combiner = AssayCombiner(
        group_key=group_key
    )
    
# %%
if action == 'consensus':
    logging.info('building consensus assays')
else:
    logging.info('removing redundant assays')

assays = combiner.remove_redundant(assays)

logging.info('redundant assays removed: {0} spectra remaining' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)
    
logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))
