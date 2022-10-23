import argparse

parser = argparse.ArgumentParser(
    description='Digest protein sequences.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input protein sequences (FASTA)'
)
parser.add_argument(
    '--out',
    help='output peptide list file'
)
parser.add_argument(
    '--fasta_rule', default='UniProt',
    help='FASTA parsing rule (default: %(default)s)'
)
parser.add_argument(
    '--protease', default='Trypsin/P',
    help='protease (default: %(default)s)'
)
parser.add_argument(
    '--max_missed_cleavages', type=int, default=2,
    help='maximum number of missed cleavages (default: %(default)s)'
)
parser.add_argument(
    '--min_peptide_length', type=int, default=7,
    help='lower sequence length limit of peptides  (default: %(default)s)'
)
parser.add_argument(
    '--max_peptide_length', type=int, default=50,
    help='upper sequence length limit of peptides  (default: %(default)s)'
)
parser.add_argument(
    '--min_peptide_mass', type=float, default=0,
    help='lower mass limit of peptides  (default: %(default)s)'
)
parser.add_argument(
    '--max_peptide_mass', type=float, default=4000,
    help='upper mass limit of peptides  (default: %(default)s)'
)
parser.add_argument(
    '--term_window_size', type=int, default=7,
    help='length of terminal windows  (default: %(default)s)'
)

methionine_group = parser.add_mutually_exclusive_group(required=False)
methionine_group.add_argument(
    '--remove_n_term_methionine',
    dest='remove_n_term_methionine', action='store_true',
    help='remove N-terminal methionine (default: %(default)s)'
)
methionine_group.add_argument(
    '--no-remove_n_term_methionine',
    dest='remove_n_term_methionine', action='store_false',
    help='not remove N-terminal methionine (default: False)'
)
parser.set_defaults(remove_n_term_methionine=True)

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
parser.set_defaults(group_duplicated=False)

parser.add_argument(
    '--split_miss', default=False, action='store_true',
    help='split peptides to separate files according to numbers of missed cleavages  (default: %(default)s)'
)

args = parser.parse_args()
fasta_files = getattr(args, 'in')
out_file = args.out
fasta_rule = args.fasta_rule
protease = args.protease
max_missed_cleavages = args.max_missed_cleavages
min_peptide_length = args.min_peptide_length
max_peptide_length = args.max_peptide_length
min_peptide_mass = args.min_peptide_mass
max_peptide_mass = args.max_peptide_mass
remove_n_term_methionine = args.remove_n_term_methionine
term_window_size = args.term_window_size
group_duplicated = args.group_duplicated

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('fasta_files', None) is None:
    fasta_files = list_files(
        path='.',
        pattern='\\.fasta$'
    )

if len(fasta_files) == 0:
    raise ValueError('no protein sequence files')

# %%
logging.info('use fasta_rule: ' + str(fasta_rule))

logging.info('use protease: ' + str(protease))
logging.info('use max_missed_cleavages: ' + str(max_missed_cleavages))

logging.info('use min_peptide_length: ' + str(min_peptide_length))
logging.info('use max_peptide_length: ' + str(max_peptide_length))

logging.info('use min_peptide_mass: ' + str(min_peptide_mass))
logging.info('use max_peptide_mass: ' + str(max_peptide_mass))

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(fasta_files[0])[0]
    if len(fasta_files) > 1:
        out_file += '_' + str(len(fasta_files))
    out_file += '.peptide.csv'

# %%
from sequence.fasta import read_fasta

proteins = []

for fasta_file in fasta_files:
    logging.info('loading protein sequences: ' + fasta_file)

    proteins_ = read_fasta(fasta_file, parsing_rule=fasta_rule)
    proteins.extend(proteins_)

    logging.info('protein sequences loaded: {0}, {1} sequences' \
        .format(fasta_file, len(proteins_)))

logging.info('protein sequences loaded: {0} sequences totally' \
    .format(len(proteins)))

# %%
from sequence.digest import digest, group_duplicated_peptides, \
    filter_peptides_by_mass

logging.info('digesting protein sequences')

peptides = digest(
    proteins, protease=protease,
    max_missed_cleavages=max_missed_cleavages,
    min_peptide_length=min_peptide_length,
    max_peptide_length=max_peptide_length,
    remove_n_term_methionine=remove_n_term_methionine,
    term_window_size=term_window_size
)

logging.info('protein sequences digested: {0} peptides'.format(len(peptides)))

# %%
if min_peptide_mass is not None or max_peptide_mass is not None:
    logging.info('calculating peptide masses')

    peptides = filter_peptides_by_mass(
        peptides,
        min_peptide_mass=min_peptide_mass,
        max_peptide_mass=max_peptide_mass
    )

    logging.info('peptides filtered: {0} remaining peptides' \
                 .format(len(peptides)))

# %%
if group_duplicated:
    logging.info('grouping duplicated peptides')

    peptides = group_duplicated_peptides(peptides)

    logging.info('duplicated peptides grouped: {0} non-redundant peptides' \
                 .format(len(peptides)))

# %%
if globals().get('split_miss', False):
    for miss in peptides['missedCleavages'].unique():
        peptides_miss = peptides.loc[peptides['missedCleavages'] == miss]

        out_file_miss = os.path.splitext(out_file)[0]
        if out_file_miss.endswith('.peptide'):
            out_file_miss = out_file_miss[:-len('.peptide')]
        out_file_miss += '_miss' + str(miss) + '.peptide.csv'

        logging.info('saving peptide list: {0}, {1} missed cleavages' \
            .format(out_file_miss, miss))

        peptides_miss.to_csv(out_file_miss, index=False)

        logging.info('peptide list saved: {0}, {1} missed cleavages, {2} entries' \
            .format(out_file_miss, miss, len(peptides_miss)))

else:
    logging.info('saving peptide list: {0}'\
        .format(out_file))

    peptides.to_csv(out_file, index=False)

    logging.info('peptide list saved: {0}, {1} entries' \
        .format(out_file, len(peptides)))
