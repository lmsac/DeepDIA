import argparse

parser = argparse.ArgumentParser(
    description='Extract detectability values from MaxQuant report.'
)
parser.add_argument(
    '--protein', nargs='+',
    help='input MaxQuant protein report'
)
parser.add_argument(
    '--peptide', nargs='+',
    help='input MaxQuant peptide report'
)
parser.add_argument(
    '--fasta', nargs='+',
    help='input protein sequences (FASTA)'
)
parser.add_argument(
    '--out',
    help='output detectability file'
)

parser.add_argument(
    '--supplement_run_pattern',
    help='MS file name pattern to determine supplement runs'
)
parser.add_argument(
    '--supplement_runs', nargs='+',
    help='MS runs used as supplement ' + \
         '(e.g., runs with fractionation can be used as supplement runs ' + \
         'to those without fractionation)'
)

ms2_count_group = parser.add_mutually_exclusive_group(required=False)
ms2_count_group.add_argument(
    '--use_ms2_count',
    dest='use_ms2_count', action='store_true',
    help='calculate detecability based on MS/MS count (default: %(default)s)'
)
ms2_count_group.add_argument(
    '--use_intensity',
    dest='use_ms2_count', action='store_false',
    help='calculate detecability based on intensity (default: True)'
)
parser.set_defaults(use_ms2_count=False)

parser.add_argument(
    '--min_protein_coverage', type=float, default=25,
    help='minimum sequence coverage (%) of proteins (default: %(default)s)'
)
parser.add_argument(
    '--min_observed_peptides', type=int, default=2,
    help='minimum number of identified peptides of each protein (default: %(default)s)'
)

parser.add_argument(
    '--fasta_rule', default='default',
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



args = parser.parse_args()
protein_report_files = args.protein
peptide_report_files = args.peptide
fasta_files = args.fasta
out_file = args.out

supplement_run_pattern = args.supplement_run_pattern
supplement_runs = args.supplement_runs
use_ms2_count = args.use_ms2_count
fasta_rule = args.fasta_rule

report_cleaner_args = vars(args)
[report_cleaner_args.pop(k) for k in [
    'protein', 'peptide', 'fasta', 'out',
    'supplement_run_pattern', 'supplement_runs',
    'use_ms2_count',
    'fasta_rule'
]]


# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('protein_report_files', None) is None:
    protein_report_files = list_files(
        path='.',
        pattern='proteinGroups\\.txt$'
    )

if len(protein_report_files) == 0:
    raise ValueError('no protein report files')

if globals().get('peptide_report_files', None) is None:
    peptide_report_files = list_files(
        path='.',
        pattern='peptides\\.txt$'
    )

if len(peptide_report_files) == 0:
    raise ValueError('no peptide report files')

if globals().get('fasta_files', None) is None:
    fasta_files = list_files(
        path='.',
        pattern='\\.fasta$'
    )

if len(fasta_files) == 0:
    raise ValueError('no protein sequence files')

# %%
logging.info('use fasta_rule: ' + str(fasta_rule))

logging.info('use ms2 count' if use_ms2_count else 'use intensity')

logging.info(
    'use the following parameters: \n' + \
    '\n'.join((
        k + '=' + str(v)
        for k, v in report_cleaner_args.items()
        if v is not None
    ))
)

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(protein_report_files[0])[0]
    if len(protein_report_files) > 1:
        out_file += '_' + str(len(protein_report_files))
    out_file += '.detectability.csv'

# %%
import pandas as pd

logging.info('loading protein report(s): ' + '; '.join(protein_report_files))

protein_report = pd.concat(
    (
        pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
        for f in protein_report_files
    ),
    ignore_index=True
)

logging.info('protein report(s) loaded: {0} rows' \
    .format(len(protein_report)))

logging.info('loading peptide report(s): ' + '; '.join(peptide_report_files))

peptide_report = pd.concat(
    (
        pd.read_csv(f, sep=',' if f.endswith('.csv') else '\t')
        for f in peptide_report_files
    ),
    ignore_index=True
)

logging.info('peptide report(s) loaded: {0} rows' \
    .format(len(protein_report)))

# %%
if supplement_runs is not None:
    logging.info('use supplement_runs: ' + '; '.join(supplement_runs))

if supplement_run_pattern is not None:
    logging.info('use supplement_run_pattern: ' + supplement_run_pattern)

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
from formatting.generic import PeptideDetectabilityReportCleaner
from formatting.maxquant import MaxQuant_protein_report_columns, \
    MaxQuant_peptide_report_columns

cleaner = PeptideDetectabilityReportCleaner(
    protein_report_columns=MaxQuant_protein_report_columns(),
    peptide_report_columns=MaxQuant_peptide_report_columns(
        use_ms2_count=use_ms2_count
    ),
    **report_cleaner_args
)

logging.info('parsing protein and peptide reports')

if supplement_runs is not None or supplement_run_pattern is not None:
    protein_report, peptide_report, peptide_report_supplement = \
        cleaner.parse_reports(
            protein_report, peptide_report,
            supplement_runs=supplement_runs,
            supplement_run_pattern=supplement_run_pattern
        )

    logging.info(
        ('reports parsed: {0} protein entries, ' + \
         '{1} main peptide entries, {2} supplement peptide entries') \
        .format(
            len(protein_report),
            len(peptide_report),
            len(peptide_report_supplement)
        )
    )

else:
    protein_report, peptide_report = cleaner.parse_reports(
        protein_report, peptide_report
    )

    logging.info('reports parsed: {0} protein entries, {1} peptide entries' \
        .format(len(protein_report), len(peptide_report)))

# %%
logging.info('filter proteins and peptides')

if supplement_runs is not None or supplement_run_pattern is not None:
    protein_report, peptide_report, peptide_report_supplement = \
        cleaner.filter_peptides(
            protein_report, peptide_report,
            peptide_report_supplement
        )

    logging.info(
        ('proteins and peptides filtered: {0} protein entries, ' + \
         '{1} main peptide entries, {2} supplement peptide entries') \
        .format(
            len(protein_report),
            len(peptide_report),
            len(peptide_report_supplement)
        )
    )

else:
    protein_report, peptide_report = cleaner.filter_peptides(
        protein_report, peptide_report
    )

    logging.info(
        'proteins and peptides filtered: ' + \
        '{0} protein entries, {1} peptide entries' \
        .format(len(protein_report), len(peptide_report))
    )

# %%
logging.info('calculate detectability')

peptide_report = cleaner.calculate_detectability_score(
    peptide_report,
    peptide_report_supplement \
    if supplement_runs is not None else None
)

logging.info(
    'detectability calculated: {0} peptide entries' \
    .format(len(peptide_report))
)

# %%
logging.info('digest protein sequences')

positive_peptide, negative_peptide = cleaner.match_theoretical_peptides(
    peptide_report, proteins
)

logging.info(
    'protein sequences digested: ' + \
    '{0} positive peptides, {1} negative peptides'
    .format(len(positive_peptide), len(negative_peptide))
)

# %%
logging.info('saving positive peptides: {0}' \
    .format(out_file))

positive_peptide.to_csv(out_file, index=False)

logging.info('positive peptides saved: {0}, {1} entries' \
    .format(out_file, len(positive_peptide)))

out_file_negative = os.path.splitext(out_file)[0]
if out_file_negative.endswith('.detectability'):
    out_file_negative = out_file_negative[:-len('.detectability')]
out_file_negative += '_negative' + '.detectability.csv'

logging.info('saving negative peptides: {0}' \
    .format(out_file_negative))

negative_peptide.to_csv(out_file_negative, index=False)

logging.info('negative peptides saved: {0}, {1} entries' \
    .format(out_file_negative, len(negative_peptide)))

