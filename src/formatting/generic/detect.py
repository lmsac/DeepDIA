import numpy as np
import pandas as pd

from .table import PeptideReportCleaner
from sequence.digest import digest, filter_peptides_by_mass
from pepmass import PeptideMassCalculator


class PeptideDetectabilityReportCleaner:
    def __init__(self,
                 protein_report_columns=None,
                 peptide_report_columns=None,
                 mass_calculator=None,
                 min_protein_coverage=25,
                 min_observed_peptides=2,
                 protease='Trypsin',
                 max_missed_cleavages=1,
                 min_peptide_length=7,
                 max_peptide_length=50,
                 min_peptide_mass=0,
                 max_peptide_mass=4500,
                 remove_n_term_methionine=True,
                 term_window_size=7):
        if protein_report_columns is None:
            protein_report_columns = default_protein_report_columns()
        if peptide_report_columns is None:
            peptide_report_columns=default_peptide_report_columns()
        if mass_calculator is None:
            mass_calculator = PeptideMassCalculator()

        self.protein_cleaner = PeptideReportCleaner(protein_report_columns)
        self.peptide_cleaner = PeptideReportCleaner(peptide_report_columns)
        self.mass_calculator = mass_calculator

        self.min_protein_coverage = min_protein_coverage
        self.min_observed_peptides = min_observed_peptides
        self.protease = protease
        self.max_missed_cleavages = max_missed_cleavages
        self.min_peptide_length = min_peptide_length
        self.max_peptide_length = max_peptide_length
        self.min_peptide_mass = min_peptide_mass
        self.max_peptide_mass = max_peptide_mass
        self.remove_n_term_methionine = remove_n_term_methionine
        self.term_window_size = term_window_size


    def parse_reports(self, protein_report, peptide_report,
                      supplement_runs=None, supplement_run_pattern=None):
        peptide_report = self.peptide_cleaner.parse_report(peptide_report)
        peptide_report = peptide_report.loc[
            peptide_report['quantity'].gt(0)
        ]

        supplement = None
        if isinstance(supplement_runs, (list, tuple, set, pd.Series)):
            supplement = peptide_report['run'].isin(supplement_runs)
        elif supplement_runs is not None:
            supplement = peptide_report['run'].eq(supplement_runs)

        if supplement_run_pattern is not None:
            supplement_ = peptide_report['run'].str.match(supplement_run_pattern)
            if supplement is not None:
                supplement |= supplement_
            else:
                supplement = supplement_

        if supplement is not None:
            peptide_report_supplement = peptide_report.loc[supplement]
            peptide_report = peptide_report.loc[~supplement]

            peptide_report_supplement = self.peptide_cleaner.remove_duplicates(
                peptide_report_supplement
            )

        peptide_report = self.peptide_cleaner.remove_duplicates(peptide_report)

        if supplement is not None:
            peptide_report_supplement = peptide_report_supplement.loc[
                ~peptide_report_supplement['sequence']
                    .isin(peptide_report['sequence'])
            ]

        protein_report = self.protein_cleaner.parse_report(protein_report)
        protein_report = self.protein_cleaner.remove_duplicates(protein_report)

        if supplement is not None:
            return protein_report, peptide_report, peptide_report_supplement
        else:
            return protein_report, peptide_report


    def filter_peptides(self, protein_report, peptide_report,
                        peptide_report_supplement=None):
        peptide_report = self.peptide_cleaner.filter_peptides(
            peptide_report,
            min_peptide_length=self.min_peptide_length,
            max_peptide_length=self.max_peptide_length
        )
        if self.min_peptide_mass is not None or \
            self.max_peptide_mass is not None:
            peptide_report = filter_peptides_by_mass(
                peptide_report,
                min_peptide_mass=self.min_peptide_mass,
                max_peptide_mass=self.max_peptide_mass
            )

        if peptide_report_supplement is not None:
            peptide_report_supplement = self.peptide_cleaner.filter_peptides(
                peptide_report_supplement,
                min_peptide_length=self.min_peptide_length,
                max_peptide_length=self.max_peptide_length
            )
            if self.min_peptide_mass is not None or \
                self.max_peptide_mass is not None:
                peptide_report_supplement = filter_peptides_by_mass(
                    peptide_report_supplement,
                    min_peptide_mass=self.min_peptide_mass,
                    max_peptide_mass=self.max_peptide_mass
                )

        if self.min_protein_coverage is not None and \
            'coverage' in protein_report.columns:
            protein_report = protein_report.loc[
                protein_report['coverage'] >= self.min_protein_coverage
            ]

        if self.min_observed_peptides is not None:
            if peptide_report_supplement is not None:
                peptide_count = pd.concat((
                    peptide_report['protein'],
                    peptide_report_supplement['protein']
                )).value_counts()
            else:
                peptide_count = peptide_report['protein'].value_counts()

            filtered_peptide = peptide_count.loc[
                peptide_count >= self.min_observed_peptides
            ].index
            peptide_report = peptide_report.loc[
                peptide_report['protein'].isin(filtered_peptide)
            ]

            if peptide_report_supplement is not None:
                peptide_report_supplement = peptide_report_supplement.loc[
                    peptide_report_supplement['protein'] \
                        .isin(filtered_peptide)
                ]

        peptide_report = peptide_report.loc[
            peptide_report['protein'].isin(protein_report['protein'])
        ]
        if peptide_report_supplement is not None:
            peptide_report_supplement = peptide_report_supplement.loc[
                peptide_report_supplement['protein'] \
                    .isin(protein_report['protein'])
            ]

        if peptide_report_supplement is not None:
            filtered_protein = pd.concat((
                peptide_report['protein'],
                peptide_report_supplement['protein']
            )).unique()
        else:
            filtered_protein = peptide_report['protein'].unique()

        protein_report = protein_report.loc[
            protein_report['protein'].isin(filtered_protein)
        ]

        if peptide_report_supplement is not None:
            return protein_report, peptide_report, peptide_report_supplement
        else:
            return protein_report, peptide_report


    def calculate_detectability_score(self, peptide_report,
                                     peptide_report_supplement=None):
        relative_quantity = peptide_report \
            .groupby(by='protein')['quantity'] \
            .transform(lambda x: x / x.max())

        peptide_report = peptide_report.assign(
            detectability=np.maximum(
                (np.log10(relative_quantity) + 5) / 5 * 0.5,
                0
            ) + 0.5
        )

        if peptide_report_supplement is not None:
            peptide_report_supplement = peptide_report_supplement.assign(
                detectability= \
                    (np.min(peptide_report['detectability']) + 0.5) / 2
            )

            peptide_report = pd.concat(
                (peptide_report, peptide_report_supplement),
                ignore_index=True
            )

        return peptide_report


    def match_theoretical_peptides(self, peptide_report, protein_sequences):
        protein = set(peptide_report['protein'])

        digested_peptide = digest(
            filter(lambda entry: entry['id'] in protein, protein_sequences),
            protease=self.protease,
            max_missed_cleavages=self.max_missed_cleavages,
            min_peptide_length=self.min_peptide_length,
            max_peptide_length=self.max_peptide_length,
            remove_n_term_methionine=self.remove_n_term_methionine,
            term_window_size=self.term_window_size
        )
        if self.min_peptide_mass is not None or \
            self.max_peptide_mass is not None:
            digested_peptide = filter_peptides_by_mass(
                digested_peptide,
                min_peptide_mass=self.min_peptide_mass,
                max_peptide_mass=self.max_peptide_mass
            )

        if len(digested_peptide) == 0:
            raise ValueError('no proteins matched')

        positive_peptide = pd.merge(
            digested_peptide,
            peptide_report[['sequence', 'protein', 'detectability']],
            on=['sequence', 'protein']
        )
        positive_peptide.drop_duplicates(subset=['sequence'], inplace=True)

        negative_peptide = digested_peptide \
            .loc[
                ~digested_peptide['sequence'].isin(peptide_report['sequence'])
            ] \
            .drop_duplicates(subset=['sequence'])
        negative_peptide = negative_peptide.assign(detectability=0)

        return positive_peptide, negative_peptide



def default_peptide_report_columns():
    columns = {
        'key': {
            'peptideSequence': {
                'name': 'sequence'
            }
        },
        'info': {
            'rawFile': {
                'name': 'run'
            }
        },
        'agg': {
            'protein': {
                'parse': lambda x: x and x.split(';')[0],
                'action': 'first'
            },
            'quantity': {
                'action': 'mean'
            }
        }
    }

    return columns


def default_protein_report_columns():
    columns = {
        'key': {
            'protein': {
                'parse': lambda x: x and x.split(';')[0]
            }
        },
        'score': {
            'coverage': {
                'parse': lambda x: x and float(x.split(';')[0].strip('%')),
                'ascending': False
            }
        }
    }

    return columns