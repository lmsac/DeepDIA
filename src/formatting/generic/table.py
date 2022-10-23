from collections import OrderedDict
import pandas as pd

from assay.modseq import check_modifications


class PeptideReportCleaner:
    def __init__(self, columns):
        if columns is None:
            raise ValueError('invalid columns')

        self.columns = columns


    def parse_report(self, data,
                     include_decoy=False,
                     include_contaminant=False):
        result = OrderedDict()
        for k, v in self.columns.items():
            for colname, col in v.items():
                if colname in data.columns:
                    value = data[colname]

                    convert = col.get('parse', None)
                    if callable(convert):
                        value = value.map(convert)

                    result[col.get('name', colname)] = value

        result = pd.DataFrame.from_dict(result)

        if not include_decoy and 'decoy' in result.columns:
            result = result.loc[~result['decoy']]
        if not include_contaminant and 'contaminant' in result.columns:
            result = result.loc[~result['contaminant']]

        return result


    def remove_duplicates(self, data):
        result = data

        score =  self.columns.get('score', None)
        if score is not None:
            score_name = [
                col.get('name', colname)
                for colname, col in score.items()
            ]
            ascending = [
                col.get('ascending', False)
                for colname, col in score.items()
            ]
            found = [colname in result.columns for colname in score_name]
            if any(found):
                result = result.sort_values(
                    by=[colname for colname, use in zip(score_name, found) if use],
                    ascending=[asc for asc, use in zip(ascending, found) if use]
                )

        key = self.columns.get('key', None)

        if key is not None:
            key_name = [
                col.get('name', colname)
                for colname, col in key.items()
            ]
            found = [colname in result.columns for colname in key_name]
            if not any(found):
                raise KeyError(key_name)
            key_name = [colname for colname, use in zip(key_name, found) if use]

            agg = self.columns.get('agg', None)
            if agg is not None:
                result = result.groupby(
                    by=key_name,
                    as_index=False
                ).aggregate({
                    col.get('name', colname): col.get('action', 'first')
                    for colname, col in agg.items()
                })

            else:
                result = result.drop_duplicates(
                    subset=key_name
                )

        return result


    def filter_peptides(self, data,
                        precursor_charge=None,
                        min_peptide_length=None,
                        max_peptide_length=None,
                        modification_action=None,
                        modification_list=None):
        result = data

        if precursor_charge is not None and 'charge' in data.columns:
            if isinstance(precursor_charge, list) or \
                isinstance(precursor_charge, tuple) or \
                isinstance(precursor_charge, set):
                drop = ~result['charge'].isin(precursor_charge)
            else:
                drop = ~result['charge'].eq(precursor_charge)

            result = result.loc[~drop, :]

        sequence_length = result['sequence'].str.len()
        drop = sequence_length.le(0)
        if min_peptide_length is not None:
            drop = drop & sequence_length.lt(min_peptide_length)
        if max_peptide_length is not None:
            drop = drop & sequence_length.gt(max_peptide_length)

        result = result.loc[~drop, :]

        if modification_action is None or \
            modification_action == 'None' or \
            modification_action == 'none':
            pass

        elif 'modification' in data.columns:
            if modification_action == 'keep_if_any' or \
                modification_action == 'include':
                action = 'if_any'
            elif modification_action == 'keep_if_exclusive':
                action = 'if_exclusive'
            elif modification_action == 'remove_if_any' or \
                modification_action == 'exclude':
                action = 'not_if_any'
            elif modification_action == 'remove_if_exclusive':
                action = 'not_if_exclusive'
            else:
                action = modification_action

            drop = ~result['modification'].apply(
                check_modifications,
                action=action,
                name_list=modification_list
            )

            result = result.loc[~drop, :]

        return result


    def finalize(self, data):
        result = data

        for k, v in self.columns.items():
            drop = result.columns.intersection([
                col.get('name', colname)
                for colname, col in v.items()
                if col.get('drop', False)
            ])
            if any(drop):
                result = result.drop(
                    drop, axis=1
                )

            for colname, col in v.items():
                convert = col.get('write', None)
                colname = col.get('name', colname)
                if callable(convert) and colname in result.columns:
                    result[colname] = result[colname].map(convert)

        return result



from assay.modseq import parse_modification
from .mod import stringify_modification


def default_rt_report_columns():
    columns = {
        'key': {
            'peptideSequence': {
                'name': 'sequence'
            },
            'modification': {
                'name': 'mod_id',
                'parse': str,
                'drop': True
            }
        },
        'info': {
            'modification': {
                'name': 'modification',
                'parse': parse_modification,
                'write': stringify_modification
            },
            'rt': {
                'name': 'rt'
            },
            'irt': {
                'name': 'irt'
            }
        },
        'score': {
            'score': {
                'ascending': False,
                'drop': True
            }
        }
    }

    return columns


def default_im_report_columns():
    columns = {
        'key': {
            'peptideSequence': {
                'name': 'sequence'
            },
            'modification': {
                'name': 'mod_id',
                'parse': str,
                'drop': True
            },
            'charge': {
                'name': 'charge'
            }
        },
        'info': {
            'modification': {
                'name': 'modification',
                'parse': parse_modification,
                'write': stringify_modification
            },
            'ionMobility': {}
        },
        'score': {
            'score': {
                'ascending': False,
                'drop': True
            }
        }
    }

    return columns


