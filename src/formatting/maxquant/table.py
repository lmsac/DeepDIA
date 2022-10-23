from .mod import parse_modification_MaxQuant
from formatting.generic.mod import stringify_modification


def MaxQuant_evidence_report_columns(type='rt'):
    columns = {
        'key': {
            'Sequence': {
                'name': 'sequence'
            },
            'Modified sequence': {
                'drop': True
            }            
        },
        'info': {
            'Modified sequence': {
                'name': 'modification',
                'parse': parse_modification_MaxQuant,
                'write': stringify_modification
            },
            'Reverse': {
                'name': 'decoy',
                'parse': lambda x: x == '+'
            },
            'Potential contaminant': {
                'name': 'contaminant',
                'parse': lambda x: x == '+'
            }            
        },
        'score': {
            'Score': {
                'ascending': False,
                'drop': True
            }
        }
    }

    if type == 'rt' or type == 'retentionTime':
        columns['info'].update({
            'Calibrated retention time': {
                'name': 'rt'
            }
        })
    elif type == 'im' or type == 'ionMobility':
        columns['key'].update({
            'Charge': {
                'name': 'charge'
            }
        })
        columns['info'].update({
            'Calibrated 1/K0': {
                'name': 'ionMobility'
            },
            'Calibrated K0': {
                'name': 'ionMobility',
                'parse': lambda x: 1 / x
            }
        })

    return columns


def MaxQuant_peptide_report_columns(use_ms2_count=False):
    columns = {
        'key': {
            'Sequence': {
                'name': 'sequence'
            },
        },
        'info': {
            'Reverse': {
                'name': 'decoy',
                'parse': lambda x: x == '+'
            },
            'Potential contaminant': {
                'name': 'contaminant',
                'parse': lambda x: x == '+'
            }
        },
        'agg': {
            'Proteins': {
                'name': 'protein',
                'parse': lambda x: x and str(x).split(';')[0],
                'action': 'first'
            }
        }
    }

    if use_ms2_count:
        columns['agg'].update({
            'MS/MS Count':{
                'name': 'quantity',
                'action': 'sum'
            }
        })
    else:
        columns['agg'].update({
            'Intensity': {
                'name': 'quantity',
                'action': 'mean'
            }
        })

    return columns


def MaxQuant_protein_report_columns():
    columns = {
        'key': {
            'Protein IDs': {
                'name': 'protein',
                'parse': lambda x: x and str(x).split(';')[0]
            }
        },
        'info': {
            'Reverse': {
                'name': 'decoy',
                'parse': lambda x: x == '+'
            },
            'Potential contaminant': {
                'name': 'contaminant',
                'parse': lambda x: x == '+'
            }
        },
        'score': {
            'Sequence coverage [%]': {
                'name': 'coverage',
                'ascending': False
            }
        }
    }

    return columns