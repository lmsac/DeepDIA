from .mod import parse_modification_Spectronaut
from formatting.generic.mod import stringify_modification


def SpectroMine_rt_report_columns():
    columns = {
        'key': {
            'PEP.StrippedSequence': {
                'name': 'sequence'
            },
            'PP.PIMID': {
                'drop': True
            }
        },
        'info': {            
            'PP.PIMID': {
                'name': 'modification',
                'parse': parse_modification_Spectronaut,
                'write': stringify_modification
            },
            'PP.EmpiricalRT': {
                'name': 'rt'
            },
            'PP.iRTEmpirical': {
                'name': 'irt'
            }
        },
        'score': {
            'PSM.Score': {
                'ascending': False,
                'drop': True
            },
            'PSM.Qvalue': {
                'ascending': True,
                'drop': True
            }
        }
    }

    return columns


def SpectroMine_im_report_columns():
    columns = {
        'key': {
            'PEP.StrippedSequence': {
                'name': 'sequence'
            },
            'PP.PIMID': {
                'drop': True
            },
            'PP.Charge': {
                'name': 'charge'
            }
        },
        'info': {            
            'PP.PIMID': {
                'name': 'modification',
                'parse': parse_modification_Spectronaut,
                'write': stringify_modification
            },
            'PSM.IonMobility': {
                'name': 'ionMobility'
            }
        },
        'score': {
            'PSM.Score': {
                'ascending': False,
                'drop': True
            },
            'PSM.Qvalue': {
                'ascending': True,
                'drop': True
            }
        }
    }

    return columns


def SpectroMine_peptide_report_columns(use_evidence_count=False):
    columns = {
        'key': {
            'PEP.StrippedSequence': {
                'name': 'sequence'
            },
        },
        'info': {
            'R.FileName': {
                'name': 'run'
            }
        },
        'agg': {
            'PG.ProteinAccessions': {
                'name': 'protein',
                'parse': lambda x: x and x.split(';')[0],
                'action': 'first'
            },
            'R.FileName': {
                'name': 'run'
            }
        }
    }

    if use_evidence_count:
        columns['agg'].update({
            'PEP.RunEvidenceCount':{
                'name': 'quantity',
                'action': 'sum'
            }
        })
    else:
        columns['agg'].update({
            'PEP.Label-Free Quant': {
                'name': 'quantity',
                'action': 'mean'
            }
        })

    return columns


def SpectroMine_protein_report_columns():
    columns = {
        'key': {
            'PG.ProteinAccessions': {
                'name': 'protein',
                'parse': lambda x: x and x.split(';')[0]
            }
        },
        'score': {
            'PG.Coverage': {
                'name': 'coverage',
                'parse': lambda x: x and float(x.split(';')[0].strip('%')),
                'ascending': False
            }
        }
    }

    return columns


def Spectronaut_rt_report_columns():
    columns = {
        'key': {            
            'EG.PrecursorId': {
                'drop': True
            }
        },
        'info': {            
            'PEP.StrippedSequence': {
                'name': 'sequence'
            },
            'EG.PrecursorId': {
                'name': 'modification',
                'parse': parse_modification_Spectronaut,
                'write': stringify_modification
            },
            'EG.ApexRT': {
                'name': 'rt'
            },
            'EG.iRTEmpirical': {
                'name': 'irt'
            }
        },
        'score': {
            'EG.Cscore': {
                'ascending': False,
                'drop': True
            },
            'EG.Qvalue': {
                'ascending': True,
                'drop': True
            }
        }
    }

    return columns


def Spectronaut_im_report_columns():
    columns = {
        'key': {
            'EG.PrecursorId': {
                'drop': True
            }
        },
        'info': {
            'PEP.StrippedSequence': {
                'name': 'sequence'
            },
            'EG.PrecursorId': {
                'name': 'modification',
                'parse': parse_modification_Spectronaut,
                'write': stringify_modification
            },
            'EG.PrecursorId': {
                'name': 'charge',
                'parse': lambda x: int(x.split('.')[-1])
            },
            'EG.IonMobility': {
                'name': 'ionMobility'
            }
        },
        'score': {
            'EG.Cscore': {
                'ascending': False,
                'drop': True
            },
            'EG.Qvalue': {
                'ascending': True,
                'drop': True
            }
        }
    }

    return columns


def Spectronaut_peptide_report_columns(use_evidence_count=False):
    columns = {
        'key': {
            'PEP.StrippedSequence': {
                'name': 'sequence'
            },
        },
        'info': {
            'R.FileName': {
                'name': 'run'
            }
        },
        'agg': {
            'PG.ProteinAccessions': {
                'name': 'protein',
                'parse': lambda x: x and x.split(';')[0],
                'action': 'first'
            },
            'R.FileName': {
                'name': 'run'
            }
        }
    }

    if use_evidence_count:
        columns['agg'].update({
            'PEP.RunEvidenceCount':{
                'name': 'quantity',
                'action': 'sum'
            }
        })
    else:
        columns['agg'].update({
            'PEP.Quantity': {
                'name': 'quantity',
                'action': 'mean'
            }
        })

    return columns


Spectronaut_protein_report_columns = SpectroMine_protein_report_columns

