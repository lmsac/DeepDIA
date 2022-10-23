import pandas as pd

from .modseq import stringify_modification
from .values import get_assay_values


class AssayToDataFrameConverter:
    def __init__(self, columns=None):
        if columns is None:
            columns = default_assay_table_columns()
        self.columns = columns


    def assay_to_dataframe(self, assay, **func_args):
        d = get_assay_values(assay, self.columns, func_args)

        if not any((
            isinstance(x, list) or isinstance(x, tuple)
            for x in d.values()
        )):
            return pd.DataFrame(d, index=[0])
        else:
            return pd.DataFrame.from_dict(d)


    def assays_to_dataframe(self, assays, **func_args):
        return pd.concat(
            filter(len, (
                self.assay_to_dataframe(x, index=i, **func_args)
                for i, x in enumerate(assays)
            )), 
            ignore_index=True
        )



def default_assay_table_columns():
    info = [
        {
            'name': 'peptideSequence',
            'path': 'peptideSequence'
        },
        {
            'name': 'modification',
            'path': 'modification',
            'convert': stringify_modification,
        },
        {
            'name': 'precursorCharge',
            'path': 'precursorCharge',
        }
    ]

    precursor = [
        {
            'name': 'precursorMz',
            'path': 'precursorMZ'
        },
        {
            'name': 'retentionTime',
            'path': 'rt'
        },
    ]

    fragment = [
        {
            'name': 'fragmentMz',
            'path': ['fragments', 'fragmentMZ']
        },
        {
            'name': 'fragmentIntensity',
            'path': ['fragments', 'fragmentIntensity']
        },
        {
            'name': 'fragmentType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'name': 'fragmentNumber',
            'path': ['fragments', 'fragmentNumber']
        },
        {
            'name': 'fragmentCharge',
            'path': ['fragments', 'fragmentCharge']
        },
        {
            'name': 'fragmentAnnotation',
            'path': ['fragments', 'fragmentAnnotation']
        }
    ]

    metadata = [
        {
            'name': 'protein',
            'path': ['metadata', 'protein']
        },
        {
            'name': 'decoy',
            'path': ['metadata', 'decoy'],
            'default': False
        }
    ]

    return info + precursor + fragment + metadata



if __name__ == '__main__':
    converter = AssayToDataFrameConverter()

    assay0 = {
        'modification': None,
        'precursorCharge': 2,
        'precursorMZ': 1184.6093373850001,
        'peptideSequence': 'AAAAAAAAAPAAAATAPTTAATTAATAAQ',
        'fragments': {
            'fragmentMZ': [201.08698363000033, 214.118618035, 285.155732035, 356.192846035, 427.229960035, 498.26707403499995, 569.3041880349999, 640.3413020349999, 737.3940660349999, 808.4311800349999, 879.4682940349999, 950.5054080349998, 1021.5425220349998, 1122.5902010349996, 1193.6273150349996, 218.11353273500032, 289.1506467350005, 390.19832573500054, 461.2354397350005, 532.2725537350005, 633.3202327350006, 734.3679117350006, 805.4050257350005, 876.4421397350005, 977.4898187350005, 1078.5374977350004, 1175.5902617350005, 1246.6273757350004, 1728.8762747350002, 1799.9133887350004, 1104.5796363349996, 1175.6167503349996, 1060.5269330350004, 1157.5796970350004, 1228.6168110350004, 1710.8657100350001, 696.3675170349999, 868.4523100349999, 954.4947065349999, 816.4153933850001, 807.4101110350001],
            'fragmentNumber': [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 20, 21, 15, 16, 12, 13, 14, 20, 18, 22, 24, 19, 19],
            'fragmentLossType': ['NH3', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'noloss', 'noloss', 'noloss', 'noloss', 'H2O'],
            'fragmentType': ['y', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'b', 'b', 'y', 'y', 'y', 'y', 'b', 'b', 'b', 'y', 'y'],
            'fragmentAnnotation': ['y2-NH3^+1', 'b3^+1', 'b4^+1', 'b5^+1', 'b6^+1', 'b7^+1', 'b8^+1', 'b9^+1', 'b10^+1', 'b11^+1', 'b12^+1', 'b13^+1', 'b14^+1', 'b15^+1', 'b16^+1', 'y2^+1', 'y3^+1', 'y4^+1', 'y5^+1', 'y6^+1', 'y7^+1', 'y8^+1', 'y9^+1', 'y10^+1', 'y11^+1', 'y12^+1', 'y13^+1', 'y14^+1', 'y20^+1', 'y21^+1', 'b15-H2O^+1', 'b16-H2O^+1', 'y12-H2O^+1', 'y13-H2O^+1', 'y14-H2O^+1', 'y20-H2O^+1', 'b18^+2', 'b22^+2', 'b24^+2', 'y19^+2', 'y19-H2O^+2'],
            'fragmentIntensity': [22774.28, 157563.6, 444095.9, 993978.3, 1666265, 2408575, 3411004, 2116178, 258959.3, 298698, 278765, 255920.5, 139878.8, 119512, 82443.61, 430567.5, 52544.59, 183384.2, 105601.1, 72300.24, 56770.18, 70014.88, 18832.65, 17846.54, 113519.6, 118261.5, 488446.9, 97438.05, 273214.5, 16548.61, 117036.2, 488446.9, 25312.04, 134764.2, 24459.05, 77885.37, 107844.3, 19617.25, 18051.07, 18561.47, 18938.11],
            'fragmentCharge': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
        }
    }

    print(converter.assay_to_dataframe(assay0))