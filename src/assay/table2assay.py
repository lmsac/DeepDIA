import numpy as np

from .modseq import parse_modification
from .values import set_assay_values


class DataFrameToAssayConverter:
    def __init__(self, columns=None):
        if columns is None:
            columns = default_assay_parsing_columns()
        self.columns = columns


    def parse_psm(self, data):
        values = {}
        for col in self.columns:
            if not isinstance(col, dict):
                if col in data.columns:
                    values[col] = data[col].iloc[0]
                continue

            name = col.get('name')
            if name in data.columns:
                if col.get('list', False):
                    values[name] = data[name].tolist()
                else:
                    values[name] = data[name].iloc[0]

        for name, value in list(values.items()):
            if isinstance(value, np.generic):
                 values[name] = value.item()

        assay = {}
        set_assay_values(assay, self.columns, **values)

        return assay


    def dataframe_to_assays(self, data, return_generator=False):
        key = data.columns.intersection([
            c.get('name')
            for c in self.columns
            if isinstance(c, dict) and c.get('key', False)
        ])
        assay_id = data[key].ne(data[key].shift()).any(axis=1).cumsum()
        groups = data.groupby(by=assay_id, sort=False, as_index=False)

        results = (
            self.parse_psm(group)
            for i, group in groups
        )
        if not return_generator:
            results = list(results)
        return results


def default_assay_parsing_columns():
    columns = [
        {
            'key': True,
            'name': 'peptideSequence',
            'path': 'peptideSequence'
        },
        {
            'key': True,
            'name': 'modification',
            'path': 'modification',
            'convert': parse_modification
        },
        {
            'key': True,
            'name': 'precursorCharge',
            'path': 'precursorCharge',
            'convert': int
        },
        {
            'key': True,
            'name': 'protein',
            'path': ['metadata', 'protein']
        },
        {
            'key': True,
            'name': 'file',
            'path': ['metadata', 'file']
        },
        {
            'key': True,
            'name': 'scan',
            'path': ['metadata', 'scan'],
            'convert': int
        },
        {
            'list': True,
            'name': 'fragmentMZ',
            'path': ['fragments', 'fragmentMZ']
        },
        {
            'list': True,
            'name': 'fragmentIntensity',
            'path': ['fragments', 'fragmentIntensity']
        },
        {
            'list': True,
            'name': 'fragmentCharge',
            'path': ['fragments', 'fragmentCharge']
        },
        {
            'list': True,
            'name': 'fragmentType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'list': True,
            'name': 'fragmentNumber',
            'path': ['fragments', 'fragmentNumber']
        },
        {
            'list': True,
            'name': 'fragmentLossType',
            'path': ['fragments', 'fragmentLossType']
        },
        {
            'name': 'rt',
            'path': 'rt',
            'convert': float
        },
        {
            'name': 'iRT',
            'path': 'iRT',
            'convert': float
        },
        {
            'name': 'precursorMZ',
            'path': 'precursorMZ',
            'convert': float
        },
        {
            'name': 'ionMobility',
            'path': 'ionMobility',
            'convert': float
        },
        {
            'name': 'qvalue',
            'path': ['metadata', 'qvalue'],
            'convert': float
        },
        {
            'name': 'score',
            'path': ['metadata', 'score'],
            'convert': float
        }
    ]

    return columns