from .mod import parse_modification_MaxQuant


def MaxQuant_assay_parsing_columns():
    import re

    columns = [
        {
            'key': True,
            'name': 'Sequence',
            'path': 'peptideSequence'
        },
        {
            'key': True,
            'name': 'Modified sequence',
            'path': 'modification',
            'convert': parse_modification_MaxQuant
        },
        {
            'key': True,
            'name': 'Charge',
            'path': 'precursorCharge',
            'convert': int
        },
        {
            'key': True,
            'name': 'Proteins',
            'path': ['metadata', 'protein']
        },
        {
            'key': True,
            'name': 'Raw file',
            'path': ['metadata', 'file']
        },
        {
            'key': True,
            'name': 'Scan number',
            'path': ['metadata', 'scan']
        },
        {
            'name': 'Masses',
            'path': ['fragments', 'fragmentMZ'],
            'convert': lambda s: \
                [] if s is None or str(s) == 'nan' else \
                [s] if isinstance(s, float) else \
                [float(x) for x in str(s).split(';')]
        },
        {
            'name': 'Intensities',
            'path': ['fragments', 'fragmentIntensity'],
            'convert': lambda s: \
                [] if s is None or str(s) == 'nan' else \
                [s] if isinstance(s, float) else \
                [float(x) for x in s.split(';')]
        },
        {
            'name': 'Matches',
            'path': ['fragments', 'fragmentCharge'],
            'convert': lambda s: \
                [] if s is None or str(s) == 'nan' else \
                [
                    int(next(iter(re.findall('\\(([0-9])+\\+\\)', x)), 1)) 
                    for x in s.split(';')
                ]
        },
        {
            'name': 'Matches',
            'path': ['fragments', 'fragmentType'],
            'convert': lambda s: \
                [] if s is None or str(s) == 'nan' else \
                [
                    next(iter(re.findall('^([a-zA-Z]+)', x)), None)
                    for x in s.split(';')
                ]
        },
        {
            'name': 'Matches',
            'path': ['fragments', 'fragmentNumber'],
            'convert': lambda s: \
                [] if s is None or str(s) == 'nan' else \
                [
                    int(next(iter(re.findall('^[a-zA-Z]+([0-9]+)', x)), -1))
                    for x in s.split(';')
                ]
        },
        {
            'name': 'Matches',
            'path': ['fragments', 'fragmentLossType'],
            'convert': lambda s: \
                [] if s is None or str(s) == 'nan' else \
                    [
                    next(iter(re.findall('^[a-zA-Z]+[0-9]+-([a-zA-Z0-9]+)', x)), None)
                    for x in s.split(';')
                ]
        },
        {
            'name': 'Retention time',
            'path': 'rt'
        },
        {
            'name': 'm/z',
            'path': 'precursorMZ'
        },
        {
            'name': 'Score',
            'path': ['metadata', 'score']
        }
    ]

    return columns