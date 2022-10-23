from .mod import to_modified_sequence_MaxQuant
from pepmass.modinfo import default_fixed_modifications, \
    default_variable_modifications

def MaxQuant_assay_library_evidence_columns(
    fixed_modifications=default_fixed_modifications(),
    variable_modifications=default_variable_modifications()):

    def modified_sequence(assay, **kwargs):
        return to_modified_sequence_MaxQuant(
            sequence=assay['peptideSequence'],
            modification=assay.get('modification', None),
            fixed_modifications=fixed_modifications,
            variable_modifications=variable_modifications
        )

    def get_id(assay, **kwargs):
        return kwargs.get('index')

    columns = [
        {
            'name': 'Sequence',
            'path': 'peptideSequence'
        },
        {
            'name': 'Modified sequence',
            'function': modified_sequence,
        },
        {
            'name': 'Proteins',
            'path': ['metadata', 'protein']
        },
        {
            'name': 'MS/MS IDs',
            'function': get_id
        },
        {
            'name': 'Raw file',
            'path': ['metadata', 'file']
        },
        {
            'name': 'Type',
            'path': ['metadata', 'msmsType']
        },
        {
            'name': 'Reverse',
            'path': ['metadata', 'decoy'],
            'convert': lambda x: '+' if x else ''
        },
        {
            'name': 'Charge',
            'path': 'precursorCharge',
        },
        {
            'name': 'm/z',
            'path': 'precursorMZ'
        },
        {
            'name': 'Calibrated retention time',
            'path': 'iRT'
        },
        {
            'name': 'Calibrated retention time start',
            'path': 'iRT',
            'convert': lambda x: x - 0.3
        },
        {
            'name': 'Calibrated retention time finish',
            'path': 'iRT',
            'convert': lambda x: x + 0.3
        },
        {
            'name': 'Retention time',
            'path': 'iRT'
        },
        {
            'name': '1/K0',
            'path': 'ionMobility'
        },
        {
            'name': 'Calibrated 1/K0',
            'path': 'ionMobility'
        },
        {
            'name': 'Intensity',
            'path': ['metadata', 'ms1Intensity'],
            'default': 10000
        }
    ]

    return columns


def MaxQuant_assay_library_msms_columns():

    def fragment_name(assay, **kwargs):
        fragments = assay['fragments']
        fragment_type = fragments['fragmentType']
        fragment_number = fragments['fragmentNumber']

        result = [
            str(t) + str(n) for t, n in zip(fragment_type, fragment_number)
        ]

        fragment_losstype = fragments.get('fragmentLossType', None)
        if fragment_losstype is not None:
            for i, l in enumerate(fragment_losstype):
                if l and l != 'noloss' and l != 'None':
                    result[i] += '-' + str(l)

        fragment_charge = fragments.get('fragmentCharge', None)
        if fragment_charge is not None:
            for i, l in enumerate(fragment_charge):
                if l and l != 1:
                    result[i] += '^' + str(l) + '+'

        return ';'.join(result)

    def get_id(assay, **kwargs):
        return kwargs.get('index')

    columns = [
        {
            'name': 'Fragmentation',
            'path': ['metadata', 'fragmentation']
        },
        {
            'name': 'Mass analyzer',
            'path': ['metadata', 'massAnalyzer']
        },
        {
            'name': 'Retention time',
            'path': 'iRT'
        },
        {
            'name': 'Retention time',
            'path': 'iRT'
        },
        {
            'name': 'PEP',
            'path': ['metadata', 'pep'],
            'default': 0.0
        },
        {
            'name': 'Score',
            'path': ['metadata', 'score'],
            'default': 150.0
        },
        {
            'name': 'Matches',
            'function': fragment_name
        },
        {
            'name': 'Intensites',
            'path': ['fragments', 'fragmentIntensity'],
            'convert': lambda x: ';'.join(map(str, x))
        },
        {
            'name': 'Masses',
            'path': ['fragments', 'fragmentMZ'],
            'convert': lambda x: ';'.join(map(str, x))
        },
        {
            'name': 'id',
            'function': get_id
        }
    ]

    return columns


def MaxQuant_assay_library_peptide_columns():
    columns = [
        {
            'name': 'Sequence',
            'path': 'peptideSequence'
        },
        {
            'name': 'Unique (Proteins)',
            'path': ['metadata', 'protein'],
            'convert': lambda x: 'no' if len(x.split(';')) > 1 else 'yes'
        },
        {
            'name': 'Leading razor protein',
            'path': ['metadata', 'protein']
        },
        {
            'name': 'Start position',
            'path': ['metadata', 'proteinStartPosition']
        },
        {
            'name': 'Missed cleavages',
            'path': ['metadata', 'missedCleavages']
        }
    ]

    return columns


