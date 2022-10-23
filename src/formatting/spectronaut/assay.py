from .mod import parse_modification_Spectronaut, \
    to_modified_sequence_Spectronaut
 

def Spectronaut_assay_library_columns():
    
    def modified_sequence(assay, **kwargs):
        return to_modified_sequence_Spectronaut(
            sequence=assay['peptideSequence'],
            modification=assay.get('modification', None)
        )

    columns = [
        {
            'name': 'PrecursorMz',
            'path': 'precursorMZ'
        },
        {
            'name': 'IonMobility',
            'path': 'ionMobility'
        },
        {
            'name': 'FragmentMz',
            'path': ['fragments', 'fragmentMZ']
        },
        {
            'name': 'iRT',
            'path': 'iRT'
        },
        {
            'name': 'RelativeFragmentIntensity',
            'path': ['fragments', 'fragmentIntensity']
        },
        {
            'name': 'StrippedSequence',
            'path': 'peptideSequence'
        },
        {
            'name': 'PrecursorCharge',
            'path': 'precursorCharge',
        },
        {
            'name': 'FragmentType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'name': 'FragmentNumber',
            'path': ['fragments', 'fragmentNumber']
        },
        {
            'name': 'FragmentCharge',
            'path': ['fragments', 'fragmentCharge']
        },
        {
            'name': 'FragmentLossType',
            'path': ['fragments', 'fragmentLossType']
        },
        {
            'name': 'ModifiedPeptide',
            'function': modified_sequence,
        },
        {
            'name': 'ProteinId',
            'path': ['metadata', 'protein']
        }
    ]
    return columns


def Spectronaut_library_parsing_columns():
    columns = [
        {
            'key': True,
            'name': 'StrippedPeptide',
            'path': 'peptideSequence'
        },
        {
            'key': True,
            'name': 'ModifiedPeptide',
            'path': 'modification',
            'convert': parse_modification_Spectronaut
        },
        {
            'key': True,
            'name': 'PrecursorCharge',
            'path': 'precursorCharge',
            'convert': int
        },
        {
            'key': True,
            'name': 'ProteinGroups',
            'path': ['metadata', 'protein']
        },
        {
            'key': True,
            'name': 'ReferenceRun',
            'path': ['metadata', 'file']
        },
        {
            'list': True,
            'name': 'FragmentMz',
            'path': ['fragments', 'fragmentMZ']
        },
        {
            'list': True,
            'name': 'RelativeIntensity',
            'path': ['fragments', 'fragmentIntensity']
        },
        {
            'list': True,
            'name': 'FragmentCharge',
            'path': ['fragments', 'fragmentCharge']
        },
        {
            'list': True,
            'name': 'FragmentType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'list': True,
            'name': 'FragmentNumber',
            'path': ['fragments', 'fragmentNumber']
        },
        {
            'list': True,
            'name': 'FragmentLossType',
            'path': ['fragments', 'fragmentLossType']
        },
        {
            'name': 'iRT',
            'path': 'iRT'
        },
        {
            'name': 'PrecursorMz',
            'path': 'precursorMZ'
        },
        {
            'name': 'IonMobility',
            'path': 'ionMobility'
        },
        {
            'name': 'ReferenceRunQvalue',
            'path': ['metadata', 'qvalue']
        }
    ]

    return columns


def SpectroMine_assay_parsing_columns():
    columns = [
        {
            'key': True,
            'name': 'PEP.StrippedSequence',
            'path': 'peptideSequence'
        },
        {
            'key': True,
            'name': 'PP.PIMID',
            'path': 'modification',
            'convert': parse_modification_Spectronaut
        },
        {
            'key': True,
            'name': 'PP.Charge',
            'path': 'precursorCharge',
            'convert': int
        },
        {
            'key': True,
            'name': 'PG.ProteinAccessions',
            'path': ['metadata', 'protein']
        },
        {
            'key': True,
            'name': 'R.FileName',
            'path': ['metadata', 'file']
        },
        {
            'key': True,
            'name': 'PSM.MS2ScanNumber',
            'path': ['metadata', 'scan']
        },
        {
            'list': True,
            'name': 'FI.CalibratedMZ',
            'path': ['fragments', 'fragmentMZ']
        },
        {
            'list': True,
            'name': 'FI.Intensity',
            'path': ['fragments', 'fragmentIntensity']
        },
        {
            'list': True,
            'name': 'FI.Charge',
            'path': ['fragments', 'fragmentCharge']
        },
        {
            'list': True,
            'name': 'FI.FrgType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'list': True,
            'name': 'FI.FrgNum',
            'path': ['fragments', 'fragmentNumber']
        },
        {
            'list': True,
            'name': 'FI.LossType',
            'path': ['fragments', 'fragmentLossType']
        },
        {
            'name': 'PP.EmpiricalRT',
            'path': 'rt'
        },
        {
            'name': 'PP.iRTEmpirical',
            'path': 'iRT'
        },
        {
            'name': 'PSM.CalibratedMS1MZ',
            'path': 'precursorMZ'
        },
        {
            'name': 'PSM.IonMobility',
            'path': 'ionMobility'
        },
        {
            'name': 'PSM.Qvalue',
            'path': ['metadata', 'qvalue']
        },
        {
            'name': 'PSM.Score',
            'path': ['metadata', 'score']
        }
    ]

    return columns

