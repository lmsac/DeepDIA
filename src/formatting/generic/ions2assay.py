from assay import AssayBuilder
from .mod import parse_modification


class IonsToAssayConverter():
    def __init__(self, assay_builder=None,
                 ions_labels=None,
                 skip_zeros=True,
                 **kwargs):
        if assay_builder is None:
            assay_builder = AssayBuilder()
        self.assay_builder = assay_builder

        self.skip_zeros = skip_zeros

        if ions_labels is None:
            ions_labels = {
                'loss': {
                    'NH3': 'n',
                    'H2O': 'o',
                    'H3PO4': 'p'
                }
            }
        self.ions_labels = ions_labels


    def get_ions_label(self, fragment_type, charge, loss=None):
        def _get_ions_label(key, value):
            d = self.ions_labels.get(key, None)
            if d is None:
                return value
            x = d.get(value, None)
            if x is not None:
                return x
            else:
                return value

        frgtype = _get_ions_label('fragment_type', fragment_type)
        frgloss = _get_ions_label('loss', loss) \
            if loss is not None and loss != 'noloss' else ''
        frgcharge = str(charge)

        return frgtype + frgloss + frgcharge


    def theoretical_fragments(self, ions):
        sequence = ions['peptide']
        modification = parse_modification(
            ions.get('modification', None)
        )

        return self.assay_builder.theoretical_fragments(
            sequence=sequence,
            modification=modification
        )


    def get_intensity(self, ions,
                      fragment_type, fragment_number, charge, loss=None):
        label = self.get_ions_label(
            fragment_type=fragment_type,
            loss=loss,
            charge=charge
        )

        frgintensity = ions['ions'].get(label, None)
        if isinstance(frgintensity, list):
            return frgintensity[fragment_number - 1]
        else:
            None


    def ions_entry_to_assay(self, ions):
        result = self.theoretical_fragments(
            ions
        )

        precursor_charge = ions.get('charge', None)
        if precursor_charge is not None:
            precursor_charge = int(precursor_charge)
            result['precursorCharge'] = precursor_charge
            result = self.assay_builder.update_precursor_mz(result)

        fragment_intensity = []
        for i, _ in enumerate(result['fragments']['fragmentMZ']):
            intensity = self.get_intensity(
                ions,
                fragment_type=result['fragments']['fragmentType'][i],
                fragment_number=result['fragments']['fragmentNumber'][i],
                loss=result['fragments']['fragmentLossType'][i],
                charge=result['fragments']['fragmentCharge'][i]
            )
            if intensity is not None:
                fragment_intensity.append(intensity)
            else:
                fragment_intensity.append(0)

        result['fragments']['fragmentIntensity'] = fragment_intensity

        if self.skip_zeros:
            from sys import float_info
            result = self.assay_builder.filter_fragments_by_intensity(
                result, absolute_intensity=float_info.min, copy=False
            )

        metadata = ions.get('metadata', None)
        if isinstance(metadata, dict):
            result['metadata'] = metadata.copy()

        return result


    def ions_to_assays(self, ions,
                       return_generator=False):
        result = (
            self.ions_entry_to_assay(entry)
            for entry in ions
        )
        if not return_generator:
            result = list(result)
        return result




if __name__ == '__main__':
    converter = IonsToAssayConverter()

    s0 = '''{
  "peptide": "AAAAAAAAAPAAAATAPTTAATTAATAAQ",
  "charge": 2,
  "ions":{
    "b1":[0, 0, 157563.6, 444095.9, 993978.3, 1666265, 2408575, 3411004, 2116178, 258959.3, 298698, 278765, 255920.5, 139878.8, 119512, 82443.61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "bn1":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "bo1":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 117036.2, 488446.9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "b2":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 107844.3, 0, 0, 0, 19617.25, 0, 18051.07, 0, 0, 0, 0],
    "bn2":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "bo2":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "y1":[0, 430567.5, 52544.59, 183384.2, 105601.1, 72300.24, 56770.18, 70014.88, 18832.65, 17846.54, 113519.6, 118261.5, 488446.9, 97438.05, 0, 0, 0, 0, 0, 273214.5, 16548.61, 0, 0, 0, 0, 0, 0, 0],
    "yn1":[0, 22774.28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "yo1":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25312.04, 134764.2, 24459.05, 0, 0, 0, 0, 0, 77885.37, 0, 0, 0, 0, 0, 0, 0, 0],
    "y2":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18561.47, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "yn2":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "yo2":[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18938.11, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  },
  "qvalue": 0
}'''

    s1 = '''{
  "peptide": "ESKSSPRPTAEK",
  "ions": {
    "b2": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.11816668510437012],
    "y2": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3882259130477905, 0.45884573459625244],
    "b1": [0.09649229049682617, 0.0, 0.007291555404663086, 0.09996521472930908, 0.04090571403503418, 0.0, 0.07087504863739014, 0.0, 0.0, 0.03869640827178955, 0.24884772300720215],
    "y1": [0.48307085037231445, 0.0, 0.0, 0.0, 0.19371294975280762, 0.07264506816864014, 0.5678353309631348, 0.4612734317779541, 1.0, 0.12032663822174072, 0.009947061538696289]
  },
  "modification": "S2(Phospho)",
  "charge": 2
}'''

    import json
    print(converter.ions_entry_to_assay(ions=json.loads(s0)))
