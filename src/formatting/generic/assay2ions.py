from .mod import stringify_modification


class AssayToIonsConverter():
    def __init__(self,
                 fragment_types=None,
                 fragment_charges=None,
                 fragment_loss_types=None,
                 ions_labels=None,
                 **kwargs):
        if fragment_types is not None:
            self.fragment_types = fragment_types
        else:
            self.fragment_types = ['b', 'y']

        if fragment_charges is not None:
            self.fragment_charges = fragment_charges
        else:
            self.fragment_charges = [1, 2]

        if fragment_loss_types is not None:
            self.fragment_loss_types = fragment_loss_types
        else:
            self.fragment_loss_types = ['noloss', 'NH3', 'H2O']

        if ions_labels is not None:
            self.ions_labels = ions_labels
        else:
            self.ions_labels = {
                'loss': {
                    'NH3': 'n',
                    'H2O': 'o',
                    'H3PO4': 'p'
                }
            }


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


    def assay_to_ions_entry(self, assay):
        peptide = assay['peptideSequence']
        ions = {}
        fragments = assay['fragments']
        for i, x in enumerate(fragments['fragmentIntensity']):
            frgtype = fragments['fragmentType'][i]
            frgcharge = fragments['fragmentCharge'][i]
            frgloss = fragments['fragmentLossType'][i]
            if frgloss is None or frgloss == '':
                frgloss = 'noloss'

            if frgtype not in self.fragment_types or \
                frgcharge not in self.fragment_charges or \
                frgloss not in self.fragment_loss_types:
                continue

            label = self.get_ions_label(
                fragment_type=frgtype,
                charge=frgcharge,
                loss=frgloss
            )

            frgnumber = int(fragments['fragmentNumber'][i])

            ions.setdefault(label, [0] * (len(peptide) - 1)) \
                [frgnumber - 1] = x

        modification = stringify_modification(
            assay.get('modification', None)
        )

        result = {
            'peptide': peptide,
            'charge': assay.get('precursorCharge', None),
            'modification': modification,
            'ions': ions,
        }

        metadata = assay.get('metadata', None)
        if isinstance(metadata, dict):
            result['metadata'] = metadata.copy()

        return result


    def assays_to_ions(self, assays,
                            return_generator=False):
        result = (
            self.assay_to_ions_entry(assay)
            for assay in assays
        )
        if not return_generator:
            result = list(result)
        return result



if __name__ == '__main__':
    converter = AssayToIonsConverter()

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

    print(converter.assay_to_ions_entry(assay0))