from copy import deepcopy
import numpy as np
import itertools

from pepmass import ModifiedPeptideMassCalculator
from .modseq import check_modifications

class AssayBuilder():
    def __init__(self,
                 mass_calculator=None,
                 fragment_types=None,
                 fragment_charges=None,
                 fragment_loss_types=None,
                 **kwargs):
        if mass_calculator is None:
            mass_calculator = ModifiedPeptideMassCalculator()
        self.mass_calculator = mass_calculator

        if fragment_types is None:
            fragment_types = ['b', 'y']
        self.fragment_types = fragment_types

        if fragment_charges is None:
            fragment_charges = [1, 2]
        self.fragment_charges = fragment_charges

        if fragment_loss_types is None:
            fragment_loss_types = ['noloss', 'NH3', 'H2O']
        self.fragment_loss_types = fragment_loss_types


    def assay(self, sequence, charge=None, modification=None, fragments=None,
              **kwargs):
        result = {
            'peptideSequence': sequence,
            'modification': modification
        }
        if charge is not None:
            result['precursorCharge'] = charge

        if fragments is None:
            fragments = {
                'fragmentMZ': [],
                'fragmentAnnotation': [],
                #'fragmentIntensity': [],
                'fragmentType': [],
                'fragmentNumber': [],
                'fragmentCharge': [],
                'fragmentLossType': []
            }
        result.update({
            'fragments': fragments
        })
        result.update(kwargs)
        return result


    def theoretical_fragments(self, sequence, modification=None, **kwargs):
        def get_fragment_annotation(fragment_type, fragment_number,
                                    charge, loss=None, **kwargs):
            return fragment_type + \
                (str(fragment_number) if fragment_number is not None else '') + \
                ('-' + loss if loss is not None and loss != 'noloss' else '') + \
                '^+' + str(charge)

        fragment_mz = []
        fragment_type = []
        fragment_number = []
        fragment_charge = []
        fragment_loss_type = []
        fragment_annotation = []

        frag_type = kwargs.pop('fragment_type', None)
        if frag_type is None:
            frag_type = list(self.fragment_types)
        else:
            if isinstance(frag_type, str):
                frag_type = [frag_type]
            frag_type = [
                x for x in frag_type
                if x in self.fragment_types
            ]

        frag_loss_type = kwargs.pop('fragment_loss_type', None)
        if frag_loss_type is None:
            frag_loss_type = list(self.fragment_loss_types)
        else:
            if isinstance(frag_loss_type, str):
                frag_loss_type = [frag_loss_type]
            frag_loss_type = [
                x for x in frag_loss_type
                if x in self.fragment_loss_types
            ]

        frag_charge = kwargs.pop('fragment_charge', None)
        if frag_charge is None:
            frag_charge = list(self.fragment_charges)
        else:
            if isinstance(frag_charge, int):
                frag_charge = [frag_charge]
            frag_charge = [
                x for x in fragment_charge
                if x in self.fragment_charges
            ]

        peptide_fragments = self.mass_calculator.fragment_mz(
            sequence=sequence,
            modification=modification,
            fragment_type=frag_type,
            loss=frag_loss_type,
            charge=frag_charge,
            **kwargs
        )

        for x in peptide_fragments:
            frgtype = x['fragment_type']
            frglossTpye = x['loss']
            frgcharge = x['charge']
            frgmz = x['fragment_mz']
            frgnum = [i + 1 for i, x in enumerate(frgmz) if x is not None]
            frgmz = [frgmz[i - 1] for i in frgnum]
            frgannot = [
                get_fragment_annotation(
                    fragment_type=frgtype,
                    fragment_number=i,
                    loss=frglossTpye,
                    charge=frgcharge
                )
                for i in frgnum
            ]

            fragment_mz.extend(frgmz)
            fragment_number.extend(frgnum)
            fragment_type.extend([frgtype] * len(frgnum))
            fragment_charge.extend([frgcharge] * len(frgnum))
            fragment_loss_type.extend([frglossTpye] * len(frgnum))
            fragment_annotation.extend(frgannot)

        result = {
            'peptideSequence': sequence,
            'modification': modification,

            'fragments': {
                'fragmentMZ': fragment_mz,
                'fragmentType': fragment_type,
                'fragmentNumber': fragment_number,
                'fragmentCharge': fragment_charge,
                'fragmentLossType': fragment_loss_type,
                'fragmentAnnotation': fragment_annotation
            }
        }

        return result


    def update_precursor_mz(self, assay):
        sequence = assay['peptideSequence']
        modification = assay.get('modification', None)
        precursor_charge = int(assay['precursorCharge'])

        precursor_mz = self.mass_calculator.precursor_mz(
            sequence=sequence,
            modification=modification,
            charge=precursor_charge
        )

        assay.update({
            'precursorMZ': precursor_mz
        })
        return assay


    def update_fragment_mz(self, assay, **kwargs):
        sequence = assay['peptideSequence']
        modification = assay.get('modification', None)

        fragments = deepcopy(assay['fragments'])
        fragment_type = fragments['fragmentType']
        fragment_number = fragments['fragmentNumber']
        fragment_charge = fragments['fragmentCharge']
        fragment_loss_type = fragments['fragmentLossType']

        fragment_types = list(set(fragment_type) \
            .intersection(self.fragment_types))
        if len(fragment_types) > 0:
            peptide_fragments = self.mass_calculator.fragment_mz(
                sequence=sequence,
                modification=modification,
                fragment_type=fragment_types,
                loss=list(set(fragment_loss_type)),
                charge=list(set(fragment_charge)),
                **kwargs
            )
        else:
            peptide_fragments = []

        fragment_mz = []
        for i, _ in enumerate(fragment_type):
            if fragment_type[i] not in self.fragment_types:
                raise ValueError('fragment not found: ' + fragment_type[i])

            mz = None
            for x in peptide_fragments:
                if x['fragment_type'] == fragment_type[i] and \
                    x['charge'] == fragment_charge[i] and \
                    (x['loss'] == fragment_loss_type[i] or \
                     x['loss'] == 'noloss' and fragment_loss_type[i] is None):
                    mz = x['fragment_mz'][fragment_number[i] - 1]
                    fragment_mz.append(mz)
                    if mz is None:
                        mz = 0
                    break
            if mz is None:
                raise ValueError('fragment not found: ' + fragment_type[i])

        assay['fragments'].update({
            'fragmentMZ': fragment_mz
        })
        return assay


    def filter_fragments_by_index(self, assay, fragment_index,
                                  invert=False, copy=True):
        if copy:
            assay = deepcopy(assay)

        if invert:
            fragment_index = set(fragment_index)
            for k, v in assay['fragments'].items():
                assay['fragments'][k] = [
                    x for i, x in enumerate(v)
                    if i not in fragment_index
                ]
        else:
            for k, v in assay['fragments'].items():
                assay['fragments'][k] = [v[i] for i in fragment_index]

        return assay


    def filter_fragments_by_type(self, assay, fragment_type=None,
                                 return_index=False, copy=True):
        if fragment_type == None:
            fragment_type = self.fragment_types

        fragment_index = [
            i
            for i, x in enumerate(assay['fragments']['fragmentType'])
            if x is not None and (\
            isinstance(fragment_type, str) and x == fragment_type or \
            x in fragment_type)
        ]

        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=fragment_index,
                copy=copy
            )

    def filter_fragments_by_amino_acid_number(self, assay,
                                              min_amino_acid_number,
                                              return_index=False,
                                              copy=True):
        fragment_index = [
            i
            for i, x in enumerate(assay['fragments']['fragmentNumber'])
            if x is not None and x >= min_amino_acid_number
        ]

        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=fragment_index,
                copy=copy
            )


    def filter_fragments_by_charge(self, assay, fragment_charge=None,
                                   return_index=False, copy=True):
        if fragment_charge == None:
            fragment_charge = self.fragment_charges

        fragment_index = [
            i
            for i, x in enumerate(assay['fragments']['fragmentCharge'])
            if x is not None and (\
            isinstance(fragment_charge, int) and x == fragment_charge or \
            x in fragment_charge)
        ]

        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=fragment_index,
                copy=copy
            )


    def filter_fragments_by_loss_type(self, assay, fragment_loss_type=None,
                                      return_index=False, copy=True):
        if fragment_loss_type == None:
            fragment_loss_type = self.fragment_loss_types

        fragment_index = [
            i
            for i, x in enumerate(assay['fragments']['fragmentLossType'])
            if (x is None or x == 'None' or x == '') and \
                (fragment_loss_type == 'noloss' or 'noloss' in fragment_loss_type) or \
               (x == fragment_loss_type or x in fragment_loss_type)
        ]

        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=fragment_index,
                copy=copy
            )


    def filter_fragments_by_mz(self, assay, min_mz=None, max_mz=None,
                               return_index=False, copy=True):
        fragment_index = [
            i
            for i, x in enumerate(assay['fragments']['fragmentMZ'])
            if x is not None and \
            (min_mz is None or x >= min_mz) and \
            (max_mz is None or x <= max_mz)
        ]

        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=fragment_index,
                copy=copy
            )


    def exclude_fragments_in_isolation_window(
            self, assay, swath_windows, return_index=False, copy=True):
        precursor_mz = assay['precursorMZ']
        isolation_window_index = np.where(
            (swath_windows['start'] < precursor_mz) & \
            (swath_windows['end'] > precursor_mz)
        )[0]
        if len(isolation_window_index) == 0:
            if not return_index:
                return assay
            else:
                return [
                    i for i, x in enumerate(assay['fragments']['fragmentMZ'])
                ]

        exclude = set(itertools.chain.from_iterable((
            self.filter_fragments_by_mz(
                assay,
                min_mz=swath_windows['start'][i],
                max_mz=swath_windows['end'][i],
                return_index=True
            )
            for i in isolation_window_index
        )))

        if not return_index:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=list(exclude),
                invert=True, copy=copy
            )
        else:
            return [
                i for i, x in enumerate(assay['fragments']['fragmentMZ'])
                if i not in exclude
            ]


    def filter_fragments_by_intensity(self, assay,
                                      absolute_intensity=None,
                                      relative_intensity=None,
                                      top_n=None,
                                      return_index=False, copy=True):
        if len(assay['fragments']['fragmentIntensity']) > 0:
            if relative_intensity is not None:
                intensity = relative_intensity * max(
                    x for x in assay['fragments']['fragmentIntensity']
                    if x is not None
                )
                if absolute_intensity is None:
                    absolute_intensity = intensity
                else:
                    absolute_intensity = max(absolute_intensity, intensity)

            fragment_intensity = assay['fragments']['fragmentIntensity']
            fragment_index = [
                i
                for i, x in enumerate(fragment_intensity)
                if x is not None and \
                (absolute_intensity is None or x >= absolute_intensity)
            ]

            if top_n is not None:
                fragment_index = sorted(
                    fragment_index,
                    key=lambda k: fragment_intensity[k],
                    reverse=True
                )[:top_n]
                fragment_index.sort()

        else:
            fragment_index = []

        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay,
                fragment_index=fragment_index,
                copy=copy
            )


    def filter_fragments(self, assay,
                         max_fragment_number=None,
                         fragment_type=None,
                         fragment_charge=None,
                         fragment_loss_type=None,
                         min_fragment_amino_acid_number=None,
                         min_fragment_mz=None,
                         max_fragment_mz=None,
                         swath_windows=None,
                         min_relative_fragment_intensity=None,
                         return_index=False, copy=True):
        fragment_index = None
        if fragment_type is not None:
            fragment_index_1 = self.filter_fragments_by_type(
                assay, fragment_type,
                return_index=True
            )
            if fragment_index is None:
                fragment_index = set(fragment_index_1)
            else:
                fragment_index = fragment_index \
                    .intersection(fragment_index_1)

        if fragment_charge is not None:
            fragment_index_1 = self.filter_fragments_by_charge(
                assay, fragment_charge,
                return_index=True
            )
            if fragment_index is None:
                fragment_index = set(fragment_index_1)
            else:
                fragment_index = fragment_index \
                    .intersection(fragment_index_1)

        if fragment_loss_type is not None:
            fragment_index_1 = self.filter_fragments_by_loss_type(
                assay, fragment_loss_type,
                return_index=True
            )
            if fragment_index is None:
                fragment_index = set(fragment_index_1)
            else:
                fragment_index = fragment_index \
                    .intersection(fragment_index_1)

        if min_fragment_amino_acid_number is not None:
            fragment_index_1 = self.filter_fragments_by_amino_acid_number(
                assay,
                min_amino_acid_number=min_fragment_amino_acid_number,
                return_index=True
            )
            if fragment_index is None:
                fragment_index = set(fragment_index_1)
            else:
                fragment_index = fragment_index \
                    .intersection(fragment_index_1)

        if min_fragment_mz is not None or \
            max_fragment_mz is not None:
            fragment_index_1 = self.filter_fragments_by_mz(
                assay, min_mz=min_fragment_mz, max_mz=max_fragment_mz,
                return_index=True
            )
            if fragment_index is None:
                fragment_index = set(fragment_index_1)
            else:
                fragment_index = fragment_index \
                    .intersection(fragment_index_1)

        if swath_windows is not None:
            fragment_index_1 = self.exclude_fragments_in_isolation_window(
                assay, swath_windows=swath_windows,
                return_index=True
            )
            if fragment_index is None:
                fragment_index = set(fragment_index_1)
            else:
                fragment_index = fragment_index \
                    .intersection(fragment_index_1)

        if fragment_index is not None:
            fragment_index = list(fragment_index)
            assay_1 = self.filter_fragments_by_index(assay, fragment_index)
        else:
            assay_1 = assay

        if max_fragment_number is not None or \
            min_relative_fragment_intensity is not None:
            fragment_index_1 = self.filter_fragments_by_intensity(
                assay_1,
                relative_intensity=min_relative_fragment_intensity,
                top_n=max_fragment_number,
                return_index=True
            )

            if fragment_index is None:
                fragment_index = fragment_index_1
            else:
                fragment_index = [fragment_index[i] for i in fragment_index_1]

        if return_index:
            if fragment_index is None:
                fragment_index = list(
                    range(0, len(assay['fragments']['fragmentType']))
                )
            return fragment_index
        else:
            if fragment_index is not None:
                assay = self.filter_fragments_by_index(
                    assay, fragment_index,
                    copy=copy
                )
            return assay


    def select_quantifying_transitions(self, assay, copy=True, **kwargs):
        if copy:
            assay = deepcopy(assay)

        fragment_index = self.filter_fragments(
            assay, return_index=True,
            **kwargs
        )

        assay['fragments']['quantifyingTransition'] = [
            i in fragment_index
            for i in range(0, len(assay['fragments']['fragmentType']))
        ]
        return assay


    def filter_assay_by_peptide_info(self, assay,
                                     precursor_charge=None,
                                     min_peptide_length=None,
                                     max_peptide_length=None,
                                     modification_action=None,
                                     modification_list=None):
        if precursor_charge is not None:
            charge = assay.get('precursorCharge', None)
            if isinstance(precursor_charge, list) or \
                isinstance(precursor_charge, tuple) or \
                isinstance(precursor_charge, set):
                if charge not in precursor_charge:
                    return None
            else:
                if charge != precursor_charge:
                    return None

        sequence_length = len(assay['peptideSequence'])

        if min_peptide_length is not None:
            if sequence_length < min_peptide_length:
                return None

        if max_peptide_length is not None:
            if sequence_length > max_peptide_length:
                return None

        modifications = assay.get('modification', None)

        if modification_action is None or \
            modification_action == 'None' or \
            modification_action == 'none':
            pass

        else:
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

            if not check_modifications(
                modifications,
                action=action,
                name_list=modification_list
            ):
                return None

        return assay


    def filter_assay_by_precursor_mz(self, assay,
                                     min_precursor_mz=None,
                                     max_precursor_mz=None,
                                     swath_windows=None):
        precursor_mz = assay.get('precursorMZ', None)

        if min_precursor_mz is not None:
            if precursor_mz is None or \
                precursor_mz < min_precursor_mz:
                return None

        if max_precursor_mz is not None:
            if precursor_mz is None or \
                precursor_mz > max_precursor_mz:
                return None

        if swath_windows is not None:
            if precursor_mz is None:
                return None

            isolation_window_index = np.where(
                (swath_windows['start'] < precursor_mz) & \
                (swath_windows['end'] > precursor_mz)
            )[0]
            if len(isolation_window_index) == 0:
                return None

        return assay


    def filter_assay(self, assay,
                     min_fragment_number=6,
                     precursor_charge=None,
                     min_peptide_length=None,
                     max_peptide_length=None,
                     modification_action=None,
                     modification_list=None,
                     min_precursor_mz=None,
                     max_precursor_mz=None,
                     swath_windows=None,
                     quantifying_transition_criteria=None,
                     **kwargs):
        if min_fragment_number is not None and \
            len(assay['fragments']['fragmentType']) < min_fragment_number:
                return None

        if precursor_charge is not None or \
            min_peptide_length is not None or \
            max_peptide_length is not None or \
            modification_action is not None or \
            modification_list is not None:
            assay = self.filter_assay_by_peptide_info(
                assay,
                precursor_charge=precursor_charge,
                min_peptide_length=min_peptide_length,
                max_peptide_length=max_peptide_length,
                modification_action=modification_action,
                modification_list=modification_list
            )
            if assay is None:
                return assay

        if min_precursor_mz is not None or \
            max_precursor_mz is not None or \
            swath_windows is not None:
            assay = self.filter_assay_by_precursor_mz(
                assay,
                min_precursor_mz=min_precursor_mz,
                max_precursor_mz=max_precursor_mz,
                swath_windows=swath_windows
            )
            if assay is None:
                return assay

        assay = self.filter_fragments(
            assay, swath_windows=swath_windows, **kwargs
        )

        if min_fragment_number is not None and \
            len(assay['fragments']['fragmentType']) < min_fragment_number:
            return None

        if quantifying_transition_criteria is not None and \
            len(quantifying_transition_criteria) > 0:
            assay = self.select_quantifying_transitions(
                assay, **quantifying_transition_criteria
            )

        return assay


    def filter_assays(self, assays, return_generator=False, **kwargs):
        result = filter(
            lambda x: x is not None,
            (self.filter_assay(assay, **kwargs)
             for assay in assays)
        )

        if not return_generator:
            result = list(result)

        return result






if __name__ == '__main__':
    assay = AssayBuilder()

    print(assay.theoretical_fragments(
        sequence='ESKSSPRPTAEK',
        modification={'name': 'Phospho', 'position': 2}
    ))