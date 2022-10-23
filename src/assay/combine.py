import itertools
import numpy as np

from .modseq import stringify_modification
from .values import get_assay_values


class AssayCombiner():
    def __init__(self, group_key=None):
        if group_key is None:
            group_key = peptide_group_key()
        self.group_key = group_key


    def combine(self, *assays, return_generator=False):
        return self.remove_redundant(
            itertools.chain.from_iterable(assays),
            return_generator=return_generator
        )


    def remove_redundant(self, assays, return_generator=False):
        result = (
            self.combine_replicates(spectra)
            for spectra in self.group_replicates(assays)
        )
        result = (x for x in result if x is not None)

        if not return_generator:
            result = list(result)

        return result


    def group_replicates(self, assays):
        def get_key(assay):
            d = get_assay_values(assay, self.group_key)
            return tuple(map(str, d.values()))

        return (
            list(v)
            for k, v in itertools.groupby(
                sorted(assays, key=get_key),
                key=get_key
            )
        )


    def combine_replicates(self, spectra):
        if len(spectra) == 0:
            return None

        return spectra[0]


class BestReplicateAssayCombiner(AssayCombiner):
    def __init__(self,
                 group_key=None,
                 score='score', higher_score_better=True):
        super(BestReplicateAssayCombiner, self) \
            .__init__(group_key=group_key)
        self.score = score
        self.higher_score_better = higher_score_better


    def combine_replicates(self, spectra):
        if len(spectra) == 0:
            return None

        score = [
            spec['metadata'][self.score]
            for spec in spectra
        ]
        if self.higher_score_better:
            index = np.argmax(score)
        else:
            index = np.argmin(score)

        return spectra[index]


def peptide_group_key(within_run=False):
    group_key = [
        'peptideSequence',
        {
            'name': 'modification',
            'path': 'modification',
            'convert': stringify_modification
        },
        'precursorCharge'
    ]

    if within_run:
        group_key.insert(0, {
            'name': 'filename',
            'path': ['metadata', 'file']
        })

    return group_key

