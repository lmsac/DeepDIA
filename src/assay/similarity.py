import itertools
import math


def dot_product(x, y):
    prod1 = sum(a * a for a in x)
    prod2 = sum(b * b for b in y)
    if prod1 == 0 or prod2 == 0:
        return 0
    
    return sum(a * b for a, b in zip(x, y)) / math.sqrt(prod1 * prod2)


def align_fragments_by_annotation(spectra, ignore_none_annotation=True):
    fragment_annotation = set(itertools.chain.from_iterable(
        spec['fragments']['fragmentAnnotation']
        for spec in spectra
    ))
    if ignore_none_annotation:
        fragment_annotation = filter(
            lambda x: x is not None, 
            fragment_annotation
        )
    
    fragment_annotation = list(fragment_annotation)
    
    return [
        [
            next(
                (i for i, x in \
                 enumerate(spec['fragments']['fragmentAnnotation'])
                 if x == annot), 
                None
            )
            for spec in spectra
        ]
        for annot in fragment_annotation
    ]
    
    
class SimilarityScorer:
    def __init__(self, 
                 alignment_func=align_fragments_by_annotation,
                 similarity_func=dot_product):
        if alignment_func is None:
            alignment_func=align_fragments_by_annotation
        self.alignment_func = alignment_func
        
        if similarity_func is None:
            similarity_func = dot_product
        self.similarity_func = similarity_func
    
    
    def align_fragments(self, spectra):
        return self.alignment_func(spectra)
    
    
    def pairwise_similarity(self, spectra):
        index = self.align_fragments(spectra)
        
        intensity = [
            [
                (spectra[i]['fragments']['fragmentIntensity'][t[i]] \
                 if t[i] is not None \
                 else 0)
                for t in index
            ]
            for i in range(len(spectra))
        ]
               
        return [
            [
                self.similarity_func(intensity[i], intensity[j])
                for j in range(i + 1, len(spectra))
            ]
            for i in range(len(spectra) - 1)
        ]
            
    
    def similarity(self, spectrum1, spectrum2):
        return self.pairwise_similarity([spectrum1, spectrum2])[0][0]
    
    