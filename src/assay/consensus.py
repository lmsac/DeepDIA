import numpy as np
import copy

from .combine import AssayCombiner
from .similarity import SimilarityScorer

class ConsensusAssayCombiner(AssayCombiner):
    def __init__(self, 
                 similarity_scorer=None,
                 group_key=None,
                 replicate_similarity_threshold=0.5,                 
                 replicate_weight=None,
                 maximum_replicates_number=100,
                 peak_quorum=0.6):
        super(ConsensusAssayCombiner, self) \
            .__init__(group_key=group_key)
            
        if similarity_scorer is None:
            similarity_scorer = SimilarityScorer()
        self.similarity_scorer = similarity_scorer
            
        self.replicate_similarity_threshold = replicate_similarity_threshold
        self.replicate_weight = replicate_weight
        self.maximum_replicates_number = maximum_replicates_number
        self.peak_quorum = peak_quorum
        
    
    def combine_replicates(self, spectra):              
        if len(spectra) == 0:
            return None
        
        if len(spectra) == 1:
            return spectra[0]
        
        if self.replicate_similarity_threshold is not None:        
            filtered_index = self.remove_dissimilar_replicates(spectra)
            spectra = [spectra[i] for i in filtered_index]        
            if len(spectra) == 0:
                return None
        
        if self.replicate_weight is not None:
            weight = [
                spec['metadata'][self.replicate_weight]
                for spec in spectra
            ]
            sorted_index = np.argsort(weight)[::-1]
        else:
            weight = [1 for x in range(len(spectra))]
            sorted_index = np.array(list(range(len(spectra))))
        
        if self.maximum_replicates_number is not None:
            sorted_index = sorted_index[:self.maximum_replicates_number]
        
        spectra = [spectra[i] for i in sorted_index]     
        weight = [weight[i] for i in sorted_index]     
        if len(spectra) == 0:
            return None        
        
        fragment_index = self.similarity_scorer.align_fragments(spectra)
        if self.peak_quorum is not None:
            fragment_index = [
                t
                for t in fragment_index
                if sum(map(lambda x: x is not None, t)) / len(spectra) \
                    > self.peak_quorum
            ]
            
        if len(fragment_index) == 0:
            return None
        
        fragments = {
            k: []
            for k in spectra[0]['fragments'].keys()
        }
        fragment_intensity = fragments['fragmentIntensity']
        fragment_mz = fragments.get('fragmentMZ', None)
        
        for t in fragment_index:
            intensity = sum((
                spectra[i]['fragments']['fragmentIntensity'][x] * weight[i]
                for i, x in enumerate(t)
                if x is not None
            )) / sum((
                weight[i]
                for i, x in enumerate(t)
                if x is not None
            ))
            fragment_intensity.append(intensity)
            
            if fragment_mz is not None:
                mz = sum((
                    spectra[i]['fragments']['fragmentMZ'][x] * weight[i]
                    for i, x in enumerate(t)
                    if x is not None
                )) / sum((
                    weight[i]
                    for i, x in enumerate(t)
                    if x is not None
                ))
                fragment_mz.append(mz)
    
            for k in fragments.keys():
                if k == 'fragmentIntensity' or k == 'fragmentMZ':
                    continue
                
                for i, x in enumerate(t):
                    if x is not None:
                        array = spectra[i]['fragments'].get(k, None)
                        if array is not None:
                            fragments[k].append(array[x])
                        else:
                            fragments[k].append(None)
                        break
                
                if len(fragments[k]) < len(fragment_intensity):
                    fragments[k].append(None)
                  
        result = copy.deepcopy(spectra[0])
        result.update({
            'fragments': fragments
        })
        return result
    
    
    def remove_dissimilar_replicates(self, spectra):
        similarity = self.similarity_scorer.pairwise_similarity(spectra)
        
        index = []
        for i in range(len(spectra)):
            score = np.median([
                (similarity[i][j - i - 1] \
                 if j > i \
                 else similarity[j][i - j - 1])
                for j in range(len(spectra))
                if j != i
            ])
            if score > self.replicate_similarity_threshold:
                index.append(i)
        
        return index
    
    