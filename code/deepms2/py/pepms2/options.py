class PeptideMS2Options:
    
    def __init__(self, max_sequence_length=50):
        self.max_sequence_length = max_sequence_length

        self.amino_acids = 'ARNDCEQGHILKMFPSTWYV'
        self.modifications = {}
        
        self.fragments = [
            ('b1', False), ('b2', False), 
            ('bn1', False), ('bn2', False), 
            ('bo1', False), ('bo2', False), 
            ('y1', True), ('y2', True), 
            ('yn1', True), ('yn2', True), 
            ('yo1', True), ('yo2', True)
        ]


    def amino_acid_size(self):
        return 20 + sum([
            len(mods) for mods in self.modifications.values()
        ])        
        

    def intensity_size(self):
        return len(self.fragments)
            

    @staticmethod
    def default():
        return PeptideMS2Options()

