class SequenceDatabase:

    def __init__(self, sequences = None):
        
        if sequences == None:
            self.sequences = []
        else:
            self.sequences = sequences

    def add_sequence(self, sequence):
        self.sequences.append(sequence)
    
    def get_sequences(self):
        return self.sequences