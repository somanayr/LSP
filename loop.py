'''
Created on Oct 28, 2013

'''
from pdb_reader import read_pdb

class Loop:
    def __init__(self, seq, atoms, start=0, end=None):
        if end is None:
            self.seq = seq[start:]
            self.atoms = atoms[start:]
        else:
            self.seq = seq[start:end]
            self.atoms = atoms[start:end]
            
    def closeness(self, other):
        """Determines a score for how closely related two loops are. Returns a decimal
        number [0,inf], 0 is equal, inf is completely different
        """
        #todo
        return 0.0
            
def get_loops(pdb_file):
    """Strips out the alpha helices and beta pleated sheets and returns an array
    of the leftover loops
    
    pdb_file -- the name of the PDB file containing the protein
    """
    pdb = read_pdb(pdb_file)
    #todo
    return [] #returns the loops in the pdb file
