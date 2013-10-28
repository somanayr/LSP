'''
Created on Oct 28, 2013

@author: ramos
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
            
def get_loops(pdb_file):
    """Strips out the alpha helices and beta pleated sheets and returns an array
    of the leftover loops
    
    pdb_file -- the name of the PDB file containing the protein
    """
    pdb = read_pdb(pdb_file)
    return [] #returns the loops in the pdb file
