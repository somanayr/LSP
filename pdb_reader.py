from math import sqrt

'''
@author: CBK
'''
from loop import Loop

class Atom:
    """An atom from a pdb file: residue number, residue type, atom type,
and x,y,z coordinates."""
    def __init__(self, resnum, restype, atomtype, x, y, z):
        self.resnum=resnum; self.restype=restype; self.atomtype=atomtype
        self.x=x; self.y=y; self.z=z

    def __str__(self): 
        return self.restype+str(self.resnum)+'_'+self.atomtype

    def dist(self, a2):
        """Euclidean distance between self and atom a2."""
        if self.x==None or a2.x==None: return None
        return sqrt((self.x-a2.x)*(self.x-a2.x)+(self.y-a2.y)*(self.y-a2.y)+(self.z-a2.z)*(self.z-a2.z))

class SSE:
    """A secondary structure element from a PDB file: HELIX or SHEET type,
and start and end residue numbers."""
    def __init__(self, type, start, end):
        self.type=type; self.start=start; self.end=end

    def __str__(self):
        return self.type + str(self.start) + '-' + str(self.end)

def read_pdb(filename):
    """Get a list of SSEs and a list of atoms from a PDB file.
Very minimal support ---
  assumes a single chain in the PDB file
  assumes valid format for SSEs and atoms
  sorts SSEs by initial residue, but assumes atoms are in proper order"""
    sses = []
    atoms = []
    for line in open(filename,'r').readlines():
        if line[0:5]=='HELIX':
            sses.append(SSE('HELIX', int(line[21:25]), int(line[34:37])))
        elif line[0:5]=='SHEET':
            sses.append(SSE('SHEET', int(line[22:26]), int(line[34:38])))
        elif line[0:4]=='ATOM' or line[0:6]=='HETATM':
            atoms.append(Atom(int(line[23:26].strip()),line[17:20].strip(),line[13:16].strip(),float(line[31:39]),float(line[39:47]),float(line[47:55])))

    sses.sort(key=lambda sse: sse.start)
    return sses, atoms


def get_loops(pdb_file):
    """Strips out the alpha helices and beta pleated sheets and returns an array
    of the leftover loops
    
    pdb_file -- the name of the PDB file containing the protein
    @author: Travis Peters
    """
    
    """Placeholder code for testing by Ryan Amos"""
    sses, atoms = read_pdb(pdb_file)
    resoff = atoms[0].resnum - 0
    loop_start = atoms[0].resnum
    loops = []
    for sse in sses:
        if(loop_start < sse.start):
            subset = atoms[loop_start-resoff:sse.start-resoff]
            loops.append(Loop([atom.restype for atom in subset], [[atom.x, atom.y, atom.z] for atom in subset]))
        loop_start = sse.end + 1
    
#     subset = atoms[loop_start - resoff:]
#     if(len(subset) != 0):
#         loops.append(Loop([atom.restype for atom in subset], [[atom.x, atom.y, atom.z] for atom in subset]))
    #todo
    return loops #returns the loops in the pdb file
