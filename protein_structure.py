'''
Created on Oct 28, 2013
'''
from numpy.lib.scimath import sqrt

class Loop:
    """A loop structure element from a PDB file: seq, atoms, start and end residue numbers."""
    def __init__(self, seq, atoms, start, end, l_anchor=None, r_anchor=None):
        self.start = start
        self.end = end
        self.seq = seq
        self.atoms = atoms
        self.l_anchor = l_anchor[0] # list of CAs
        self.l_type = l_anchor[1]   # type of l_anchor SSE
        self.r_anchor = r_anchor[0] # list of CAs
        self.r_type = r_anchor[1]   # type of r_anchor SSE
        
    def __str__(self):
        # Left anchor(s)
        la_str = ""
        if len(self.l_anchor) > 0:
            for a in self.l_anchor:
                la_str += str(a) + ","
            la_str = la_str[0:-1]
        else:
            la_str = "-"

        # Right anchor(s)
        ra_str = ""
        if len(self.r_anchor) > 0:
            for a in self.r_anchor:
                ra_str += str(a) + ","
            ra_str = ra_str[0:-1]
        else:
            ra_str = "-"

        return "Loop" + str(self.start) + "-" + str(self.end) + ", Anchors[" + la_str + " # " + ra_str + "], Types" + str((self.l_type,self.r_type))

class Atom:
    """
    An atom from a pdb file: residue number, residue type, atom type, and x,y,z coordinates.
    @author: CBK
    """
    def __init__(self, resnum, restype, atomtype, x, y, z):
        self.resnum=resnum; 
        self.restype=restype; 
        self.atomtype=atomtype
        self.x=x; self.y=y; self.z=z

    def __str__(self): 
        return self.restype+str(self.resnum)+'_'+self.atomtype

    def dist(self, a2):
        """Euclidean distance between self and atom a2."""
        if self.x==None or a2.x==None: return None
        return sqrt((self.x-a2.x)*(self.x-a2.x)+(self.y-a2.y)*(self.y-a2.y)+(self.z-a2.z)*(self.z-a2.z))

class SSE:
    """
    A secondary structure element from a PDB file: HELIX or SHEET type, and start and end residue numbers.
    @author: CBK
    """
    def __init__(self, ssetype, start, end):
        self.type=ssetype; self.start=start; self.end=end

    def __str__(self):
        return self.type + str(self.start) + '-' + str(self.end)

def rmsd(p1, p2):
    if len(p1) != len(p2):
        raise Exception("Degree of points must be the same")
    total = 0
    for t in range(len(p1)):
        d = (p1[t] - p2[t])
        total += d * d
    return sqrt(total / len(p1)) #todo remove sqrt?