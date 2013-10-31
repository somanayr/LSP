from math import sqrt
from string import lower

AA_MAP = { 'ala' : 'A', 
           'cys' : 'C', 
           'asp' : 'D', 
           'glu' : 'E', 
           'phe' : 'F', 
           'gly' : 'G', 
           'his' : 'H', 
           'ile' : 'I', 
           'lys' : 'K', 
           'leu' : 'L',
           'met' : 'M', 
           'asn' : 'N', 
           'pro' : 'P', 
           'gln' : 'Q', 
           'arg' : 'R', 
           'ser' : 'S', 
           'thr' : 'T', 
           'val' : 'V', 
           'trp' : 'W', 
           'tyr' : 'Y'}

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
    def __init__(self, type, start, end):
        self.type=type; self.start=start; self.end=end

    def __str__(self):
        return self.type + str(self.start) + '-' + str(self.end)

class Loop:
    """A loop structure element from a PDB file: seq, atoms, start and end residue numbers."""
    def __init__(self, seq, atoms, start=0, end=None):
        self.start = start
        self.end   = end
        self.seq = seq
        self.atoms = atoms
        
    # Test...   
    def __str__(self):
        return "Loop" + str(self.start) + "-" + str(self.end)

def read_pdb(filename):
    """Get a list of SSEs and a list of atoms from a PDB file.
    Very minimal support ---
        - assumes a single chain in the PDB file
        - assumes valid format for SSEs and atoms
        - sorts SSEs by initial residue, but assumes atoms are in proper order
    @author: CBK
    """
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
    """
    Strips out the alpha helices and beta pleated sheets and returns an array of the leftover loops.
    
    @author: Travis Peters
    """
    
    # Collect SSEs & atoms from a given PDB file
    sses, atoms = read_pdb(pdb_file)

    # Group each atom with its appropriate residue
    residues = get_residues(atoms)
    
    loops     = []                # List of Loops
    total_res = atoms[-1].resnum  # Total num. of residues in PDB file    
    lstart    = sses[-1].end + 1  # Record the start residue of the next loop
    lend      = total_res         # Record the end residue of the next loop        
    
    # Insert fake SSE to extract the first loop (boundary condition)
    sses.insert(0, SSE("Term", 0, 0))

    # Iterate & extract loops from the given structure 
    prev_sse = len(sses)-2
    for sse in reversed(sses):
        if (lend - lstart) + 1 > 0:
            seq, loop_atoms = get_context(residues, lstart, lend)
            loops.append(Loop(seq, loop_atoms, lstart, lend))
        
        # Test: display SSE(s) & Loop(s) ------------ #
        #print "                Loop", lstart, lend
        #print sse
        # ------------------------------------------- #

        # Update bounds of next loop start/end for extraction
        lstart = sses[prev_sse].end + 1
        lend = sse.start - 1
        
        # Update prev SSE position 
        prev_sse -= 1
    
    # Remove the fake SSE
    del sses[0]
    
    return loops[::-1] #returns the loops in the PDB file

def get_residues(atoms):
    """
    Given a list of atoms, return a list of lists (residues list contains a residue at 
    index i & residue i contains a list of all atoms belonging to that residue.
    
    @author: Travis Peters
    """
    
    # Create an empty list where a residue will be stored at some index i.
    residues = [None] * (atoms[-1].resnum)
    
    # Temp list to store atoms for a given residue
    res_atoms = []
    
    # Begin construction: start w/ resnum of first atom
    resnum = atoms[0].resnum

    # Assign each atom to a residue
    for atom in atoms:
        # Add atoms to list for this residue
        if atom.resnum == resnum:
            res_atoms.append(atom)
        else:
            # New residue encountered - store atoms to current residue & update resnum            
            residues[(resnum-1)] = res_atoms
            resnum = atom.resnum
            res_atoms = [atom]
    
    # Make sure the final res gets stored
    residues[(resnum-1)] = res_atoms
    
    return residues

def get_context(residues, start, end):
    """
    Construct a list of atoms belonging to the given residue(s) & the corresponding sequence.
    @author: Travis Peters
    """
    # TODO: I'm not sure this is right... should there be so many of the same AAs consecutively in a sequence?
    # Given a list of residues, and a start & end residue number collect all relevant atoms
    loop_atoms = []
    for i in range(start, end+1):
        loop_atoms.extend(residues[i-1])
    
    # Construct the sequence from the atoms
    seq = ""
    for a in loop_atoms:
        seq += AA_MAP[lower(a.restype)]

    return seq, loop_atoms

if __name__ == '__main__':
    pdb_filename = "pdb_files/112L.pdb" #argv[1]
    loops = get_loops(pdb_filename)
    
    # Test: Display loops and their bounds
    def test_display_loops():
        print "\nTotal Loops: " + str(len(loops))
        for l in loops:
            print l

    # Test: Display each loop's sequence
    def test_display_loop_seqs():
        print "\nLoop Structure Sequences:"
        for l in loops:
            print l.get_seq()

    # Test: Display each loop's residues/atoms
    def test_display_loop_residues():
        print "\nLoop Residues:"
        for l in loops:
            resnum = l.atoms[0].resnum-1
            print l

            for a in l.atoms:
                if a.resnum == resnum:
                    print a,
                else:
                    print "\nRes. " + str(resnum+1) + ": " + str(a),
                    resnum += 1
            print "\n"

    #################################
    #         Test Methods          #
    #################################
    test_display_loops()
    #test_display_loop_seqs()
    #test_display_loop_residues()
    

# def get_loops(pdb_file):
#     """Strips out the alpha helices and beta pleated sheets and returns an array
#     of the leftover loops
#     
#     pdb_file -- the name of the PDB file containing the protein
#     @author: Travis Peters
#     """
#     
#     """Placeholder code for testing by Ryan Amos"""
#     sses, atoms = read_pdb(pdb_file)
# #     resoff = atoms[0]
#     loop_start = 1
#     loops = []
#     for sse in sses:
#         if(loop_start < sse.start):
#             subset = atoms[loop_start - 1:sse.start - 1]
# #             loops.append(Loop([atom.restype for atom in subset], [[atom.x, atom.y, atom.z] for atom in subset]))
#             loop = Loop(subset)
#             if(len(loop.seq) >= 2):
#                 loops.append(loop)
#         loop_start = sse.end + 1
#     
# #     subset = atoms[loop_start - resoff:]
# #     if(len(subset) != 0):
# #         loops.append(Loop([atom.restype for atom in subset], [[atom.x, atom.y, atom.z] for atom in subset]))
#     #todo
#     return loops #returns the loops in the pdb file
