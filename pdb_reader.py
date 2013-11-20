from math import sqrt
from string import lower
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import sys
import os
from protein_structure import SSE, Atom, Loop

# A minimum number of consecutive residues must have the same sse_code
# in order for a group of residues to be declared a HELIX or SHEET.
# NOTE: if some contiguous sequence of residues is not declared to be a
#  HELIX or a SHEET it will later be considered a Loop.
REQUIRED_SSE_RES_NUM = 3

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

def read_pdb(filename):
    """
    Get a list of SSEs and a list of atoms from a PDB file.
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

def get_ca_chain(atoms):
    """
    Get just the CA atoms from the list.
    Try to make them a legal chain by filling in placeholder CAs for missing
    residues, and choosing just one alternative when there are multiple possible
    CAs for a residue.
    @author: CBK
    """
    cas = [a for a in atoms if a.atomtype=='CA']
    r1 = cas[0].resnum
    i = 0
    while i<len(cas):
        if cas[i].resnum > r1+i:
            #print('missing ', r1+i)
            cas.insert(i,Atom(i,'---','CA',None,None,None))
            i += 1
        elif cas[i].resnum < r1+i:
            #print('duplicate at ', r1+i-1)
            del cas[i]
        else:
            i += 1
    return cas

def dssp_sse_extract_from_pdb(filename):
    """
    Construct a list of SSEs from an un-annotated PDB file.
    @author: Travis Peters
    
    # Test... ##################################
    #aa_type = dssp[0][key][0]
    #sse_code = dssp[0][key][1]        
    #x = float(dssp[0][key][2])
    #y = float(dssp[0][key][3])
    #z = float(dssp[0][key][4])

    #print(str(key) + ":" + str(dssp[0][key]))
    #print "  AA Type  = " + str(aa_type),
    #print "\n  SSE Code = " + str(sse_code),
    #print "\n  Location = " + str((x,y,z))
    ############################################
    """

    # Assumes dssp executable is located in root of project directory
    if sys.platform[0:2] == "win":
        DSSP_EXEC = "dssp"
    else:
        DSSP_EXEC = "./dssp"

    # NOTE: The DSSP codes for secondary structure used here are: 
    # - H        Alpha helix (4-12) 
    # - G        3-10 helix 
    # - I        pi helix 
    # - B        Isolated beta-bridge residue 
    # - E        Strand 
    # - T        Turn 
    # - S        Bend 
    # - -        None 
    HELIX = ['H', 'G', 'I']
    SHEET = ['B', 'E']

    # DSSP call returns a dictionary that maps (chainid, resid) to 
    # (amino acid type, secondary structure code, and accessibility).
    dssp = dssp_dict_from_pdb_file(filename, DSSP_EXEC)
        
    sses = []
    sse_start = None
    sse_type = None
    res_count = 0
    for key in dssp[1]:

        # Extract the residue number
        resnum = key[1][1]

        # Extract SSE code for a residue
        sse_code = dssp[0][key][1]
        
        # Record SSEs by examining sse_codes of consecutive residues
        if sse_code in HELIX:
            if sse_type == 'HELIX':
                res_count += 1
            else:
                # Did we just detect an SSE
                if res_count >= REQUIRED_SSE_RES_NUM:
                    sses.append( SSE(sse_type, sse_start, resnum-1) )

                # Start recording a new SSE
                res_count = 0
                sse_start = resnum
                sse_type = 'HELIX'
            
        elif sse_code in SHEET:
            if sse_type == 'SHEET':
                res_count += 1
            else:
                # Did we just detect an SSE
                if res_count >= REQUIRED_SSE_RES_NUM:
                    sses.append( SSE(sse_type, sse_start, resnum-1) )

                # Start recording a new SSE
                res_count = 0
                sse_start = resnum
                sse_type = 'SHEET'
            
        else:
            if not (sse_type == None):
                # Did we just detect an SSE
                if res_count >= REQUIRED_SSE_RES_NUM:
                    sses.append( SSE(sse_type, sse_start, resnum-1) )

            # sse_code suggests we are not detecting an SSE
            res_count = 0
            sse_start = resnum
            sse_type = None
            
    return sses

def get_loops(pdb_file, numAnchors=1):
    """
    Strips out alpha helices and beta pleated sheets and returns an array of the leftover loops.
    @author: Travis Peters
    """
    
    # Collect SSEs & atoms from a given PDB file
    sses, atoms = read_pdb(pdb_file)

    # Check: If no SSEs were found, file must not be annotated - use BioPython/DSSP to extract SSEs
    if len(sses) == 0:
        sses = dssp_sse_extract_from_pdb(pdb_file)
    
    # Extract the c-alphas 
    cas = get_ca_chain(atoms)

    # Group each atom with its appropriate residue
    residues = get_residues(cas)
    
    loops     = []                # List of Loops
    total_res = atoms[-1].resnum  # Total num. of residues in PDB file    
    lstart    = sses[-1].end + 1  # Record the start residue of the next loop
    lend      = total_res         # Record the end residue of the next loop        
    
    # Insert fake SSE to extract the first loop (boundary condition)
    sses.insert(0, SSE("Term", 0, 0))

    # Iterate & extract loops from the given structure 
    current_sse = len(sses)-1
    for sse in reversed(sses):
        # Determine the type of SSEs to the left/right of loop
        rtype = ltype = None
        # Set right anchor/anchor type
        if not (current_sse >= len(sses)-1):
            rtype = sses[current_sse+1].type
            
        # Set left anchor/anchor type
        if not (current_sse <= 0):
            ltype = sses[current_sse].type
        
        # Construct loop structure        
        if (lend - lstart) + 1 > 0:
            seq, loop_res, la, ra = get_context(residues, lstart, lend, numAnchors)
            loops.append(Loop(seq, loop_res, lstart, lend, (la,ltype), (ra,rtype)))
        
        # Test: display SSE(s) & Loop(s) ------------ #            
        #print "                Loop", lstart, lend, (ltype,rtype)
        #print sse
        # ------------------------------------------- #

        # Update bounds of next loop start/end for extraction
        lstart = sses[current_sse-1].end + 1
        lend = sse.start - 1
        
        # Update prev SSE position 
        current_sse -= 1
    
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

def get_context(residues, start, end, anchors=1):
    """
    Construct a list of atoms belonging to the given residue(s) & the corresponding sequence.
    @author: Travis Peters
    """

    # Given a list of residues, and a start & end residue number, define the loop anchors    
    l_anchor = []
    a = 1
    while residues[start-a][0].resnum > 1 and a <= anchors:
        l_anchor.append(residues[start-(a+1)][0])
        a += 1
        
    r_anchor = []
    a = 0
    while residues[end-1][0].resnum < len(residues) and a < anchors:
        r_anchor.append(residues[end+a][0])
        a += 1
    
    # Given a list of residues, and a start & end residue number collect all relevant atoms
    loop_res = []
    for i in range(start, end+1):
        loop_res.extend(residues[i-1])
    
    # Construct the sequence from the atoms
    seq = ""
    for a in loop_res:
        seq += AA_MAP[lower(a.restype)]

    return seq, loop_res, l_anchor[::-1], r_anchor

###############################################################################
#                                   Tests                                     #
###############################################################################

# Test: Display each loop's sequence
def test_display_all(loops):
    print "\nLoop Structures:"
    for l in loops:
        print str(l) + "; Sequence: " + l.seq

# Test: Display loops and their bounds
def test_display_loops(loops):
    print "\nTotal Loops: " + str(len(loops))
    for l in loops:
        print l

# Test: Display each loop's sequence
def test_display_loop_seqs(loops):
    print "\nLoop Structure Sequences:"
    for l in loops:
        print l.seq

# Test: Display each loop's residues/atoms
def test_display_loop_residues(loops):
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

# Test: Display the SSEs that were extracted from an un-annotated PDB file
def test_dssp_pdb_extract(numFilesToRead=5):
    pdb_dir = "pdb"
    files = os.listdir(pdb_dir)
    processed = 1

    # For each file in the pdb directory - extract SSEs & report
    for f in files:
        if ".pdb" in f:            
            sses = dssp_sse_extract_from_pdb(pdb_dir+"/"+f)
                
            print "There were " + str(len(sses)) + " SSEs found in '" + f + "':"
            for s in sses:
                print s
            print
            
            if processed >= numFilesToRead:
                break
        
            processed += 1

# Test: Display the anchor atoms for each loop in the list of loops.
def test_get_loop_anchors(loops):
    for loop in loops:
        print loop
        
        print "  Left:",
        if len(loop.l_anchor) > 0:
            for a in loop.l_anchor:
                print str(a),
            print
        else:
            print " <None> "
            
        print "  Right:",
        if len(loop.r_anchor) > 0:
            for a in loop.r_anchor:
                print str(a),
            print
        else:
            print " <None> "
            
if __name__ == '__main__':

#    pdb_filename = argv[1]
    pdb_dir = "pdb/"
    files = os.listdir(pdb_dir)
    pdb_filename = pdb_dir + files[0]
    
    loops = get_loops(pdb_filename, numAnchors=3)

    #################################
    #         Test Methods          #
    #################################
    #test_display_all(loops)
    #test_display_loops(loops)
    #test_display_loop_seqs(loops)
    #test_display_loop_residues(loops)
    #test_dssp_pdb_extract(1)
    test_get_loop_anchors(loops)
     
