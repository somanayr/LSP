'''
Created on Oct 28, 2013
'''
from pdb_reader import read_pdb
import sys
from SubstitutionMatrix import blosum62
from numpy.lib.scimath import sqrt
import numpy

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
        number [0,inf), 0 is equal, the higher the number, the less similar
        @author Ryan Amos
        """
        
#         v1 = [self.atoms[0], [self.atoms[-1][i] - self.atoms[0][i] for i in range(len(self.atoms[0]))]]
#         v2 = [other.atoms[0], [other.atoms[-1][i] - other.atoms[0][i] for i in range(len(other.atoms[0]))]]
#         v2_to_origin = numpy.matrix([1, 0, 0, -v2[0][0]], [0, 1, 0, -v2[0][1]], [0, 0, 1, -v2[0][2]], [0, 0, 0, 1])
#         
#         a,b,nw_score = nw(other.seq, other.seq, blosum62, -4)
#         selfIdx = 0
#         otherIdx = 0
#         score = 0.0
#         for idx in range(len(a)):
#             self_point = other_point = None
#             if(idx != '-'):
#                 self_point = self.atoms[selfIdx]
#                 selfIdx += 1
#             if idx != '-':
#                 other_point = other.atoms[otherIdx]
#                 otherIdx += 1
#             if other_point != None and self_point != None:
#                 score += rmsd(self_point, other_point)
         
        return rmsd([rmsd(self.atoms[0],self.atoms[-1])],[rmsd(other.atoms[0], other.atoms[-1])])
    
def rmsd(p1, p2):
    if len(p1) != len(p2):
        raise Exception("Degree of points must be the same")
    total = 0
    for t in len(p1):
        d = (p1[t] - p2[t])
        total += d * d
    return sqrt(total / len(p1)) #todo remove sqrt?
    
    
    

def make_2d_list(nrows, ncols):
    """Create a list of lists, of the given dimensions, filled with None."""
    return [[None]*ncols for i in range(nrows)]
    
def score(x, y, substScore, matrix, gap, sw):
    """Determines the score at x,y in matrix given substitution score substScore and
    sw determining if we should use Smith-Waterman scoring or Needleman-Wunsch scoring
    (sw is true for Smith-Waterman, false for Needleman-Wunsch)"""
    if x is 0 and y is 0: #no score can be made!
        return
    
    #scoreMatch - score for diagonal match/substitution
    #scoreXSkip - score for a horizontal gap
    #scoreYSkip - score for vertical gap
    #use -sys.maxint - 1 is the lowest possible integer value in python -(2^31) as an exaggerated low score
    scoreMatch = scoreXSkip = scoreYSkip = -sys.maxint - 1
    if x > 0 : #if horizontal gap possible
        scoreXSkip = matrix[x - 1][y][0] + gap
        if y > 0 : #if horizontal & vertical gap possible implies diagonal match/sub possible
            scoreMatch = matrix[x - 1][y - 1][0] + substScore
    
    if y > 0: #if vertical ap possible
        scoreYSkip = matrix[x][y-1][0] + gap
        
    if scoreMatch >= scoreXSkip and scoreMatch >= scoreYSkip: #if diagonal has highest score
        matrix[x][y] = (scoreMatch, 'm')
    elif scoreXSkip >= scoreYSkip and scoreXSkip >= scoreMatch: #if horizontal has highest score
        matrix[x][y] = (scoreXSkip, 'x')
    elif scoreYSkip >= scoreXSkip and scoreYSkip >= scoreMatch: #if vertical has highest score
        matrix[x][y] = (scoreYSkip, 'y')
    else: #this state should not be possible!
        raise Exception('Invalid state', scoreMatch, scoreXSkip, scoreYSkip)
        
    #if doing smith-waterman, need to prevent from going negative
    if(sw and matrix[x][y][0] < 0): 
        matrix[x][y] = (0, matrix[x][y][1])
        
def makeScoreMatrix(a, b, subst, gap, sw=False):
    """Makes the score matrix for sequences a,b with math/substitution value/penalty matrix subst,
    gap penalty gap, and Smith-Waterman/Needleman-Wunsch state sw (sw is true for Smith-Waterman and
    false for Needleman-Wunsch)
    
    Returns the matrix (2d list)"""
    
    # 1 vertical slot for each element in b, 1 horizontal slot for each element in a, +1 slot for a starting score
    scores = make_2d_list(len(a)+1,len(b)+1) 
    a_shorter = len(a) < len(b)
    
    v = 0
    #seed score of 0
    scores[0][0] = (v, 'e')
    
    #fill in sides
    #horizontal side
    for x in range(1, len(a) + 1):
        if(not sw): 
            v += gap
        scores[x][0] = (v, 'x')
    
    #vertical side
    v = 0
    for y in range(1, len(b) + 1):
        if(not sw):
            v += gap
        scores[0][y] = (v, 'y')
    
    #fill in matrix by moving along the diagonal with slope 1 and filling it everything below (y >= d) and to the right (x >= d) of that diagonal
    #d is the diagonal position, which starts at 0, but we already have the 0 layer filled out, so +1 to all positions in the matrix
    for d in range(len(a if a_shorter else b)):
        v2 = b[d] #vertical protein/DNA/etc value
        for x in range(d, len(a)) :
            v1 = a[x] #horizontal value
            score(x+1, d+1, subst[v1,v2], scores, gap, sw)
        v1 = a[d] #horizontal value
        for y in range(d + 1, len(b)):
            v2 = b[y] #vertical value
            score(d+1, y+1, subst[v1,v2], scores, gap, sw)
    
    return scores;

def nw(a, b, subst, gap):
    """Global alignment of strings a and b, with scoring under the given substitution matrix and gap penalty (a negative number)."""
    # TODO: your code here
    
    scores = makeScoreMatrix(a, b, subst, gap)
    
    x = len(a)
    y = len(b)
    
    aSub, bSub = findBestPath(x, y, scores, a, b)
    
    return aSub, bSub, scores[len(a)][len(b)][0]  # aligned versions of a and b

def sw(a, b, subst, gap):
    """Local alignment of strings a and b, with scoring under the given substitution matrix and gap penalty (a negative number)."""
    # TODO: your code here
    scores = makeScoreMatrix(a, b, subst, gap, sw=True)
    
    #find max score and position of max score
    x = 0
    y = 0
    maxS = 0
    for i in range(len(a) + 1):
        for j in range(len(b) + 1):
            score = scores[i][j][0]
            if maxS < score : 
                x = i
                y = j
                maxS = score
                
    aSub, bSub = findBestPath(x, y, scores, a, b, sw=True)
    
    return aSub, bSub, maxS  # aligned versions of a and b

def findBestPath(x, y, scores, a, b, sw=False):
    """"Traces backwards from x,y through scores to generate equal length sequences from a,b to describe the collision. If sw is true, stops when the score is 0. If nw is true, stops when out of instructions"""
    #A,B list/subset/where the letters for a go, where the letters for b go
    aSub = []
    bSub = []
    
    #Stop when score is 0 (sw) or when instructions end (i.e. second argument of tuple is 'e')
    while ((scores[x][y][0] is not 0) if sw else (scores[x][y][1] is not 'e')):
        scoreType = scores[x][y][1] #where the score came from
        if scoreType is 'm': #diagonal match/substitution
            aSub.insert(0, a[x - 1])
            bSub.insert(0, b[y - 1])
            x -= 1
            y -= 1
        elif scoreType is 'y': #vertical gap
            aSub.insert(0, '-')
            bSub.insert(0, b[y - 1])
            y -= 1
        elif scoreType is 'x': #horizontal gap
            aSub.insert(0, a[x - 1])
            bSub.insert(0, '-')
            x -= 1
        pass
    
    return aSub, bSub
