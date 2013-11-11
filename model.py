'''
Created on Oct 28, 2013
'''
import numpy
from numpy.matrixlib.defmatrix import matrix
from numpy.core.fromnumeric import mean
from numpy.core.numeric import cross
from numpy.dual import norm
from transform import TransformFrame

class Model:
    def __init__(self, models=None, loop=None):
        '''
        Generates a model from list of similar loops, loops
        '''
        if models is None and loop is None:
            raise Exception("Must provide models to merge or loop")
        self.loop = loop
        self.models = models
        if models is not None:
            self.sses = None
        else:
            self.sses = (loop.l_anchor, loop.r_anchor)
        self.ssesSignature = models[0].ssesSignature if models is not None else sorted([loop.l_type, loop.r_type])
        self.phi
        self.theta
        self.positions
        self.seq = self.gen_seq() if loop is not None else None
        #todo
        
    def score(self, loop):
        """Scores how well the loop matches the Model"""
        #todo
        if(len(self.loops) == 0):
            return float("inf")
        
        total = 0
        for l in self.loops:
            total += l.closeness(loop)
        return total / len(self.loops)
    
    def compare(self, other):
        """Scores how close two Models are."""
        #todo
        if (len(other.loops) == 0):
            return float("inf")
        
        total = 0
        for loop in other.loops:
            total += self.score(loop)
        return total / len(other.loops)
    
    def compareNEW(self, other):
        sOffsetV = {c: self.sses[1][0][c] - self.sses[0][0][c] for c in 'xyz'}
        sSSE0V = {c: mean([atom[c] for atom in self.sses[0]]) for c in 'xyz'}
        
#         n = norm(sOffsetV)
#         sX = [c / n for c in sOffsetV]
#         n = norm(sSSE0V)
#         sY = [c / n for c in sSSE0V]
#         sZ = cross(sX, sY)
#         sY = cross(sX, sZ)
#         sSSE1V = {c: mean([atom[c] for atom in self.sses[1]]) for c in 'xyz'}
        
        oOffsetV = {c: other.sses[1][0][c] - other.sses[0][0][c] for c in 'xyz'}
        oSSE0V = {c: mean([atom[c] for atom in other.sses[0]]) for c in 'xyz'}
        
        
        toFrame = TransformFrame.createFromVectors(self.sses[0][0], sOffsetV, sSSE0V)
        fromFrame = TransformFrame.createFromVectors(other.sses[0][0], oOffsetV, oSSE0V)
        
        for i in range(len(self.positions)):
            point = fromFrame.transformTo(toFrame, other.positions[i])
        
#         n = norm(oOffsetV)
#         oX = [c / n for c in oOffsetV]
#         n = norm(oSSE0V)
#         oY = [c / n for c in oSSE0V]
#         oZ = cross(oX, oY)
#         oY = cross(oX, oZ)
#         oSSE1V = {c: mean([atom[c] for atom in other.sses[1]]) for c in 'xyz'}
        
    #helper method to generate a probabilistic sequence representing
    #a cluster of loops (loops must be the same length)
    def gen_seq(self):     
        num_loops = len(self.loop)
        seq = [{}]*len(self.loop[0])
        
        #iterate over each of the loop sequences, adding/updating
        #dictionary entries for each amino acid at each sequence position
        for curr_loop in self.loop:
            curr_seq = curr_loop.seq
            for i in range(len(curr_seq)):
                aa = curr_seq[i] #retrieve the current amino acid
                
                #if the amino acid's already represented at that sequence
                #location, just increment the count
                if seq[i].has_key(aa):
                    seq[i][aa] += 1.0
                
                #otherwise, create a new entry, initialized with a count
                #of 1.0
                else:
                    seq[i][aa] = 1.0
        
        #normalize the sequence (so probabilities sum to 1), then return
        return self.normalize_seq(seq, num_loops)
        
    
    #helper method to normalize a probabilistic sequence,
    #given the total number of entries per sequence location
    def normalize_seq(self, seq, n):
        for pos in seq:
            for aa in pos.keys():
                pos[aa] /= n
    
    
    #helper method to merge two probabilistic sequences,
    #given two models (each with their own representative sequence)
    #and their respective weights (i.e. number of loops represented)
    def merge_seqs(self, other, n1, n2):
        merged_seq = [{}]*len(self.seq)
        weight1 = float(n1/(n1 + n2))
        weight2 = float(n2/(n1 + n2))
        
        for i in range(len(merged_seq)):
            #add in first model's weighted sequence info
            for aa in self.seq[i].keys():
                merged_seq[i][aa] = self.seq[i].get(aa)*weight1
            
            #add in the second model's weighted sequence info
            for aa in other.seq[i].keys():
                #if there's already an entry from the first model,
                #just update the score
                if merged_seq[i].has_key(aa):
                    merged_seq[i][aa] += other.seq[i].get(aa)*weight2
                
                #if it's a new amino acid, make a new entry
                else:
                    merged_seq[i][aa] = other.seq[i].get(aa)*weight2
        
        return merged_seq
        
    def __str__(self):
        mean_displacement = numpy.mean([loop.displacement() for loop in self.loops])
        return ("Model %f%s" % (mean_displacement, "".join(["\n\t%s" % (loop) for loop in self.loops])))
