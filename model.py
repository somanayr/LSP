'''
Created on Oct 28, 2013
'''
import numpy
from numpy.matrixlib.defmatrix import matrix
from numpy.core.fromnumeric import mean
from numpy.core.numeric import dot
from numpy.dual import norm
from transform import TransformFrame, Vec
from loop import rmsd
from numpy.lib.scimath import arccos, sqrt
import transform
import sys
from numpy.ma.core import negative

class Model:
    def __init__(self, parents, positions, sses, ssesSignature, seq, size):
        '''
        Generates a model from list of similar loops, loops
        '''
        self.parents = parents
        self.positions = positions
        self.sses = sses
        self.ssesSignature = ssesSignature
        self.seq = seq
        self.size = size
        self.__loops = None #lazy loaded, NOT SAFE
        
    @classmethod
    def fromLoop(cls, loop):
        """Returns a Model representing the loop"""
        sses = sorted([(loop.l_type, loop.l_anchor), (loop.r_type, loop.r_anchor)])
        return Model([loop], loop.atoms, [sses[0][1], sses[1][1]], [sses[0][0], sses[1][0]], Model.gen_seq([loop.seq]), 1);
        
    @classmethod
    def fromModels(cls, m1, m2):
        """Returns a Model representing a merge of m1 and m2"""
        #find necessary vectors
        
        
        if m1.ssesSignature != m2.ssesSignature:
            return Exception("Loop SSE signature mismatch!")
        if len(m1.positions) != len(m2.positions):
            raise Exception("Loop length mismatch!")

        sOffsetV = [m1.sses[1][0].__dict__[c] - m1.sses[0][0].__dict__[c] for c in 'xyz']
#         sOffsetV = [m1.positions[-1].__dict__[c] - m1.positions[0].__dict__[c] for c in 'xyz']
#         sSSEV = [m1.positions[0].__dict__[c] - mean([atom.__dict__[c] for atom in m1.sses[0]]) for c in 'xyz']
        sSSEV = Model.get_sse_vector(m1.sses[0], m1.positions[0])
        
        oOffsetV = [m2.sses[1][0].__dict__[c] - m2.sses[0][0].__dict__[c] for c in 'xyz']
#         oOffsetV = [m2.positions[-1].__dict__[c] - m2.positions[0].__dict__[c] for c in 'xyz']
#         oSSE0V = [m2.positions[0].__dict__[c] - mean([atom.__dict__[c] for atom in m2.sses[0]]) for c in 'xyz']
#         oSSE1V = [m2.positions[-1].__dict__[c] - mean([atom.__dict__[c] for atom in m2.sses[1]]) for c in 'xyz']
        oSSE0V = Model.get_sse_vector(m1.sses[0], m1.positions[0])
        oSSE1V = Model.get_sse_vector(m1.sses[1], m1.positions[-1])
        
        oSSEV = oSSE0V
        
        sOrigin = m1.sses[0][0]
        oOrigin = m2.sses[0][0]
        
        otherPositions = m2.positions
        
        if m1.ssesSignature[0] == m1.ssesSignature[1]: #if we have helix/helix or sheet/sheet, we don't know which sheet is "first", so decide that based on which has a better score
            scoreNormal = m1.__compute_scores(m2.sses[0], m2.sses[1], m2.positions)
            revPos = [k for k in reversed(m2.positions)]
            scoreReversed = m1.__compute_scores(m2.sses[1], m2.sses[0], revPos)
            if(scoreReversed < scoreNormal):
                otherPositions = revPos
                oSSEV = oSSE1V
                oOffsetV = negative(oOffsetV)
                oOrigin = m2.sses[1][0]
            revPos = None
        
        sFrame = TransformFrame.createFromVectors(sOrigin, transform.Vec.from_array(sOffsetV), transform.Vec.from_array(sSSEV))
        oFrame = TransformFrame.createFromVectors(oOrigin, transform.Vec.from_array(oOffsetV), transform.Vec.from_array(oSSEV))
        
        positions = []
        for i in range(len(m1.positions)):
            sPoint = sFrame.transformInto(m1.positions[i]) #we transform from global space to loop space so that we have a relative points... xyz from the SSE, not the origin
            oPoint = oFrame.transformInto(otherPositions[i])
#             positions[i] = Vec({
#                                 'x': sqrt((sPoint.x + oPoint.x) * (sPoint.x + oPoint.x)),
#                                 'y': sqrt((sPoint.y + oPoint.y) * (sPoint.y + oPoint.y)),
#                                 'z': sqrt((sPoint.z + oPoint.z) * (sPoint.z + oPoint.z))
#                                 })

            
            positions.append(Vec({'xyz'[j]: m1.positions[i].__dict__['xyz'[j]] + sqrt((sPoint[j] + oPoint[j]) * (sPoint[j] + oPoint[j])) for j in range(3)}))
        
        return Model([m1,m2], positions, m1.sses, m1.ssesSignature, m1.merge_seqs(m2, m1.size, m2.size), m1.size + m2.size)
        
    @classmethod
    def get_sse_vector(cls, sse_atoms, loop_anchor):
        if(len(sse_atoms) == 1):
            return [loop_anchor.__dict__[c] - sse_atoms[0].__dict__[c] for c in 'xyz']
        else:
            return [sse_atoms[int((len(sse_atoms) * .25))].__dict__[c] - sse_atoms[int((len(sse_atoms).__dict__[c] * .75))] for c in 'xyz']
        
    def score(self, loop):
        """Scores how well the loop matches the Model"""
        
        sys.stderr.write("Score is deprecated, please do not use!\n")
        
        if(len(self.loops) == 0):
            return float("inf")
        
        total = 0
        for l in self.loops:
            total += l.closeness(loop)
        return total / len(self.loops)
    
    def compareOLD(self, other):
        """Scores how close two Models are."""
        #todo
        if (len(other.loops) == 0):
            return float("inf")
        
        total = 0
        for loop in other.loops:
            total += self.score(loop)
        return total / len(other.loops)
    
    def compare(self, other):
        """
        Compares two models to each other. A higher score is worse. Both models MUST have:
        1) The same SSE identifier
        2) The same number of elements in the loop
        """
        #check validity
        if self.ssesSignature != other.ssesSignature:
            return float('inf')
        if len(self.positions) != len(other.positions):
            return float('inf')
        
        if self.ssesSignature[0] == self.ssesSignature[1]: #if ends are helix/helix or sheet/sheet, then you don't know which end aligns with which
            return max(self.__compute_scores(other.sses[0], other.sses[1], other.positions),
                       self.__compute_scores(other.sses[1], other.sses[0], [k for k in reversed(other.positions)])
                       )
        else:
            return self.__compute_scores(other.sses[0], other.sses[1], other.positions)
        
        
        
    def __compute_scores(self, other_sses_0, other_sses_1, other_positions):
        """Private method do not call unless you know what you're doing! Computes how well our model matches up against the given data"""
        #get necessary vectors
        print(self.sses)
        sOffsetV = [self.sses[1][0].__dict__[c] - self.sses[0][0].__dict__[c] for c in 'xyz']
#         sOffsetV = [self.positions[-1].__dict__[c] - self.positions[0].__dict__[c] for c in 'xyz']
#         sOffsetV = 
#         print [self.positions[-1].__dict__[c] for c in 'xyz'], [self.positions[0].__dict__[c] for c in 'xyz']
#         sSSE0V = [self.positions[0].__dict__[c] - mean([atom.__dict__[c] for atom in self.sses[0]]) for c in 'xyz']
#         sSSE1V = [self.positions[-1].__dict__[c] - mean([atom.__dict__[c] for atom in self.sses[1]]) for c in 'xyz']
        sSSE0V = Model.get_sse_vector(self.sses[0], self.positions[0])
        sSSE1V = Model.get_sse_vector(self.sses[1], self.positions[-1]) 

        
        oOffsetV = [other_sses_1[0].__dict__[c] - other_sses_0[0].__dict__[c] for c in 'xyz']
#         oOffsetV = [other_positions[-1].__dict__[c] - other_positions[0].__dict__[c] for c in 'xyz']
#         oSSE0V = [other_positions[0].__dict__[c] - mean([atom.__dict__[c] for atom in other_sses_0]) for c in 'xyz']
#         oSSE1V = [other_positions[-1].__dict__[c] - mean([atom.__dict__[c] for atom in other_sses_1]) for c in 'xyz']
        oSSE0V = Model.get_sse_vector(other_sses_0, other_positions[0])
        oSSE1V = Model.get_sse_vector(other_sses_1, other_positions[-1])
        
        sFrame = TransformFrame.createFromVectors(self.sses[0][0], transform.Vec.from_array(sOffsetV), transform.Vec.from_array(sSSE0V))
        oFrame = TransformFrame.createFromVectors(other_sses_0[0], transform.Vec.from_array(oOffsetV), transform.Vec.from_array(oSSE0V))
        
        total = 0
        for i in range(len(self.positions)):
            sPoint = sFrame.transformInto(self.positions[i]) #we transform from global space to loop space so that we have a relative points... xyz from the SSE, not the origin
            oPoint = oFrame.transformInto(other_positions[i])
            total += rmsd(sPoint, oPoint)
        total /= len(self.positions)
        
        if(total > 2): #must be at most 2 angstroms apart
            return float('inf')
        
        s_theta = arccos(dot(sSSE0V, negative(sOffsetV)) / (norm(sSSE0V) * norm(sOffsetV)))
        s_phi = arccos(dot(sSSE1V, sOffsetV) / (norm(sSSE1V) * norm(sOffsetV)))
        
        o_theta = arccos(dot(oSSE0V, negative(oOffsetV)) / (norm(oSSE0V) * norm(oOffsetV)))
        o_phi = arccos(dot(oSSE1V, oOffsetV) / (norm(oSSE1V) * norm(oOffsetV)))
        
        s_anchor_d = norm(sOffsetV)
        o_anchor_d = norm(oOffsetV)
        
        #Using eculdian distance... is there a better way to do this?
        #I figure anchor d is the biggest value, and that's the one we want weighted the most heavily
        #Perhaps there's some kind of correlation constant or something that compares how close two things are based on how close they are to the average of the two... like variance or something? This will overestimate how similar small things are
        return (s_anchor_d - o_anchor_d) * (s_anchor_d - o_anchor_d) + (s_phi - o_phi) * (s_phi - o_phi) + (s_theta - o_theta) * (s_theta - o_theta) + total;
    
    def get_loops(self, loops):
        """Finds loops from this and all parent sequences and dumps them into loops. Lazy loads loops"""
        if(self.__loops == None):
            l = []
            if(len(self.parents) == 1):
                l = [self.parents[0]]
            else:
                self.parents[0].get_loops(l)
                self.parents[1].get_loops(l)
            self.__loops = l
        loops += self.__loops        
   
    @classmethod
    def gen_seq(cls, loops):
        """helper method to generate a probabilistic sequence representing a cluster of loops (loops must be the same length)"""

        num_loops = len(loops)
        seq = [{} for i in range(len(loops[0]))]
        
        #iterate over each of the loop sequences, adding/updating
        #dictionary entries for each amino acid at each sequence position
        for curr_loop in loops:
            for i in range(len(curr_loop)):
                aa = curr_loop[i] #retrieve the current amino acid
                
                #if the amino acid's already represented at that sequence
                #location, just increment the count
                if seq[i].has_key(aa):
                    seq[i][aa] += 1.0
                
                #otherwise, create a new entry, initialized with a count
                #of 1.0
                else:
                    seq[i][aa] = 1.0
        
        #normalize the sequence (so probabilities sum to 1), then return
        return Model.normalize_seq(seq, num_loops)
    
    @classmethod
    def normalize_seq(cls, seq, n):
        """helper method to normalize a probabilistic sequence, given the total number of entries per sequence location"""
        for pos in seq:
            for aa in pos.keys():
                pos[aa] /= n
        return seq
    
    def merge_seqs(self, other, n1, n2):
        """helper method to merge two probabilistic sequences, given two models (each with their own representative sequence) and their respective weights (i.e. number of loops represented)"""
        merged_seq = [{} for i in range(len(self.seq))]
        weight1 = (n1/float(n1 + n2))
        weight2 = (n2/float(n1 + n2))
        
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
        loops = []
        self.get_loops(loops)
        mean_displacement = self.seq
        loops = []
        self.get_loops(loops)
        return ("Model %s%s" % (mean_displacement, "".join(["\n\t%s" % (loop) for loop in loops])))
