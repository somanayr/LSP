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

class Model:
    def __init__(self, parents, positions, sses, ssesSignature):
        '''
        Generates a model from list of similar loops, loops
        '''
        self.parents = parents
        self.positions = positions
        self.sses = sses
        self.ssesSignature = ssesSignature
        #todo
        
    @classmethod
    def fromLoop(cls, loop):
        """Returns a Model representing the loop"""
        sses = sorted([(loop.l_type, loop.l_anchor), (loop.r_type, loop.r_anchor)])
        return Model([loop], loop.atoms, [sses[0][1], sses[1][1]], "".join([sses[0][0], sses[1][0]]));
        
    @classmethod
    def fromModels(cls, m1, m2):
        """Returns a Model representing a merge of m1 and m2"""
        #find necessary vectors
        sOffsetV = [m1.sses[1][0][c] - m1.sses[0][0][c] for c in 'xyz']
        sSSEV = [mean([atom[c] for atom in m1.sses[0]]) for c in 'xyz']

        oOffsetV = [m2.sses[1][0][c] - m2.sses[0][0][c] for c in 'xyz']
        oSSE0V = [mean([atom[c] for atom in m2.sses[0][0]]) for c in 'xyz']
        oSSE1V = [mean([atom[c] for atom in m2.sses[1]]) for c in 'xyz']
        
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
                oOffsetV = -oOffsetV
                oOrigin = m2.sses[1][0]
            revPos = None
        
        sFrame = TransformFrame.createFromVectors(sOrigin, transform.Vec.from_array(sOffsetV), transform.Vec.from_array(sSSEV))
        oFrame = TransformFrame.createFromVectors(oOrigin, transform.Vec.from_array(oOffsetV), transform.Vec.from_array(oSSEV))
        
        positions = [0]*len(m1.positions)
        for i in range(len(m1.positions)):
            sPoint = sFrame.transformInto(m1.positions[i]) #we transform from global space to loop space so that we have a relative points... xyz from the SSE, not the origin
            oPoint = oFrame.transformInto(otherPositions[i])
            positions[i] = Vec({
                                'x': sqrt((sPoint.x + oPoint.x) * (sPoint.x + oPoint.x)),
                                'y': sqrt((sPoint.y + oPoint.y) * (sPoint.y + oPoint.y)),
                                'z': sqrt((sPoint.z + oPoint.z) * (sPoint.z + oPoint.z))
                                })
        
        return Model([m1,m2], positions, m1.sses, m1.ssesSignature)
        
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
        """
        Compares two models to each other. A higher score is worse. Both models MUST have:
        1) The same SSE identifier
        2) The same number of elements in the loop
        """
        #check validity
        if self.ssesSignature != other.ssesSignature:
            return float('-inf')
        if len(self.positions) != len(other.positions):
            return float('-inf')
        
        if self.ssesSignature[0] == self.ssesSignature[1]: #if ends are helix/helix or sheet/sheet, then you don't know which end aligns with which
            return max(self.__compute_scores(other.sses[0], other.sses[1], other.positions),
                       self.__compute_scores(other.sses[1], other.sses[0], [k for k in reversed(other.positions)])
                       )
        else:
            return self.__compute_scores(other.sses[0], other.sses[1], other.positions)
        
        
        
    def __compute_scores(self, other_sses_0, other_sses_1, other_positions):
        """Private method do not call unless you know what you're doing! Computes how well our model matches up against the given data"""
        #get necessary vectors
        sOffsetV = [self.sses[1][0][c] - self.sses[0][0][c] for c in 'xyz']
        sSSE0V = [mean([atom[c] for atom in self.sses[0]]) for c in 'xyz']
        sSSE1V = [mean([atom[c] for atom in self.sses[1]]) for c in 'xyz']

        
        oOffsetV = [other_sses_1[0][c] - other_sses_0[0][c] for c in 'xyz']
        oSSE0V = [mean([atom[c] for atom in other_sses_0[0]]) for c in 'xyz']
        oSSE1V = [mean([atom[c] for atom in other_sses_1]) for c in 'xyz']
        
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
        
        s_theta = arccos(dot(sSSE0V, -sOffsetV) / (norm(sSSE0V) * norm(sOffsetV)))
        s_phi = arccos(dot(sSSE1V, sOffsetV) / (norm(sSSE1V) * norm(sOffsetV)))
        
        o_theta = arccos(dot(sSSE0V, -sOffsetV) / (norm(oSSE0V) * norm(oOffsetV)))
        o_phi = arccos(dot(oSSE1V, oOffsetV) / (norm(oSSE1V) * norm(oOffsetV)))
        
        s_anchor_d = norm(sOffsetV)
        o_anchor_d = norm(oOffsetV)
        
        #Using eculdian distance... is there a better way to do this?
        #I figure anchor d is the biggest value, and that's the one we want weighted the most heavily
        #Perhaps there's some kind of correlation constant or something that compares how close two things are based on how close they are to the average of the two... like variance or something? This will overestimate how similar small things are
        return (s_anchor_d - o_anchor_d) * (s_anchor_d - o_anchor_d) + (s_phi - o_phi) * (s_phi - o_phi) + (s_theta - o_theta) * (s_theta - o_theta) + total;
    
    def gen_seq(self):
        """helper method to generate a probabilistic sequence representing a cluster of loops (loops must be the same length)"""
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
    
    def normalize_seq(self, seq, n):
        """helper method to normalize a probabilistic sequence, given the total number of entries per sequence location"""
        for pos in seq:
            for aa in pos.keys():
                pos[aa] /= n
    
    def merge_seqs(self, other, n1, n2):
        """helper method to merge two probabilistic sequences, given two models (each with their own representative sequence) and their respective weights (i.e. number of loops represented)"""
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
