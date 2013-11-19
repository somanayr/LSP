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
        self.__loops = None #lazy loaded, not safe to call, use get_loops instead
        
    @classmethod
    def fromLoop(cls, loop):
        """Returns a Model representing the loop"""
#         sses = sorted([(loop.l_type, loop.l_anchor), (loop.r_type, loop.r_anchor)])
        return Model([loop], loop.atoms, [loop.l_anchor, loop.r_anchor], [loop.l_type, loop.r_type], Model.__gen_seq([loop.seq]), 1)
        
    @classmethod
    def fromModels(cls, m1, m2):
        """Returns a Model representing a merge of m1 and m2"""
        if m1.ssesSignature != m2.ssesSignature:
            return Exception("Loop SSE signature mismatch!")
        if len(m1.positions) != len(m2.positions):
            raise Exception("Loop length mismatch!")

        #find necessary vectors
        sOffsetV = [m1.sses[1][0].__dict__[c] - m1.sses[0][0].__dict__[c] for c in 'xyz']
        sSSEV = Model.__get_sse_vector(m1.sses[0], m1.positions[0])
        
        oOffsetV = [m2.sses[1][0].__dict__[c] - m2.sses[0][0].__dict__[c] for c in 'xyz']
        oSSE0V = Model.__get_sse_vector(m1.sses[0], m1.positions[0])
#         oSSE1V = Model.__get_sse_vector(m1.sses[1], m1.positions[-1])
        
        oSSEV = oSSE0V
        
        sOrigin = m1.sses[0][0]
        oOrigin = m2.sses[0][0]
        
        otherPositions = m2.positions
        
#         if m1.ssesSignature[0] == m1.ssesSignature[1]: #if we have helix/helix or sheet/sheet, we don't know which sheet is "first", so decide that based on which has a better score
#             scoreNormal = m1.__compute_scores(m2.sses[0], m2.sses[1], m2.positions)
#             revPos = [k for k in reversed(m2.positions)]
#             scoreReversed = m1.__compute_scores(m2.sses[1], m2.sses[0], revPos)
#             
#             #if we get a better score with the reversed orientation, flip the orientation
#             if(scoreReversed < scoreNormal):
#                 otherPositions = revPos
#                 oSSEV = oSSE1V
#                 oOffsetV = negative(oOffsetV)
#                 oOrigin = m2.sses[1][0]
#                 
#             revPos = None #clear the memory
        
        sFrame = TransformFrame.createFromVectors(sOrigin, transform.Vec.from_array(sOffsetV), transform.Vec.from_array(sSSEV))
        oFrame = TransformFrame.createFromVectors(oOrigin, transform.Vec.from_array(oOffsetV), transform.Vec.from_array(oSSEV))
        
        #Determine the position
        positions = []
        for i in range(len(m1.positions)):
            #transform
            sPoint = sFrame.transformInto(m1.positions[i]) #we transform from global space to loop space so that we have a relative points... xyz from the SSE, not the origin
            oPoint = oFrame.transformInto(otherPositions[i])
            
            #Position is m1.positions[i] - sPoint + (sPoint + oPoint) / 2, in other words, the average, displaced back to the position of m1
            v = Vec.from_array(sFrame.transformOutOf(Vec({'xyz'[j]: (sPoint[j] + oPoint[j]) / 2 for j in range(3)})))
            
            positions.append(v)
        
        return Model([m1,m2], positions, m1.sses, m1.ssesSignature, m1.__merge_seqs(m2, m1.size, m2.size), m1.size + m2.size)
        
    @classmethod
    def __get_sse_vector(cls, sse_atoms, loop_anchor):
        """Get the vector along an SSE. If the SSE contains only 1 atom, it will return the vector from that atom to the loop anchor,
        which is the first element in the loop. Otherwise, it will return the vector between the atom at 25% and 75%"""
        if(len(sse_atoms) == 1):
            return [loop_anchor.__dict__[c] - sse_atoms[0].__dict__[c] for c in 'xyz']
        else:
            return [sse_atoms[int((len(sse_atoms) * .25))].__dict__[c] - sse_atoms[int((len(sse_atoms).__dict__[c] * .75))] for c in 'xyz']
    
    def compare(self, other, max_rmsd=2):
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
        
#         if self.ssesSignature[0] == self.ssesSignature[1]: #if ends are helix/helix or sheet/sheet, then you don't know which end aligns with which
#             #Try both, return the best
#             return max(self.__compute_scores(other.sses[0], other.sses[1], other.positions, max_rmsd),
#                        self.__compute_scores(other.sses[1], other.sses[0], [k for k in reversed(other.positions)], max_rmsd))
#         else:


        return self.__compute_scores(other.sses[0], other.sses[1], other.positions)
        
        
        
    def __compute_scores(self, other_sses_0, other_sses_1, other_positions, max_rmsd=2):
        """Private method do not call unless you know what you're doing! Computes how well our model matches up against the given data"""
        #get necessary vectors
        sOffsetV = [self.sses[1][0].__dict__[c] - self.sses[0][0].__dict__[c] for c in 'xyz']
        sSSE0V = Model.__get_sse_vector(self.sses[0], self.positions[0])
        sSSE1V = Model.__get_sse_vector(self.sses[1], self.positions[-1]) 

        
        oOffsetV = [other_sses_1[0].__dict__[c] - other_sses_0[0].__dict__[c] for c in 'xyz']
        oSSE0V = Model.__get_sse_vector(other_sses_0, other_positions[0])
        oSSE1V = Model.__get_sse_vector(other_sses_1, other_positions[-1])
        
        sFrame = TransformFrame.createFromVectors(self.sses[0][0], transform.Vec.from_array(sOffsetV), transform.Vec.from_array(sSSE0V))
        oFrame = TransformFrame.createFromVectors(other_sses_0[0], transform.Vec.from_array(oOffsetV), transform.Vec.from_array(oSSE0V))
        
        #Average RMSD
        avg_rmsd = 0
        for i in range(len(self.positions)):
            sPoint = sFrame.transformInto(self.positions[i]) #we transform from global space to loop space so that we have a relative points... xyz from the SSE, not the origin
            oPoint = oFrame.transformInto(other_positions[i])
            avg_rmsd += rmsd(sPoint, oPoint)
        avg_rmsd /= len(self.positions)
        
        if(max_rmsd > 0 and avg_rmsd > max_rmsd): #must be at most 2 angstroms apart
            return float('inf')
        
        #Theta and phi are the angles between the SSE and anchor-anchor vector
        s_theta = arccos(dot(sSSE0V, negative(sOffsetV)) / (norm(sSSE0V) * norm(sOffsetV)))
        s_phi = arccos(dot(sSSE1V, sOffsetV) / (norm(sSSE1V) * norm(sOffsetV)))
        
        o_theta = arccos(dot(oSSE0V, negative(oOffsetV)) / (norm(oSSE0V) * norm(oOffsetV)))
        o_phi = arccos(dot(oSSE1V, oOffsetV) / (norm(oSSE1V) * norm(oOffsetV)))
        
        s_anchor_d = norm(sOffsetV)
        o_anchor_d = norm(oOffsetV)
        
        #Using eculdian distance... is there a better way to do this?
        #I figure anchor d is the biggest value, and that's the one we want weighted the most heavily
        #Perhaps there's some kind of correlation constant or something that compares how close two things are based on how close they are to the average of the two... like variance or something? This will overestimate how similar small things are
        return (s_anchor_d - o_anchor_d) * (s_anchor_d - o_anchor_d) + (s_phi - o_phi) * (s_phi - o_phi) + (s_theta - o_theta) * (s_theta - o_theta) + avg_rmsd;
    
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
    def __gen_seq(cls, loops):
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
        return Model.__normalize_seq(seq, num_loops)
    
    @classmethod
    def __normalize_seq(cls, seq, n):
        """helper method to normalize a probabilistic sequence, given the total number of entries per sequence location"""
        for pos in seq:
            for aa in pos.keys():
                pos[aa] /= n
        return seq
    
    def __merge_seqs(self, other, n1, n2):
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
        
        return ("Model[\n\tsize=%d,\n\tsignature=%s,\n\tsses=[\n\t\t%s,\n\t\t%s\n\t]\n\tloop_len=%d,\n\tseq=%s,\n\tpositions=%s\n]" % (
                                                                                        self.size,
                                                                                        str(self.ssesSignature),
                                                                                        str([Vec(sses) for sses in self.sses[0]]),
                                                                                        str([Vec(sses) for sses in self.sses[1]]),
                                                                                        len(self.positions),
                                                                                        "[\n\t\t" + "\n\t\t".join([str(elem) for elem in self.seq]) + "\n\t]",
                                                                                        "[\n\t\t" + "\n\t\t".join([str(pos if type(pos) is Vec else Vec(pos)) for pos in self.positions]) + "\n\t]"
                                                                                        )
                )
