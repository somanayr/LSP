'''
Created on Oct 28, 2013
'''
from numpy.core.numeric import dot
from numpy.dual import norm
from transform import TransformFrame, Vec
from protein_structure import rmsd
from numpy.lib.scimath import arccos
import transform
from numpy.ma.core import negative

class Model:
    
    __phi = [float("inf"), float("-inf")]
    __theta = [float("inf"), float("-inf")]
    __anchor_dist = [float("inf"), float("-inf")]
    
#     def __init__(self, parents, positions, sses, ssesSignature, seq, size):
#         '''
#         Generates a model from list of similar loops, loops
#         '''
#         self.parents = parents
#         self.positions = positions
#         self.sses = sses
#         self.ssesSignature = ssesSignature
#         self.seq = seq
#         self.size = size
#         self.__loops = None #lazy loaded, not safe to call, use get_loops instead
#         

    def __init__(self, parents, positions, theta, phi, anchor_dist, ssesSignature, seq, size):
        self.parents = parents
        self.positions = positions
        self.theta = theta
        self.phi = phi
        self.anchor_dist = anchor_dist
        self.ssesSignature = ssesSignature
        self.seq = seq
        self.size = size
        self.__loops = None

    @classmethod
    def fromLoop(cls, loop):
        """Returns a Model representing the loop"""
        #get necessary vectors
        offset_v = [loop.r_anchor[0].__dict__[c] - loop.l_anchor[0].__dict__[c] for c in 'xyz']
        sse0_v = Model.__get_sse_vector(loop.l_anchor, loop.atoms[0])
        sse1_v = Model.__get_sse_vector(loop.r_anchor, loop.atoms[-1]) 

        
        sFrame = TransformFrame.createFromVectors(loop.l_anchor[0], transform.Vec.from_array(offset_v), transform.Vec.from_array(sse0_v))
        
        #Theta and phi are the angles between the SSE and anchor-anchor vector
        theta = arccos(dot(sse0_v, negative(offset_v)) / (norm(sse0_v) * norm(offset_v)))
        phi = arccos(dot(sse1_v, offset_v) / (norm(sse1_v) * norm(offset_v)))
        
        anchor_dist = norm(offset_v)
        
        update_min_max(phi, Model.__phi)
        update_min_max(theta, Model.__theta)
        update_min_max(anchor_dist, Model.__anchor_dist)
        
        return Model([loop], [Vec.from_array(sFrame.transformInto(atom)) for atom in loop.atoms], theta, phi, anchor_dist, [loop.l_type, loop.r_type], Model.__gen_seq([loop.seq]) , 1)

        
#         sses = sorted([(loop.l_type, loop.l_anchor), (loop.r_type, loop.r_anchor)])
#         return Model([loop], loop.atoms, [loop.l_anchor, loop.r_anchor], [loop.l_type, loop.r_type], Model.__gen_seq([loop.seq]), 1)
        
    @classmethod
    def fromModels(cls, m1, m2):
        """Returns a Model representing a merge of m1 and m2"""
        #check validity
        if m1.ssesSignature != m2.ssesSignature:
            print "SSE signature mistmatch!"
            raise Exception("SSE signature mismatch!")
        if len(m1.positions) != len(m2.positions):
            print "Position size mismatch!"
            raise Exception("Position size mismatch!")
        if abs(m1.anchor_dist - m2.anchor_dist) > 2:
            print "Anchor distance mismatch!"
            raise Exception("Anchor distance mismatch!")
        
        positions = [
                      weighted_average(
                                       m1.positions[pos],
                                       m2.positions[pos],
                                       m1.size,
                                       m2.size
                                       )
                     for pos in range(len(m1.positions))
                     ]
        
        
        print "\nMerging! (%f)" % m1.compare(m2, verbose=True)
        print m1.size, m2.size
        print m1.positions
        print m2.positions
        print positions
        
        
        
        return Model(
                     [m1, m2], #parents
                     positions, #positions
                     weighted_average(m1.theta, m2.theta, m1.size, m2.size), #theta
                     weighted_average(m1.phi, m2.phi, m1.size, m2.size), #phi
                     weighted_average(m1.anchor_dist, m2.anchor_dist, m1.size, m2.size), #anchor dist
                     m1.ssesSignature, #sse signature
                     m1.__merge_seqs(m2, m1.size, m2.size), #seq
                     m1.size + m2.size #size
                     )
        
    @classmethod
    def __get_sse_vector(cls, sse_atoms, loop_anchor):
        """Get the vector along an SSE. If the SSE contains only 1 atom, it will return the vector from that atom to the loop anchor,
        which is the first element in the loop. Otherwise, it will return the vector between the atom at 25% and 75%"""
        if(len(sse_atoms) == 1):
            return [loop_anchor.__dict__[c] - sse_atoms[0].__dict__[c] for c in 'xyz']
        else:
            return [sse_atoms[int((len(sse_atoms) * .25))].__dict__[c] - sse_atoms[int((len(sse_atoms).__dict__[c] * .75))] for c in 'xyz']
    
    def compare(self, other, max_rmsd=2, verbose=False):
        """
        Compares two models to each other. A higher score is worse. Both models MUST have:
        1) The same SSE identifier
        2) The same number of elements in the loop
        """
        #check validity
        if self.ssesSignature != other.ssesSignature:
            print "SSE signature mistmatch!"
            return float('inf')
        if len(self.positions) != len(other.positions):
            print "Position size mismatch!"
            return float('inf')

        return self.__compute_scores(other, max_rmsd=max_rmsd, verbose=verbose)
        
        
        
    def __compute_scores(self, other, max_rmsd=2, verbose=False):
        """Private method do not call unless you know what you're doing! Computes how well our model matches up against the given data"""
        #Average RMSD
        avg_rmsd = 0
        for i in range(len(self.positions)):
            sPoint = [self.positions[i].x, self.positions[i].y, self.positions[i].z] 
            oPoint = [other.positions[i].x, other.positions[i].y, other.positions[i].z]
            avg_rmsd += rmsd(sPoint, oPoint)
        avg_rmsd /= len(self.positions)
        
        avg_rmsd = 0
        
        
        #Cutoffs
        if(max_rmsd > 0 and avg_rmsd > max_rmsd): #must be at most 2 angstroms apart
            if(verbose): print "Structure mismatch! %f > %f" % (avg_rmsd, max_rmsd) 
            return float('inf')
        
        #Compute score
        anchor_dist_score = perc_diff(self.anchor_dist, other.anchor_dist, *Model.__anchor_dist)
        phi_score = perc_diff(self.phi, other.phi, *Model.__phi)
        theta_score = perc_diff(self.theta, other.theta, *Model.__theta)
        
        #Result cutoffs. Must be 95% similarity
        if max_rmsd > 0 and anchor_dist_score > .05:
            if(verbose): print "Anchor distance mismatch! abs(%f-%f) > %f" % (self.anchor_dist, other.anchor_dist, max_rmsd) 
            return float('inf')
         
        if max_rmsd > 0 and phi_score > .05:
            if(verbose): print "Phi mismatch! abs(%f-%f) > %f" % (self.phi, other.phi, max_rmsd) 
            return float('inf')
         
        if max_rmsd > 0 and theta_score > .05:
            if(verbose): print "Theta mismatch! abs(%f-%f) > %f" % (self.theta, other.theta, max_rmsd) 
            return float('inf')
        
        result = anchor_dist_score + phi_score + theta_score
        
#         print "Score %f + %f + %f = %f" % (anchor_dist_score, phi_score, theta_score, result)
        
        return result
        
        #Using eculdian distance... is there a better way to do this?
        #I figure anchor d is the biggest value, and that's the one we want weighted the most heavily
        #Perhaps there's some kind of correlation constant or something that compares how close two things are based on how close they are to the average of the two... like variance or something? This will overestimate how similar small things are
#         return (self.anchor_dist - other.anchor_dist) * (self.anchor_dist - other.anchor_dist) + (self.phi - other.phi) * (self.phi - other.phi) + (self.theta - other.theta) * (self.theta - other.theta) + avg_rmsd;
    
    def get_loops(self, loops):
        """Finds loops from this and all parent sequences and dumps them into loops. Lazy loads loops"""
        if(self.__loops == None): #determine if we have a loop cache
            #if not, generate one recursively
            l = []
            if(len(self.parents) == 1): #base case
                l = [self.parents[0]]
            else: #recursion case
                self.parents[0].get_loops(l)
                self.parents[1].get_loops(l)
                
            #save result for next time so we don't have to recompute
            self.__loops = l
        
        #add our loops to caller's loops. Will not allow caller to modify our cache
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


def weighted_average(a,b,w1,w2):
    """
    Returns a weighted average of a and b with weights w1, w2. w1,w2 do not need to be normalized
    """
    #normalize weights
    total = w1 + w2
    w1 /= float(total)
    w2 /= float(total)
    
    #compute avg
    return a*w1 + b*w2

def perc_diff(a,b,mn,mx):
    """Converts a,b to percentage from mx to mn, then returns the difference."""
    if(a > mx or a < mn):
        raise Exception("a (%s) out of bounds [%s,%s]" % (a, mn, mx))
    if(b > mx or b < mn):
        raise Exception("b (%s) out of bounds [%s,%s]" % (b, mn, mx))
    
    dif = mx-mn
    
    #compute percents
    a_perc = (a - mn) / dif
    b_perc = (b - mn) / dif
    
    return abs(a_perc - b_perc)

def update_min_max(x, min_max_ar):
    """Takes a min/max array [min,max] and a value to update the array with"""
    min_max_ar[0] = min(x, min_max_ar[0])
    min_max_ar[1] = max(x, min_max_ar[1])