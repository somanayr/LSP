'''
Created on Oct 28, 2013
'''
from model import Model
from random import sample
from heapq import heappop, heappush
import numpy
import pprint

EPSILON = .0000001

def cluster(loops):
    k = 1
    prev_score = None
    prevCentroids = None
    while True:
        centroids, clusters = kmeans(loops, k)
        score = evaluateClusters(centroids, clusters)
        print(score, prev_score, k)
        print("\n".join([str(centroid) for centroid in centroids]) + "\n" + "\n".join([str([str(c) for c in cluster]) for cluster in clusters]))
        if(prev_score is not None and prev_score - score < EPSILON):
            return prevCentroids
        
        prev_score = score
        prevCentroids = centroids
        
        k += 1
        
def evaluateClusters(centroids, clusters):
    return numpy.mean([numpy.mean([centroids[i].score(loop) for loop in clusters[i]]) for i in range(len(centroids))])
    

def kmeans(loops, k):
    """Clusters loops, an array of Loops, into k clusters
    @author: Serena Liu
    """
    if len(loops) < k:
        raise Exception("Not enough loops to form " + str(k) + " clusters.")
    
    #get k random loops for initial centroids, converting each
    #loop into Model form (for consistency's sake)
    init_loops = sample(loops, k)
    centroids = [Model([loop]) for loop in init_loops]
    
    #keep reclustering and calculating new centroids until
    #the centroid lists (and clusters) converge
    old_score = None
    score = 0
    while (old_score != score):
        clusters = cluster_loops(loops, centroids)
        centroids = get_models(clusters)
        score = evaluateClusters(centroids, clusters)
        old_score = score
    
    #return the final clustering and centroid list
#     clusters = cluster_loops(loops, centroids)
    return (centroids, clusters)
    
    #todo
#     return [loops[0:len(loops)/2], loops[len(loops)/2:-1], loops[-1:]]

def get_models(clusters):
    """Returns a list of representative models for each kmeans
    """
    return [Model(loop = cluster) for cluster in clusters]

def predict(seq, models):
    """Predicts which class the sequence belongs to based on the models
    Returns a class ID corresponding to the index of the best model match
    """ 
    #todo
    return 0

###########################################
"""Helper methods for k-means clustering"""

def find_closest(centroids, loop):
    """Finds the closest match to loop in the list of models, and 
    returns the index of the matching model."""
    min_dist = float("inf")
    min_index = i = 0
    
    #iterate over list, finding distance between each pair
    while i < len(centroids):
        dist = centroids[i].score(loop)
        
        #if new minimum, update min_dist and min_index
        if dist < min_dist:
            min_dist = dist
            min_index = i
        i += 1
    
    return min_index

def cluster_loops(loops, centroids):
    """Completes one iteration of k-means clustering: sorts each loop in
    the list of loops into a centroid kmeans, then returns the resulting
    clustering as a list of lists."""
    
    clusters = [[] for i in range(len(centroids))] #empty list of lists to hold clusters
    
    #find closest centroid for each loop, and add loop to
    #the appropriate kmeans list
    for loop in loops:
        ind = find_closest(centroids, loop)
#         print ind; print [[str(elem) for elem in cluster] for cluster in clusters]; print [str(centroid) for centroid in centroids]; print loop
        clusters[ind].append(loop)
    
    return clusters


##################################################
"""Code for hierarchical clustering"""

#given a list of loops, cluster by length (the number
#of alpha carbons in the loop)
#loops - the initial list of Loops
#returns a dictionary of lists, where each list contains
#loops of the same length
def length_cluster(loops):
    """"Cluster loops by length (number of alpha carbons)."""
    clusters = {}
    
    for loop in loops:
        l = len(loop.atoms)
        
        #if there are no loops of length l in clusters,
        #create a new entry
        if not clusters.has_key(l):
            clusters[l] = [loop]
        
        #otherwise, update existing entry
        else:
            clusters[l].append(loop)
    
    return clusters

def sse_similarity_cluster(loops_by_size):
    """
    Given a cluster of loops (clustered by loop length), return clusters that:
      (1) have the same length, AND
      (2) have the same flanking SSEs.
    
    Note: loops_by_size is a dictionary where the KEY is an integer number representing the
    size of loops in that cluster, and the VALUE is a list of those loops. This function
    will return a a list of tuples with this form:

                            ( loop_length, flank_type, [loops] )
    
    Where loop_length is as described above, flank_type is one of "HH", "HS", "SH", or "SS", and
    the list of loops is in fact a list of the loops belonging to that "bin".
    """
    # Final container of all bins based on length & flanking SSEs
    bins = []
    
    # For each loop length, create bins that discriminate on flanking SSE type
    for bin_id in loops_by_size:

        # Flanking SSE bins "Helix/Helix", "Helix/Sheet", "Sheet/Helix", "Sheet/Sheet"
        HH = []
        HS = []
        SH = []
        SS = []

        # For a given loop length, look at all loops present in the bin and create new bins
        # where we put loops with the same left/right SSE in the same bins        
        for loop in loops_by_size[bin_id]:
            
            if loop.l_type[0] == "H" and loop.r_type[0] == "H":
                HH.append(loop)
            elif loop.l_type[0] == "H" and loop.r_type[0] == "S":
                HS.append(loop)
            elif loop.l_type[0] == "S" and loop.r_type[0] == "H":
                SH.append(loop)
            elif loop.l_type[0] == "S" and loop.r_type[0] == "S":
                SS.append(loop)
            else:
                # This should be an error until we change the implementation
                # to handle loops with None Type flanking SSEs (i.e. no left/right anchor)
                # because it is the first/lasp "loop".
                print "warning: didn't add loop to bin:", loop, "due to flanking SSE abnormality"
        
        # Only create a real bin if there are loops in the bin
        if len(HH) > 0:
            b = (bin_id, "HH", HH)
            bins.append( b )
        if len(HS) > 0:
            b = (bin_id, "HS", HS)
            bins.append( b )
        if len(SH) > 0:
            b = (bin_id, "SH", SH)
            bins.append( b )
        if len(SS) > 0:
            b = (bin_id, "SS", SS)
            bins.append( b )

    return bins

#given a list of loop Models, return a tree (nested lists)
#encoding a hierarchical clustering of the models
#models - the initial list of loop Models
def hierarchical(models):
    """A tree (nested lists) for the loop Models."""
    if not type(models[0]) is Model:
        models = [Model.fromLoop(model) for model in models]
    
    if len(models)<=1: return models #if there's one or fewer models, return as is
                                        #there's no clustering to be done
    
    dist_q = PriQ() #priority queue, sorting by minimum pair distance
    active_subtrees = make_treelist(models) #nothing is merged yet, so all are active
    
    #populate the queue with all unique subtree pairs, keyed by distance apart
    initialize_queue(dist_q, active_subtrees)
    
    #keep merging until there's only one active tree left
    while (len(active_subtrees) > 1):
        min_pair = dist_q.get()
        (sub1, sub2) = min_pair #retrieve both subtrees
        
        #if both subtrees are still active, merge and update queue/lists;
        #otherwise, ignore
        if (sub1 in active_subtrees) and (sub2 in active_subtrees):
            #if the subtree pair is too different to be mergeable, then
            #all other pairs will also be too different (since the PriQ is
            #keyed by distance). In this case, return final cluster list.
            if not (mergeable(sub1, sub2)):
                return [sub[1] for sub in active_subtrees]
            
            else:
                #remove subtrees from active list (they will now be merged)
                active_subtrees.remove(sub1); active_subtrees.remove(sub2)
                
                #merge the subtrees together, then update dist_q and active_subtrees
                #to include the new subtree
                print("Merging %s, %s, %f" % (str(sub1[1]), str(sub2[1]), sub1[1].compare(sub2[1])))
                merged = merge(sub1, sub2) 
                update_queue(dist_q, merged, active_subtrees)
                active_subtrees.append(merged)
    
    #if we arrive at this point, then everything clustered together;
    #just return the overarching representative model in this case
    return [active_subtrees[0][1]]


#helper function to convert a list of Models into a list
#of subtrees
#p_list - the initial list of models
#returns the corresponding list of subtree tuples of the format
#(num_loops, rep_model, subtree)
def make_treelist(m_list):
    tree_list = []
    for i in range(len(m_list)):
        model = m_list[i]
        tree_list.append((1, model,[model])) #only one model, so overall model
                                        #representation and subtree are identical
    return tree_list


#helper function to initialize a priority queue from a list of Models
#by populating it with entries for each unique model pair
#q - the priority queue to be populated
#model_list - the list of subtrees
def initialize_queue(q, model_list):
    #iterate over all unique profile pairs
    for i in range(len(model_list)):
        sub1 = model_list[i]
        for j in range(i+1, len(model_list)):
            sub2 = model_list[j]
            
            #find pair-distance and format as a tuple, then store in queue
            dist = sub1[1].compare(sub2[1]) #score how similar the representative models are
            entry = (sub1, sub2)
            q.put(entry, dist)


#helper function to update a priority queue, with all new unique pairs
#given a new subtree, a list of still active subtrees, and a reference
#to the priority queue to be updated
def update_queue(q, new_tree, active_subtrees):
    new_model = new_tree[1] #retrieve location info for the new subtree
    
    #iterate over each active subtree
    for tree in active_subtrees:
        model = tree[1]
        dist = new_model.compare(model) #find distance to new subtree
        q.put((tree, new_tree), dist) #update the queue with new pair entry


#helper function to determine whether two clusters are similar
#enough to be merged
def mergeable(sub1, sub2):
    """Returns true if two clusters are similar enough to be merged.
    Depends on Model's compare() function returning infinity if the RMSD
    between representative loop structures is above the given threshold."""
    return sub1[1].compare(sub2[1]) != float("inf")


#helper function to merge two subtree entries, given two tuples
#of the format: (num_loops, rep_model, subtree), where num_loops is
#the total number of loops represented by the subtree cluster, 
#rep_model is the representative model for that subtree
#and model is the list of lists representation of the subtree
#sub1, sub2 - the two tuples describing the subtrees
#returns a tuple with the merged subtree information
def merge(sub1, sub2):
    (n1, m1, tree1) = sub1
    (n2, m2, tree2) = sub2
    
    merge_n = n1 + n2
#     tot_loops = m1.loops + m2.loops #masterlist of loops
    rep_model = Model.fromModels(m1, m2)
    merged_tree = [tree1, tree2] #list of lists format
    
    return (merge_n, rep_model, merged_tree)


def print_tree(tree):
    """Print the tree."""
    pprint.PrettyPrinter(indent=2).pprint(tree)


################################################################
"""CBK's priority queue code"""

class QEntry (tuple):
    """A priority queue entry,
heapq library will put the smallest priority value at the root."""
    def __lt__(self, other):
        return self[1] < other[1]

    def __le__(self, other):
        return self[1] <= other[1]
    
    def __gt__(self, other):
        return self[1] > other[1]

    def __ge__(self, other):
        return self[1] >= other[1]

class PriQ:
    """A simple priority queue, with the ability to insert an item with
a particular priority (but not to update the priority) and get the
item with the lowest priority."""
    def __init__(self):
        self.items = []
    
    def __len__(self):
        return len(self.items)
    
    def put(self, item, pri):
        """Insert item into the queue at the given priority."""
        ent = QEntry((item, pri))
        heappush(self.items, ent)
    
    def get(self):
        """Remove and return the item with the lowest priority."""
        return heappop(self.items)[0]
    
    def is_empty(self):
        return len(self) == 0
