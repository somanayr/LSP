'''
Created on Oct 28, 2013
'''
from pdb_reader import get_loops
#from cluster import cluster, hierarchical
from cluster import length_cluster, hierarchical
from model import Model
import os

def loop_extraction(pdb_dir, maxLoopCount=200):
    """Given a directory, locate PDB files and extract loops from each PDB file."""
    pass

def create_loop_bins(unordered_loops):
    """
    Given a list of loops, group loops into "bins" based on (1) size, and (2) flanking SSE types
    > NOTE: we currently require loops to be the exact same size if we are to compare them. 
        [Yoon Joo told us to consider bins with ranges of lengths (1-3, 4-7, 8-12, etc.)]
    > NOTE: flanking SSEs should be exactly the same when we compare loops.
     """
    # Create length bins
    loop_bins_by_length = length_cluster(unordered_loops)
    
    # TODO: Create sub-bins by discriminating between flanking SSEs
    #loop_bins_by_anchor_types = []

    bins = loop_bins_by_length

    return bins

def classify(in_sequence, predictive_model):
    """
    Given an unknown (loop) sequence and a predictive model, determine which cluster
    best represents the input sequence. 
    """
    pass

if __name__ == '__main__':
#     loops = get_loops("112L.pdb") + get_loops("1A7W.pdb") + get_loops("1B6W.pdb") + get_loops("1QlQ.pdb")
#     print("----LOOPS----")
#     print("\n".join([str(loop) for loop in loops]))
#     centroids = cluster(loops)
#     print("\n----MODELS----")
#     print("\n".join([str(centroid) for centroid in centroids]))

    ###############
    # DEBUG FLAGS #
    ###############
    # Indicates whether or not to write extracted loop info to console
    display_loop_debug = False
    
    # Indicates whether or not to write Model representations of loops to console
    display_loop_model_debug = True
        
    ###########################################################################
    # Initial loop extraction
    ###########################################################################
    
    # General container for all extracted Loop Structures
    loops = []
    
    # Iterate over each file in the PDB directory
    for f in os.listdir("pdb"):
        
        # If the file is in fact a PDB file, attempt to extract the Loop Structures
        if(f.endswith(".pdb")):
            print "Loading loops from file: " + f + "..."
            
            # Get the loop candidates from a PDB file
            l_cand = get_loops("pdb/" + f)

            # Check each extract loop and determine if it is a valid loop (i.e. has left/right anchor)            
            for loop in l_cand:

                # Debug: Write loop info to console
                if display_loop_debug:
                    la = ""
                    if len(loop.l_anchor) > 0:
                        la = "["
                        for a in loop.l_anchor:
                            la += str(a) + ","
                        la = la[0:-1] + "]"
                    else:
                        la = "[]"
    
                    ra = ""
                    if len(loop.r_anchor) > 0:
                        ra = "["
                        for a in loop.r_anchor:
                            ra += str(a) + ","
                        ra = ra[0:-1] + "]"
                    else:
                        ra = "[]"

                    print "[" + la + "," + ra + "]"

                # If both the left & right anchor exist, we will consider this loop to be
                # a valid Loop Structure (i.e. this loop is not an "end point")
                if(len(loop.l_anchor) >= 1 and len(loop.r_anchor) >= 1):
                    loops.append(loop)
            
            # Debug: After each file, print the updated amount of loops
            print " > Total Loops:", len(loops)
            
            # ? (why would we break after finding 200-ish loops?)
            if(len(loops) >= 200): break

    # Debug: Display the Model representation of each Loop Structure extracted from the PDB extraction
#     if display_loop_model_debug:
#         print "\nModel Representations for Loop Structures:"
#         loop_num = 1
#         for loop in loops:
#             model = Model.fromLoop(loop)
#             print loop_num, ":", "<" + str(model.get_loops()) + ">" , model
#             loop_num += 1

# Ryan:        
#    results = hierarchical(loops)
#        
#    print "\n\n\n\n----Results----"
#    for model in results:
#        print model
#        

    ###########################################################################
    # "Bin up" Loop Structures (their Models, now...) by:
    # (1) length --> all loops should be the same size
    # (2) SSE Anchors --> only compare loops that are the same size and that 
    #     have the same left/right anchors
    ###########################################################################
    bins = create_loop_bins(loops)
    
    # Test... (bins are built only by loop size)
    for bin_id in bins:
        print bin_id, bins[bin_id]
    
    ###########################################################################
    # Use hierarchical clustering to form loop similarity clusters
    ###########################################################################
    #clusters = hierarchical(bins)
    
    ###########################################################################
    # Classify - given an input (loop) sequence, match it to some loop cluster
    # or inform user that the loop cannot be characterized
    ###########################################################################
    
    
    ###########################################################################
    # Cross-Validation - input known loop sequences and determine how 
    # successfully we can map the loop sequence back to the correct
    # cluster/model.
    ###########################################################################
    