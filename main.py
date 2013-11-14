'''
Created on Oct 28, 2013
'''
from pdb_reader import get_loops
from cluster import length_cluster, sse_similarity_cluster, hierarchical
from model import Model
from classify import classify_loop_seq
from SubstitutionMatrix import blosum62
import os

def loop_extraction(pdb_dir, maxLoopCount=200):
    """Given a directory, locate PDB files and extract loops from each PDB file."""
    # TODO: take code from main and make it a function so that main consists of
    # primarily driver code (cleaner, simpler, shorter...).
    pass

def create_loop_bins(unordered_loops):
    """
    Given a list of loops, group loops into "bins" based on (1) size, and (2) flanking SSE types
    > NOTE: we currently require loops to be the exact same size if we are to compare them. 
        [YoonJoo told us to consider bins with ranges of lengths (1-3, 4-7, 8-12, etc.)]
    > NOTE: flanking SSEs should be exactly the same when we compare loops.
     """
    # Create bins by loop length
    loop_bins_by_length = length_cluster(unordered_loops)

    # Debug: (bins are currently only built only by loop size)
#     for bin_id in loop_bins_by_length:
#         print "Bin w/ Loops of Size:", bin_id
#         for loop in loop_bins_by_length[bin_id]:
#             print loop
#     
    # Create bins by using loop length info & discriminating between flanking SSEs
    bins = sse_similarity_cluster(loop_bins_by_length)

    # Debug: (bins are built by both loop size & flanking SSE)
#     s = 0
#     for bin in bins:
#         print "Bin Length:", bin[0], " & Bin Type:", bin[1], " & # Loops:", len(bin[2])
#         #print " > Loops:", [l. for l in bin[2]]
#         s += len(bin[2])
#         
#     print "\nTotal Bins: " + str(len(bins)) + " - Total Loops: " + str(s)
    
    return bins

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
    display_file_debug = True

    # Indicates whether or not to write extracted loop info to console
    display_loop_debug = False

    # Indicates whether or not to write Model representations of loops to console
    display_loop_model_debug = False
        
    ###########################################################################
    # Initial loop extraction
    ###########################################################################
    
    # General container for all extracted Loop Structures
    loops = []
    
    # Iterate over each file in the PDB directory
    for f in os.listdir("pdb"):
        
        # If the file is in fact a PDB file, attempt to extract the Loop Structures
        if(f.endswith(".pdb")):
            
            if display_file_debug: 
                print "Loading loops from file: " + f + "...";   
            
            # Get the loop candidates from a PDB file
            try:
                l_cand = get_loops("pdb/" + f)
            except:
                continue

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
            if display_file_debug: 
                print " > Total Loops:", len(loops)

            # ? (why would we break after finding 200-ish loops?)
#             if(len(loops) >= 200): break

    # Debug: Display the Model representation of each Loop Structure extracted from the PDB extraction
    if display_loop_model_debug:
        print "\nModel Representations for Loop Structures:"
        loop_num = 1
        for loop in loops:
            model = Model.fromLoop(loop)
            print loop_num, ":", "<" + str(model.get_loops()) + ">" , model
            loop_num += 1


    ###########################################################################
    # "Bin up" Loop Structures (their Models, now...) by:
    # (1) length --> all loops should be the same size
    # (2) SSE Anchors --> only compare loops that are the same size and that 
    #     have the same left/right anchors
    ###########################################################################
    bins = create_loop_bins(loops)
    
    ###########################################################################
    # Use hierarchical clustering to form loop similarity clusters
    # Note: loops are clustered within their individual bins. 
    ###########################################################################
    
    # Bin_Clusters is a list of tubles where the first element is the bin tuple
    # and the second element is the clusters returned from hierarchical clustering.
    bin_clusters = []
    for bin in bins:
        #print "Bin:", (bin[0], bin[1], len(bin[2])) 
        clusters = hierarchical(bin[2])
        bin_clusters.append((bin, clusters))

#         # Debug: Display the representative model(s) resulting from clustering
#         print "\n\n\n\n----Results----"
#         for model in clusters:
#             print model
# 
#         _ = raw_input('Press enter to run hierarchical clustering on the next bin...')
    
    ###########################################################################
    # Classify - given an input (loop) sequence, match it to some loop cluster
    # or inform user that the loop cannot be characterized
    ###########################################################################
    test_seq = "NGEMFT"
    print "> Input Sequence:", test_seq, "with length:", len(test_seq)
    
    # First, eliminate bins that don't have sequences with the same length as
    # the input sequence.
    valid_bin_clusters = [bc for bc in bin_clusters if bc[0][0] == len(test_seq)]

    # Now use the valid clusters to attempt to classify the input sequence.
    prediction_scores = []
    for bc in valid_bin_clusters:
        prediction_scores.extend( classify_loop_seq(test_seq, bc[1], subst=blosum62) )

    ## Display results #################################################################
    sorted_prediction_scores = prediction_scores.sort(key=lambda x: x[1], reverse=True)
    TOP_SCORES = 5
    for modelNum, model in enumerate(prediction_scores):
        print "Model Score:", model[1], "\n", model[0]
        
        # Only display a total of TOP_SCORES of the best models/scores
        if modelNum >= TOP_SCORES-1: break
        
    ####################################################################################
    
    ##Naive testing - check if a loop gets placed into it's model
    total_score = 0
    total_clusters = len(bin_clusters)
    for bc in bin_clusters:
        bin_data, models = bc
        
        #Find loops
        loop_set = []
        for model in models:
            temp_loops = []
            model.get_loops(temp_loops)
            if(len(temp_loops) > 2): #reject exact and close to exact matches, why test what we know is going to hit 100%?
#             if(True):
                loop_set += temp_loops
        
        #Quit if we didn't find any suitable loops
        if(len(loop_set) == 0):
            total_clusters -= 1
            continue
        
        cluster_score = 0
        for loop in loop_set:
            scores = classify_loop_seq(loop.seq, models, blosum62)
            scores = sorted(scores, key=lambda x:-x[1])
            tries = 1 #start at 1 so no div by zero stuff
            
            #Iterate until we find the match
            for score in scores:
                temp_loop_set = []
                score[0].get_loops(temp_loop_set)
                if loop in temp_loop_set: #is match
#                     print "Success after %d tries" % (tries)
                    break
                tries+=1
            
            #Higher score is better!
            cluster_score += 1.0/tries
        cluster_score /= len(loop_set)
        print "Cluster score on bin %s, %d: %f (bin size models=%d, loops=%d)" % (str(bc[0][1]), len(bc[0][2][0].seq), cluster_score, len(models), len(loop_set))
        
        total_score += cluster_score
    
    total_score /= total_clusters
    print("Total score: %f" % total_score)
   
    ###########################################################################
    # Cross-Validation - input known loop sequences and determine how 
    # successfully we can map the loop sequence back to the correct
    # cluster/model.
    ###########################################################################
    