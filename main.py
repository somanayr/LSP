from pdb_reader import get_loops
from cluster import length_cluster, sse_similarity_cluster, hierarchical
from model import Model
from classify import classify_loop_seq
from SubstitutionMatrix import blosum62
import os
from validation import *
from time import time

def extract_loops_from_dir(pdb_dir="pdb", loopLimit=-1, fileLimit=50):
    """
    Given a directory, locate PDB files and extract loops from each PDB file.
    @author: Travis Peters
    """    

    ###########################################################################
    # DEBUG FLAGS
    ###########################################################################

    # Indicates whether or not to write extracted loop info to console
    display_file_debug = False

    # Indicates whether or not to write extracted loop info to console
    display_loop_debug = False

    ###########################################################################
    # Loop extraction
    ###########################################################################
    
    # General container for all extracted Loop Structures
    loops = []
    
    # Store a list of files in the PDB directory
    pdb_files = os.listdir(pdb_dir)

    # Keep track of the number of processed PDB files
    files = 0
    
    # Iterate over each file in the PDB directory
    for f in pdb_files:
        
        # If the file is in fact a PDB file, attempt to extract the Loop Structures
        if(f.endswith(".pdb")):
            
            if display_file_debug: 
                print "Loading loops from file: " + f + "...";   
            
            # Get the loop candidates from a PDB file
            try:
                l_cand = get_loops(pdb_dir + "/" + f)
            except:
                continue

            # Check each extracted loop and determine if it is a valid loop (i.e. has left/right anchor)            
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

                # If both the left & right anchor exist and the loop in question is between 7 and 9
                # residues, we will consider this loop to be a valid Loop Structure
                # (i.e. this loop is not an "end point")
                if(len(loop.l_anchor) >= 1 and len(loop.r_anchor) >= 1) and not (len(loop.atoms) < 7 or len(loop.atoms) > 9):
                    loops.append(loop)
            
            # We've read a file in, note that
            files += 1
            
            # Debug: After each file, print the updated amount of loops
            if display_file_debug:
                print " > Total Loops:", len(loops)
            
            # Display progress:
            file_progress = float(files) / fileLimit
            loop_progress = float(len(loops)) / loopLimit
            dir_progress = float(files) / len(pdb_files)
            
            progress = max(file_progress, loop_progress, dir_progress)
            if display_file_debug:
                print("Loading... %.2f%%" % (progress*100))

            # Break if loop limit has been exceeded
            if(loopLimit >= 0 and len(loops) >= loopLimit): break

            # Break if file limit has been exceeded
            if(fileLimit >= 0 and files >= fileLimit): break

    return loops

def create_loop_bins(unordered_loops):
    """
    Given a list of loops, group loops into "bins" based on (1) size, and (2) flanking SSE types
    > NOTE: we currently require loops to be the exact same size if we are to compare them. 
        [YoonJoo told us to consider bins with ranges of lengths (1-3, 4-7, 8-12, etc.)]
    > NOTE: flanking SSEs should be exactly the same when we compare loops.
    @author: Travis Peters
     """
    # Create bins by loop length
    loop_bins_by_length = length_cluster(unordered_loops)

    # Create bins by using loop length info & discriminating between flanking SSEs
    bins = sse_similarity_cluster(loop_bins_by_length)

    return bins


def generate_histogram_data(bin_clusters):
    data = {}
    for bc in bin_clusters:
        for model in bc[1]:
            if not model.size in data:
                data[model.size] = 0
            data[model.size] += 1
            
    for key in sorted(data.keys()):
        print(str(key) + "\t" + str(data[key]))

    
def compute_score_naive(bin_clusters, first_only=True):
    ## Naive testing - check if a loop gets placed into it's model
    total_model_score = 0
    total_structure_score = 0
    total_partial_structure_score = [[0,0],[0,0]]
    total_clusters = len(bin_clusters)
    for bc in bin_clusters:
        bin_data, models = bc
        
        #Find loops
        loop_set = []
        for model in models:
            temp_loops = []
            model.get_loops(temp_loops)
            if(len(temp_loops) > 2): #reject exact and close to exact matches, why test what we know is going to hit 100%?
                loop_set += temp_loops
        
        #Quit if we didn't find any suitable loops
        if(len(loop_set) == 0):
            total_clusters -= 1
            continue
        
        cluster_model_score = 0
        cluster_structure_score = 0
        cluster_partial_structure_score = [[0,0],[0,0]]
        for loop in loop_set:
            scores = classify_loop_seq(loop.seq, models, blosum62)
            scores = sorted(scores, key=lambda x:-x[1])
            model = Model.fromLoop(loop)
            
            #compute model score
            if not first_only:
                structure_score = 0.0
                tries = 0 #start at 1 so no div by zero stuff
                
                
                structure_score = scores[0][0].compare(model, max_rmsd=-1, verbose=False)
                
                #Iterate until we find the match
                for score in scores:
                    temp_loop_set = []
                    score[0].get_loops(temp_loop_set)
                    if loop in temp_loop_set: #is match
                        break 
                    tries+=1
                 
                #Higher score is better!
                model_score = (len(scores) - tries) / (len(scores) + 0.0)
            else:
                #Only search first result
                temp_loop_set = []
                scores[0][0].get_loops(temp_loop_set)
                model_score = 0.0
                structure_score = scores[0][0].compare(model, max_rmsd=-1, verbose=False)
                if loop in temp_loop_set:
                    model_score = 1.0
                    cluster_partial_structure_score[0][0] += structure_score
                    cluster_partial_structure_score[0][1] += 1
                else:
                    cluster_partial_structure_score[1][0] += structure_score
                    cluster_partial_structure_score[1][1] += 1
            
            cluster_model_score += model_score
            cluster_structure_score += structure_score
            
        cluster_model_score /= len(loop_set)
        cluster_structure_score /= len(loop_set)
        
        try:
            cluster_partial_structure_score[0] = cluster_partial_structure_score[0][0] / cluster_partial_structure_score[0][1]
            
            total_partial_structure_score[0][0] += cluster_partial_structure_score[0]
            total_partial_structure_score[0][1] += 1
        except:
            pass
        
        try:    
            cluster_partial_structure_score[1] = cluster_partial_structure_score[1][0] / cluster_partial_structure_score[1][1]
            
            total_partial_structure_score[1][0] += cluster_partial_structure_score[1]
            total_partial_structure_score[1][1] += 1
        except:
            pass
            
        total_model_score += cluster_model_score
        total_structure_score += cluster_structure_score

    if(total_clusters != 0):
        total_model_score /= total_clusters
        total_structure_score /= total_clusters
        try:
            ps1 = total_partial_structure_score[0][0] / total_partial_structure_score[0][1]
        except:
            ps1 = float("nan")
        try:
            ps2 = total_partial_structure_score[1][0] / total_partial_structure_score[1][1]
        except:
            ps2 = float("nan")
        print("Total score: (%f, %f (%f, %f))" % (total_model_score, total_structure_score, ps1, ps2))
    else:
        print("Insufficient data to compute score")
    

if __name__ == '__main__':
    
    ###########################################################################    
    # TODO: Implement user interaction. User should be able to enter
    # information about a loop:
    #   - loop sequence
    #   - left/right SSE anchor types
    # and have information about predicted structure returned. 
    ###########################################################################
        
    ###########################################################################
    # X-Val parameters
    ###########################################################################
#    LOOP_SET = [200, 500, 1000, 2000, 5000]
    LOOP_SET = [200, 500]
    FILES    = -1
    LOOPS    = -1
    FOLDS    = 5
    REPS     = 10

    for l in LOOP_SET:
    
        print "\n[Loop Set:", str(l) + "]"
         
        ###########################################################################
        # Initial loop extraction
        ###########################################################################
        start_time = time()
        loops = extract_loops_from_dir(pdb_dir="pdb", loopLimit=l, fileLimit=FILES)
        print "\nData Load Time:", time() - start_time, "seconds"
    
        print "\nTotal Loops:", len(loops), "\n"
    
        ###########################################################################
        # Run Cross-Validation
        ###########################################################################    
        do_xval_knn(loops, k=1, nfold=FOLDS, nrep=REPS)
