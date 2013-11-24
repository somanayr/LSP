from main import create_loop_bins, extract_loops_from_dir
from cluster import hierarchical
from classify import knn
from time import time
import random

###############################################################################
# Validation Helper Functions
###############################################################################
def make_rand_folds(data, nfold):
    """
    Given a list (data) split it into nfold (more-or-less) even-sized sub-lists.
    @author: Travis Peters
    """

    # Randomize the list so that when we divide it into folds, we have random samples
    data = random.sample(data, len(data))

    avg = len(data) / float(nfold)
    out = []
    last = 0.0
    while last < len(data):
        out.append(data[int(last):int(last + avg)])
        last += avg
    
    return out

def median(l):
    """Compute median of a list."""
    sorted_list = sorted(l)
    length = len(sorted_list)
    if not length % 2:
        return (sorted_list[length / 2] + sorted_list[length / 2 - 1]) / 2.0
    return sorted_list[length / 2]

def mean(l):
    """Compute mean of a list."""
    return sum(l)/float(len(l))

###############################################################################
# Validation
###############################################################################
def xval(loops, classify, nfold=5, nrep=10):
    """
    Cross-validate the classify function on the loops. For each rep (for a total of nrep), 
    split loops into nfold subsets, training on (nfold-1) of them and testing on the left out one. 
    Classify is a function taking a training set and a test set and returns the best model
    and the comparison scores between each instance in the test set & its best matching model.
    @author: Travis Peters
    """
    
    # Hold scores over each rep
    scores = []
    
    for rep in range(nrep):
        print "================================================================"
        print('rep %d\n'%rep)
        rep_start_time = time()
        
        # Split into folds & cross validate the classifier (repeatedly split and predict)    
        folds = make_rand_folds(loops, nfold)

        # Iterate over each fold, holding each one as the "test" fold
        fold_scores = [0.0]*nfold
        fold = 0
        for test_loops in folds:
            print "Fold: ", fold

            ###########################################################################
            # Use the remaining folds as a training set
            ###########################################################################
            training_loops = []
            [training_loops.extend(x) for x in folds if x != test_loops]
            
            ###########################################################################
            # "Bin up" Loop Structures (their Models, now...) by:
            # (1) length --> all loops should be the same size
            # (2) SSE Anchors --> only compare loops that are the same size and that 
            #     have the same left/right anchors
            ###########################################################################
            bins = create_loop_bins(training_loops)

            ###########################################################################
            # Use hierarchical clustering to form loop similarity clusters
            # > Note: loops are clustered within their individual bins. 
            # > Note: training_bin_clusters is a list of tubles where the first element 
            # is the bin tuple and the second element is the clusters returned from 
            # hierarchical clustering (this will be our training set).
            ###########################################################################            
            training_bin_clusters = []
            for bin in bins:
                clusters = hierarchical(bin[2])
                training_bin_clusters.append((bin, clusters))

            ###########################################################################
            # Classify - given an input (loop) sequence, match it to some loop cluster
            ###########################################################################
            classify_results = classify(training_bin_clusters, test_loops)            
            _, cmp_scores = zip(*classify_results)
                        
            ###########################################################################
            # Analyze - record performance 
            ###########################################################################

            # Record median score over each folds for this rep
            fold_scores[fold] = median(cmp_scores) #TODO: compute mean of (medians over each fold)
    
            # Update the processed fold count
            fold += 1

        # Record mean of the median values returned from x-val on the previous folds.        
        mu = mean(fold_scores)
        scores.append(mu)
        
        # Rep. Debug 
        print "Rep.", rep, "Stats:"
        print " > mean =", mu
        print " > execution time:", time() - rep_start_time, "seconds"

    print "================================================================"
    print scores
    print sum(scores)
    print len(scores)
    print sum(scores) / float(len(scores))
    
    print "average score: {:.5f}".format((sum(scores) / float(len(scores))))

def do_xval_knn(loops, k, nfold, nrep):
    """
    Main for the classification cross-validation part.
    @author: cbk - modifications by Travis Peters
    """
    xval(loops, lambda train, test: knn(train, test, k), nfold, nrep)    
    print "params : k=" + str(k) + ", nfold=" + str(nfold) + ", nrep=" + str(nrep) 
    print "================================================================"
    
if __name__ == '__main__':

    ###########################################################################
    # Initial loop extraction
    ###########################################################################
    start_time = time()
    loops = extract_loops_from_dir(pdb_dir="pdb", loopLimit=-1, fileLimit=500)
    print "Data Load Time:", time() - start_time, "seconds"

    ###########################################################################
    # Run Cross-Validation
    ###########################################################################    
    start_time = time()
    do_xval_knn(loops, k=1, nfold=5, nrep=5)
    print "Total Execution Time:", time() - start_time, "seconds"
    