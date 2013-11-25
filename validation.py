from main import *
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
def xval(loops, classify, nfold=5, nrep=10, perc_cutoff=.05):
    """
    Cross-validate the classify function on the loops. For each rep (for a total of nrep), 
    split loops into nfold subsets, training on (nfold-1) of them and testing on the left out one. 
    Classify is a function taking a training set and a test set and returns the best model
    and the comparison scores between each instance in the test set & its best matching model.
    @author: Travis Peters
    """
    
    # Hold scores over each rep
    median_scores = []
    scores = []

    # Baseline testing
    rand_median_scores = []
    rand_scores = []
    
    for rep in range(nrep):
        #print "================================================================"
        #print('rep %d\n'%rep)
        #rep_start_time = time()
        
        # Split into folds & cross validate the classifier (repeatedly split and predict)    
        folds = make_rand_folds(loops, nfold)

        # Iterate over each fold, holding each one as the "test" fold
        median_fold_scores = [0.0]*nfold
        mean_fold_scores = [0.0]*nfold 
        
        # Baseline testing
        rand_median_fold_scores = [0.0]*nfold
        rand_mean_fold_scores = [0.0]*nfold 
               
        fold = 0
        for test_loops in folds:
            #print "Fold: ", fold

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
                clusters = hierarchical(bin[2], perc_cutoff)
                training_bin_clusters.append((bin, clusters))

            ###########################################################################
            # Classify - given an input (loop) sequence, match it to some loop cluster
            ###########################################################################
            classify_results = classify(training_bin_clusters, test_loops, randFlag=False)            
            _, cmp_scores = zip(*classify_results)

            # Baseline testing - random selection
            rand_classify_results = classify(training_bin_clusters, test_loops, randFlag=True)            
            _, rand_cmp_scores = zip(*rand_classify_results)

            ###########################################################################
            # Analyze - record performance 
            ###########################################################################

            # Record median/mean score over each folds for this rep
            median_fold_scores[fold] = median(cmp_scores)
            mean_fold_scores[fold] = mean(cmp_scores)

            # Record median/mean for baseline testing
            rand_median_fold_scores[fold] = median(rand_cmp_scores)
            rand_mean_fold_scores[fold] = mean(rand_cmp_scores) 
                
            # Update the processed fold count
            fold += 1

        # Record mean of the median values & true mean returned from x-val on the previous folds.        
        median_mu = mean(median_fold_scores)
        median_scores.append(median_mu)
        
        mu = mean(mean_fold_scores)        
        scores.append(mu)

        # Record mean of the median values & true mean returned from baseline x-val on the previous folds.        
        rand_median_mu = mean(rand_median_fold_scores)
        rand_median_scores.append(rand_median_mu)
        
        rand_mu = mean(rand_mean_fold_scores)        
        rand_scores.append(rand_mu)
        
        # Rep. Debug #############################################################
        #print "Rep.", rep, "Stats:"
        #print " > mean of medians          =", median_mu
        #print " > true mean                =", mu
        #print " > [random] mean of medians =", rand_median_mu
        #print " > [random] true mean       =", rand_mu
        #print " > execution time           =", time() - rep_start_time, "seconds"
        ##########################################################################

    print "================================================================"
    print "average median score         : {:.5f}".format((sum(median_scores) / float(len(median_scores))))
    print "average score                : {:.5f}".format((sum(scores) / float(len(scores))))
    print "[random] average median score: {:.5f}".format((sum(rand_median_scores) / float(len(rand_median_scores))))
    print "[random] average score       : {:.5f}".format((sum(rand_scores) / float(len(rand_scores))))

def do_xval_knn(loops, k=1, nfold=5, nrep=50, perc_cutoff=.02):
    """
    Main for the classification cross-validation part.
    @author: cbk - modifications by Travis Peters
    """
    start_time = time()
    xval(loops, lambda train, test, randFlag: knn(train, test, k, randFlag), nfold, nrep, perc_cutoff)    
    print "params : k=" + str(k) + ", nfold=" + str(nfold) + ", nrep=" + str(nrep) + ", perc_cutoff=" + str(perc_cutoff) 
    print "================================================================"
    print "\nTotal Execution Time:", time() - start_time, "seconds"

    
if __name__ == '__main__':
    
    ###########################################################################
    # X-Val parameters
    ###########################################################################
    LOOP_SET = [200, 500, 1000, 2000, 5000]
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
                