from model import Model
import random
from SubstitutionMatrix import blosum62

def classify_loop_seq(input_seq, clusters, subst=blosum62):
    """
    Classify - given an input (loop) sequence and a list of representative loop
    models, match it to some loop cluster or inform user that the loop cannot be characterized. 
    This essentially assigns the input sequence to a cluster and suggests that the loop
    structure of the cluster characterizes the input sequence. 
    
    TODO: we haven't handled the case where we should not classify the loop_sequence.
    TODO: once we have identified a representative cluster for an input sequence, we
    should actually say something more meaningful about its predicted structure. 
    @author: Travis Peters
    """
    scores = []
    for cluster in clusters:

        # Compute score for each incident of an AA at a particular position. The score
        # is calculated by computing the substitution score for the input sequence and
        # ever incident AA at that position in the model sequence. We take the sum of
        # all of these scores over every position in the input sequence and that
        # determines the matching score for the input sequence and some model (cluster). 
        score = 0.0
        
        # Iterate over each element of the sequence
        for index, seq_element in enumerate(cluster.seq):
            
            # Iterate over each incident AA at some position in the sequence
            for incident in seq_element:
                score += subst[input_seq[index], incident] * seq_element[incident]
        
        scores.append( (cluster, score) )
    
    return scores # a list of tuples (model, match_score)

def rand_neighbors(instance, others, k=1):
    """
    Return k random neighbors to the test instance
    @author: Travis Peters
    """
    
    # Eliminate bins where sequence length != to length of the test instance's sequence
    valid_bin_clusters = [bc for bc in others if ((bc[0][0] == len(instance.seq)) and (bc[0][1] == instance.getSSESignature()))]

    # Now use the valid clusters to attempt to classify the test instance sequence.
    prediction_scores = []
    for bc in valid_bin_clusters:
        prediction_scores.extend( classify_loop_seq(instance.seq, bc[1]) )

    # Sort the "neighbors" by their computed classification score.
    prediction_scores.sort(key=lambda x: x[1], reverse=True)

    # Now that we have a list of neighbors from closest to farthest
    # (best scoring to worst scoring), drop the distance data and just 
    # construct a list of the neighbors themselves. 
    neighbors = [x[0] for x in prediction_scores]
    
    rand_neighbors = random.sample(neighbors, len(neighbors))

    return rand_neighbors[:k]

def neighbors(instance, others, k=1):
    """
    Return the k nearest neighbors to instance among others.
    
    Inputs:
    # instance - a given loop.
    # others   - a list of ALL the other loops/models.
    # k        - the number of closest neighbors to return.

    Output:
    # list of the k nearest neighbors.
    @author: Travis Peters
    """

    # Eliminate bins where sequence length != to length of the test instance's sequence
    valid_bin_clusters = [bc for bc in others if ((bc[0][0] == len(instance.seq)) and (bc[0][1] == instance.getSSESignature()))]

    # Now use the valid clusters to attempt to classify the test instance sequence.
    prediction_scores = []
    for bc in valid_bin_clusters:
        prediction_scores.extend( classify_loop_seq(instance.seq, bc[1]) )

    # Sort the "neighbors" by their computed classification score.
    prediction_scores.sort(key=lambda x: x[1], reverse=True)

    # Now that we have a list of neighbors from closest to farthest
    # (best scoring to worst scoring), drop the distance data and just 
    # construct a list of the neighbors themselves. 
    neighbors = [x[0] for x in prediction_scores]

    return neighbors[:k]


def knn(train, test, k=1, randFlag=False):
    """
    Return a list of comparison scores computed between each loop in test and its 
    'nearest neighbor' (highest scoring model) from the models in train.
    @author: Travis Peters
    """
    
    # Given a test instance, find its k-nearest neighbors in the training set
    knn_list = []
    for test_instance in test:
        if not randFlag:
            knn_list.append( (test_instance, neighbors(test_instance, train, k)) )
        else:
            knn_list.append( (test_instance, rand_neighbors(test_instance, train, k)) )

    # Compare each test instance and its best matching model...
    test_results = []
    for t in knn_list:
        # Un-pack test/best matching loop
        testLoop = t[0]
        bestMatchLoop = t[1][0]
        
        # Compute comparison
        testLoopModel = Model.fromLoop(testLoop)
        cmp_score = testLoopModel.compare(bestMatchLoop, max_rmsd=-1)        

        # Record (1) best matching model, and (2) comparison score
        test_results.append( (t[1], cmp_score) )

    # Tuple of ( best_model, compare_score)
    return test_results