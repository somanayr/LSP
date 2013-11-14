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