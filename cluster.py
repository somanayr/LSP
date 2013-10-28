'''
Created on Oct 28, 2013
'''

def cluster(loops):
    """Clusters loops, an array of Loops"""
    #todo
    return [loops[0:len(loops)/2], loops[len(loops)/2:-1], loops[-1:]]

def get_models(clusters):
    """Returns a list of representative models for each cluster
    """
    #todo
    return [cluster[0] for cluster in clusters]

def predict(seq, models):
    """Predicts which class the sequence belongs to based on the models
    Returns a class ID corresponding to the index of the best model match
    """ 
    #todo
    return 0