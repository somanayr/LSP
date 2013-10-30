'''
Created on Oct 28, 2013
'''
from model import Model
from random import sample

def cluster(loops, k):
    """Clusters loops, an array of Loops, into k clusters
    @author: Serena Liu
    """
    if len(loops) < k:
        raise Exception("Not enough loops to form " + str(k) + " clusters.")
    
    #get k random loops for initial centroids, converting each
    #loop into Model form (for consistency's sake)
    init_loops = sample(loops, k)
    new_centroids = [Model([loop]) for loop in init_loops]
    old_centroids = None
    
    #keep reclustering and calculating new centroids until
    #the centroid lists (and clusters) converge
    while (new_centroids != old_centroids):
        old_centroids = new_centroids
        new_centroids = get_models(cluster_loops(loops, old_centroids))
    
    #return the final clustering and centroid list
    clusters = cluster_loops(loops, new_centroids)
    return (new_centroids, clusters)
    
    #todo
#     return [loops[0:len(loops)/2], loops[len(loops)/2:-1], loops[-1:]]

def get_models(clusters):
    """Returns a list of representative models for each cluster
    """
    return [Model(cluster) for cluster in clusters]

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
    the list of loops into a centroid cluster, then returns the resulting
    clustering as a list of lists."""
    
    clusters = [[]*len(centroids)] #empty list of lists to hold clusters
    
    #find closest centroid for each loop, and add loop to
    #the appropriate cluster list
    for loop in loops:
        ind = find_closest(centroids, loop)
        clusters[ind].append(loop)
    
    return clusters

        