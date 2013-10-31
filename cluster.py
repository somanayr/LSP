'''
Created on Oct 28, 2013
'''
from model import Model
from random import sample
import numpy

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

        