'''
Created on Oct 28, 2013
'''
from pdb_reader import get_loops
from cluster import cluster

if __name__ == '__main__':
    loops = get_loops("112L.pdb")
    print([str(loop) for loop in loops])
    centroids = cluster(loops)
    print(centroids)
