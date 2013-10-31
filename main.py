'''
Created on Oct 28, 2013
'''
from pdb_reader import get_loops
from cluster import cluster

if __name__ == '__main__':
    loops = get_loops("112L.pdb") + get_loops("1A7W.pdb") + get_loops("1B6W.pdb") + get_loops("1QlQ.pdb")
    print("----LOOPS----")
    print("\n".join([str(loop) for loop in loops]))
    centroids = cluster(loops)
    print("\n----MODELS----")
    print("\n".join([str(centroid) for centroid in centroids]))
