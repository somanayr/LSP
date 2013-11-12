'''
Created on Oct 28, 2013
'''
from pdb_reader import get_loops
from cluster import cluster, hierarchical
from model import Model
import os

if __name__ == '__main__':
#     loops = get_loops("112L.pdb") + get_loops("1A7W.pdb") + get_loops("1B6W.pdb") + get_loops("1QlQ.pdb")
#     print("----LOOPS----")
#     print("\n".join([str(loop) for loop in loops]))
#     centroids = cluster(loops)
#     print("\n----MODELS----")
#     print("\n".join([str(centroid) for centroid in centroids]))

    loops = []
    for f in os.listdir("pdb"):
        if(f.endswith(".pdb")):
            print "Loading loops from file: " + f
            l_cand = get_loops("pdb/" + f) #get loop candidates
            
            for loop in l_cand:
                print loop.l_anchor
                print loop.r_anchor
                if(len(loop.l_anchor) >= 1 and len(loop.r_anchor) >= 1):
                    loops.append(loop)
            
            print len(loops)
            
            if(len(loops) >= 200): break
    for loop in loops:
        print Model.fromLoop(loop)
        
    print hierarchical(loops)
        
    
