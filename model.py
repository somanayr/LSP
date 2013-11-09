'''
Created on Oct 28, 2013
'''
import numpy

class Model:
    def __init__(self, models=None, loop=None):
        '''
        Generates a model from list of similar loops, loops
        '''
        if models is None and loop is None:
            raise Exception("Must provide models to merge or loop")
        self.loop = loop
        self.models = models
        self.sses = models[0].sses if models is not None else sorted([(loop.l_type, loop.l_anchor), (loop.r_type, loop.l_anchor)])
        self.phi
        self.theta
        self.positions
        #todo
        
    def score(self, loop):
        """Scores how well the loop matches the Model"""
        #todo
        if(len(self.loops) == 0):
            return float("inf")
        
        total = 0
        for l in self.loops:
            total += l.closeness(loop)
        return total / len(self.loops)
    
    def compare(self, other):
        """Scores how close two Models are."""
        #todo
        if (len(other.loops) == 0):
            return float("inf")
        
        total = 0
        for loop in other.loops:
            total += self.score(loop)
        return total / len(other.loops)
    
    def __str__(self):
        mean_displacement = numpy.mean([loop.displacement() for loop in self.loops])
        return ("Model %f%s" % (mean_displacement, "".join(["\n\t%s" % (loop) for loop in self.loops])))
