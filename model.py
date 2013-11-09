'''
Created on Oct 28, 2013
'''
import numpy
from numpy.matrixlib.defmatrix import matrix
from numpy.core.fromnumeric import mean
from numpy.core.numeric import cross
from numpy.dual import norm
from transform import TransformFrame

class Model:
    def __init__(self, models=None, loop=None):
        '''
        Generates a model from list of similar loops, loops
        '''
        if models is None and loop is None:
            raise Exception("Must provide models to merge or loop")
        self.loop = loop
        self.models = models
        if models is not None:
            self.sses = None
        else:
            self.sses = (loop.l_anchor, loop.r_anchor)
        self.ssesSignature = models[0].ssesSignature if models is not None else sorted([loop.l_type, loop.r_type])
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
    
    def compareNEW(self, other):
        sOffsetV = {c: self.sses[1][0][c] - self.sses[0][0][c] for c in 'xyz'}
        sSSE0V = {c: mean([atom[c] for atom in self.sses[0]]) for c in 'xyz'}
        
#         n = norm(sOffsetV)
#         sX = [c / n for c in sOffsetV]
#         n = norm(sSSE0V)
#         sY = [c / n for c in sSSE0V]
#         sZ = cross(sX, sY)
#         sY = cross(sX, sZ)
#         sSSE1V = {c: mean([atom[c] for atom in self.sses[1]]) for c in 'xyz'}
        
        oOffsetV = {c: other.sses[1][0][c] - other.sses[0][0][c] for c in 'xyz'}
        oSSE0V = {c: mean([atom[c] for atom in other.sses[0]]) for c in 'xyz'}
        
        
        toFrame = TransformFrame.createFromVectors(self.sses[0][0], sOffsetV, sSSE0V)
        fromFrame = TransformFrame.createFromVectors(other.sses[0][0], oOffsetV, oSSE0V)
        
        for i in range(len(self.positions)):
            point = fromFrame.transformTo(toFrame, other.positions[i])
        
#         n = norm(oOffsetV)
#         oX = [c / n for c in oOffsetV]
#         n = norm(oSSE0V)
#         oY = [c / n for c in oSSE0V]
#         oZ = cross(oX, oY)
#         oY = cross(oX, oZ)
#         oSSE1V = {c: mean([atom[c] for atom in other.sses[1]]) for c in 'xyz'}
        
        
        
        
        
    
    def __str__(self):
        mean_displacement = numpy.mean([loop.displacement() for loop in self.loops])
        return ("Model %f%s" % (mean_displacement, "".join(["\n\t%s" % (loop) for loop in self.loops])))
