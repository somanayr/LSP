'''
Created on Nov 9, 2013

@author: ramos
'''
import numpy

class Vec:
    def __init__(self, dict):
        self.x = dict['x']
        self.y = dict['y']
        self.z = dict['z']
        
    def __str__(self):
        return ("(%s,%s,%s)" % (str(self.x), str(self.y), str(self.z)))

class TransformFrame:
    def __init__(self, o, x, y, z):
        self.o = o
        self.x = x
        self.y = y
        self.z = z
    
    def transformInto(self, point):
        a = [point.x * c for c in self.x]
        b = [point.y * c for c in self.y]
        c = [point.z * c for c in self.z]
        
        return Vec({'xyz'[i]: self.o[i] + a[i] + b[i] + c[i] for i in range(3)})
    
    def transformOutOf(self, point):
        translated = [point.__dict__['xyz'[i]] - self.o[i] for i in range(3)]
        return Vec({c: numpy.dot(translated, self.__dict__[c]) for c in 'xyz'})
    
    def transformTo(self, other, point):
        return other.transformInto(self.transformOutOf(point))
    
    @classmethod
    def createFromVectors(cls, origin, x, y):
        origin = [origin.__dict__[c] for c in 'xyz']
        x = [x.__dict__[c] for c in 'xyz']
        y = [y.__dict__[c] for c in 'xyz']
        norm = numpy.linalg.norm(x)
        x = [c / norm for c in x]
        norm = numpy.linalg.norm(y)
        y = [c / norm for c in y]
        z = numpy.cross(x, y)
        y = numpy.cross(x, z)
        return TransformFrame(origin,x,y,z)
    
    
if __name__ == "__main__":
    xVec = Vec({
             "x": 1,
             "y": 1,
             "z": 0
             })
    yVec = Vec({
             "x": -1,
             "y": 1,
             "z": 0
             })
#     zVec = {
#             "x": 0,
#             "y": 0,
#             "z": 1
#             }
    origin = Vec({
              "x": 0,
              "y": 0,
              "z": 0
              })
    
    point = Vec({
             "x": 1,
             "y": 1,
             "z": 1
             })
    
    frame = TransformFrame.createFromVectors(origin, xVec, yVec)
    print frame.transformOutOf(point)
