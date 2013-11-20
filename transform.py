'''
Created on Nov 9, 2013

@author: Ryan Amos
'''
import numpy
import math

class Vec:
    """Represents a 3d vector. Has components x,yz"""
    def __init__(self, d):
        """Creates a Vec from a dictionary of {'x': x, 'y': y, 'z': z}"""
        if(not isinstance(d, dict)):
            d = d.__dict__
        self.x = d['x']
        self.y = d['y']
        self.z = d['z']
        
    def __str__(self):
        return ("(%s,%s,%s)" % (str(self.x), str(self.y), str(self.z)))
    
    @classmethod
    def from_array(cls, ar):
        """Generates a Vector from an array [x,y,z]"""
        dict = {}
        dict["x"] = ar[0]
        dict["y"] = ar[1]
        dict["z"] = ar[2]
        return Vec(dict)
    
    def __repr__(self):
        return self.__str__()

RET_MODE_VECTOR = 0
RET_MODE_ARRAY = 1

class TransformFrame:
    
    def __init__(self, o, x, y, z):
        self.o = o
        self.x = x
        self.y = y
        self.z = z
        self.retMode = RET_MODE_ARRAY
    
    def transformOutOf(self, point):
        """Algorithm taken from Dartmouth class CS77 framework code by Jonathan Denning, http://www.cs.dartmouth.edu/~cs77/assignments/assignment02.zip"""
        a = [point.x * c for c in self.x]
        b = [point.y * c for c in self.y]
        c = [point.z * c for c in self.z]
        
        if(self.retMode == RET_MODE_VECTOR):
            return Vec({'xyz'[i]: self.o[i] + a[i] + b[i] + c[i] for i in range(3)})
        else:
            return [self.o[i] + a[i] + b[i] + c[i] for i in range(3)]
    
    def transformInto(self, point):
        """Algorithm taken from Dartmouth class CS77 framework code by Jonathan Denning, http://www.cs.dartmouth.edu/~cs77/assignments/assignment02.zip"""
        translated = [point.__dict__['xyz'[i]] - self.o[i] for i in range(3)]
        if(self.retMode == RET_MODE_VECTOR):
            return Vec({c: numpy.dot(translated, self.__dict__[c]) for c in 'xyz'})
        else:
            return [numpy.dot(translated, self.__dict__[c]) for c in 'xyz']
    
    def transformTo(self, other, point):
        return other.transformInto(Vec.from_array(self.transformOutOf(point)))
    
    @classmethod
    def createFromVectors(cls, origin, x, y):
        origin = [origin.__dict__[c] for c in 'xyz']
        x = [x.__dict__[c] for c in 'xyz']
        y = [y.__dict__[c] for c in 'xyz']
        norm = numpy.linalg.norm(x)
        x = [c / norm for c in x]
        norm = numpy.linalg.norm(y)
        if(math.isnan(norm) or norm == 0):
            print y
            print "ERROR! is nan or 0"
            raise Exception("Is NaN or 0")
        y = [c / norm for c in y]
        z = numpy.cross(x, y)
        norm = numpy.linalg.norm(z)
        if norm == 0:
            raise "Illegal arguments! x and y cannot be parallel" 
        z = [c / norm for c in z]
        y = numpy.cross(z, x)
        norm = numpy.linalg.norm(y)
        y = [c / norm for c in y]
        return TransformFrame(origin,x,y,z)
    
    def __str__(self):
        return ("frame[o=%s,x=%s,y=%s,z=%s]" % (strVecArray(self.o), strVecArray(self.x), strVecArray(self.y), strVecArray(self.z)))
    
identity_frame = TransformFrame.createFromVectors(Vec({
             "x": 0,
             "y": 0,
             "z": 0
             }), Vec({
             "x": 1,
             "y": 0,
             "z": 0
             }), Vec({
             "x": 0,
             "y": 1,
             "z": 0
             }))
    
def strVecArray(vecArray):
    return ("(%s,%s,%s)" % (str(vecArray[0]), str(vecArray[1]), str(vecArray[2])))
    
    
############################################
"""Debug code"""

if __name__ == "__main__":
    xVec = Vec({
             "x": 1.1,
             "y": 23,
             "z": 2
             })
    yVec = Vec({
             "x": 0.34,
             "y": 1,
             "z": 4
             })
    xVec2 = Vec({
             "x": 0,
             "y": 0,
             "z": 1
             })
    yVec2 = Vec({
             "x": 0,
             "y": 1,
             "z": 0
             })
    origin = Vec({
              "x": 0,
              "y": 0,
              "z": 0
              })
    
    point = Vec({
             "x": 1,
             "y": 2,
             "z" : 1
             })
    
    frame = TransformFrame.createFromVectors(origin, xVec, yVec)
    frame2 = TransformFrame.createFromVectors(origin, xVec2, yVec2)
    print "Frame: " + str(frame)
    print "Frame2: " + str(frame2)
    print "Point: " + str(point)
    print "Transform to: " + str(identity_frame.transformTo(frame, point))
    print "Transform into: " + str(frame.transformInto(point))
    print "Transform out of: " + str(frame.transformOutOf(point))
