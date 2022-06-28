from typing import Tuple
import cython
import numpy
cimport numpy

cpdef numpy.ndarray toColorImage(numpy.ndarray[object, ndim=2] input) :
    """
    Takes a np.ndarray of Vector3 and turns them into 3-elements ndarrays
    Meant for generating color images from wind maps
    """
    cdef numpy.ndarray[object, ndim=3] newarray = numpy.empty((input.shape[0], input.shape[1], 3), object)
    for i in range(input.shape[0]):
        for j in range(input.shape[1]):
            newarray[i,j,:] = Vector3.toNdArray(input[i,j])
    return newarray


cdef class Vector3:
    def __cinit__(self, x:float, y:float, z:float):
        self.x = x
        self.y = y
        self.z = z
    
    @classmethod
    @cython.returns(Vector3)
    def fromTuple(cls, tup: Tuple) :
        return Vector3.__new__(Vector3, tup[0], tup[1], tup[2])

    @classmethod
    @cython.returns(Tuple)
    def toTuple(cls, vec: Vector3) :
        return (vec.x, vec.y, vec.z)

    @classmethod
    @cython.returns(Tuple)
    def toNdArray(cls, vec) :
        if isinstance(vec, Vector3):
            return numpy.array([vec.x, vec.y, vec.z])
        else:
            return numpy.array([0,0,0])

    @cython.returns(Vector3)
    def __add__(Vector3 self, Vector3 other) :
        return Vector3.__new__(Vector3, self.x + other.x, self.y + other.y, self.z + other.z)
    
    @cython.returns(Vector3)
    def __neg__(Vector3 self) :
        return Vector3.__new__(Vector3, -self.x, -self.y, -self.z)
    
    @cython.returns(Vector3)
    def __sub__(Vector3 self, Vector3 other) :
        return Vector3.__new__(Vector3, self.x - other.x, self.y - other.y, self.z - other.z)

    @cython.returns(Vector3)
    def __mul__(Vector3 self, float scalar) :
        return Vector3.__new__(Vector3, self.x * scalar, self.y * scalar, self.z * scalar)
    
    @cython.returns(Vector3)
    def __truediv__(Vector3 self, float scalar) :
        return Vector3.__new__(Vector3, self.x /scalar, self.y /scalar, self.z /scalar)

    @cython.returns(Vector3)
    def __floordiv__(Vector3 self, float scalar) :
        return Vector3.__new__(Vector3, self.x // scalar, self.y // scalar, self.z // scalar)
    
    def __repr__(self):
        return "("+str(self.x)+", "+str(self.y)+", "+str(self.z)+")"
    
    cpdef float dot(self, Vector3 other):
        return self.x*other.x+self.y*other.y+self.z*other.z

    cpdef Vector3 projectOnVector(self, Vector3 other):
        return self.dot(other) * other.normalized()

    cpdef Vector3 projectOnPlane(self, Vector3 planeNormal):
        return self-self.projectOnVector(planeNormal)
