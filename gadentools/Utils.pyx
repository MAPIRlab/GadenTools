from typing import Tuple
import cython
import numpy
cimport numpy

cpdef numpy.ndarray toColorImage(numpy.ndarray[object, ndim=2] input) :
    """
    Takes a np.ndarray of Vector3 and turns them into 3-elements ndarrays
    Meant for generating color images from wind maps
    """
    cdef numpy.ndarray[numpy.float64_t, ndim=3] newarray = numpy.empty((input.shape[0], input.shape[1], 3), float)
    for i in range(input.shape[0]):
        for j in range(input.shape[1]):
            newarray[i,j,:] = Vector3.toNdArray(input[i,j])
    return newarray


def block(original, image):
    """
    Takes an image and draws the obstacles on it in black. "original" must be a concentration or wind map that was generated with the "blockObstacles" option.
    """
    for i in range(original.shape[0]):
        for j in range(original.shape[1]):
            if __isBlocked(original[i,j]) :
                if len(image.shape)==3:
                    image[i,j,:] = 0
                else:
                    image[i,j] = 0

def __isBlocked(value):
    '''Returns True if value corresponds to a cell that was marked as blocked by generateConcentrationMap or generateWindMap.'''
    if isinstance(value, float) or isinstance(value, int):
        return value<0
    else:
        return False

cpdef numpy.ndarray makeSobelKernelX(int size):
    cdef int halfsize = size//2
    cdef numpy.ndarray[numpy.float_t, ndim=2] newarray = numpy.full((size, size), 0, float)
    for i in range(-halfsize, halfsize+1):
        for j in range(-halfsize, halfsize+1):
            newarray[halfsize+i, halfsize+j] = 0 if j==0 else j / (i*i + j*j)
    return newarray

cpdef numpy.ndarray makeSobelKernelY(int size):
    cdef int halfsize = size//2
    cdef numpy.ndarray[numpy.float_t, ndim=2] newarray = numpy.full((size, size), 0, float)
    for i in range(-halfsize, halfsize+1):
        for j in range(-halfsize, halfsize+1):
            newarray[halfsize+i, halfsize+j] = 0 if i==0 else i / (i*i + j*j)
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

    cpdef float magnitude(Vector3 self):
        return sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    cpdef Vector3 normalized(Vector3 self) :
        return self / sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
    
    cpdef Vector3 copy(Vector3 self)  :
        return Vector3.__new__(Vector3, self.x, self.y, self.z)
    
    cpdef Vector3 round(Vector3 self, float p) :
        return Vector3.__new__(Vector3, round(self.x,p), round(self.y,p), round(self.z,p) )
