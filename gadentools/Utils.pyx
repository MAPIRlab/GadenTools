from typing import Tuple
import cython

cdef class Vector3:
    def __cinit__(self, x:float, y:float, z:float):
        self.x = x
        self.y = y
        self.z = z
    
    @classmethod
    @cython.returns(Vector3)
    def fromTuple(cls, tup: Tuple) :
        return Vector3.__new__(Vector3, tup[0], tup[1], tup[2])

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
    
