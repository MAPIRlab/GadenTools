from libc.math cimport sqrt, floor

cdef class Vector3:
    cdef public float x, y ,z

        
    cdef inline float magnitude(Vector3 self):
        return sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    cdef inline Vector3 normalized(Vector3 self) :
        return self / sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
    
    cdef inline Vector3 copy(Vector3 self)  :
        return Vector3.__new__(Vector3, self.x, self.y, self.z)
    