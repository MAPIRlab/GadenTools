from libc.math cimport sqrt, floor

cdef class Vector3:
    cdef public float x, y ,z
    
    cpdef float dot(self, Vector3 other)
    cpdef Vector3 cross(self, Vector3 other)
    cpdef Vector3 projectOnVector(self, Vector3 other)
    cpdef Vector3 projectOnPlane(self, Vector3 planeNormal)

        
    cpdef float magnitude(Vector3 self)
    cpdef Vector3 normalized(Vector3 self)     
    cpdef Vector3 copy(Vector3 self)      
    cpdef Vector3 roundTo(Vector3 self, float p)