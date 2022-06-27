#cython: boundscheck=False
#cython: wraparound=False
import struct
from typing import Tuple
import zlib
import numpy
cimport numpy
import math
from Utils cimport Vector3
import time
cimport cython

cdef inline int indexFrom3D(float x, float y, float z, int number_of_cells_x, int number_of_cells_y) :
        return (<int>z) * number_of_cells_x * number_of_cells_y + (<int>y) * number_of_cells_x + (<int>x)

cdef class Filament:
    cdef public Vector3 position
    cdef public float stdDev 
    def __cinit__(self, value : Tuple):
        self.position = Vector3(value[0],value[1],value[2])
        self.stdDev = value[3]


cdef class Simulation:
    cdef object simulationFolder
    cdef object occupancyFile
    cdef int currentIteration
    cdef int currentWind
    cdef dict filaments
    cdef numpy.ndarray U
    cdef numpy.ndarray V
    cdef numpy.ndarray W
    cdef numpy.ndarray Environment
    cdef Vector3 env_min
    cdef Vector3 env_max
    cdef int number_of_cells_x
    cdef int number_of_cells_y
    cdef int number_of_cells_z
    cdef float cell_size
    cdef int gas_type
    cdef float total_moles_in_filament
    cdef float num_moles_all_gases_in_cm3


    def __init__(self, simulationFolder : str, occupancyFile : str):
        """
        simulationFolder: Path to the folder that contains the iteration_i files
        
        occupancyFile: Path to the OccupancyGrid3D.csv
        """
        self.simulationFolder = simulationFolder
        self.occupancyFile = occupancyFile
        self.currentIteration=-1
        self.currentWind=-1
        self.__readFile(0)
        self.__loadOccupancyFile()

    cpdef float getCurrentConcentration(self, Vector3 location):
        '''Returns concentration (in ppm) at location in the most recently loaded iteration.'''
        cdef int[:] env = self.Environment #memoryview
        return self.__getConcentration(self.currentIteration, location, env)

    cpdef float getConcentration(self, int iteration, Vector3 location):
        '''Returns concentration (in ppm) at location in the most recently loaded iteration.'''
        cdef int[:] env = self.Environment #memoryview
        return self.__getConcentration(iteration, location, env)

    cdef float __getConcentration(self, int iteration, Vector3 location, int[:] env):
        '''Returns concentration (in ppm) at location in the specified iteration. If the iteration is not loaded, it gets loaded'''
        self.__readFile(iteration) #if it's the current iteration it doesn't do anything
        cdef float total = 0.0
        cdef Filament fil
        for ind, filament in self.filaments.items() :
            fil = filament
            if (<Vector3>(fil.position-location)).magnitude() > fil.stdDev/100 * 5:
                continue
            if self.__checkPath(location, fil.position, env) :
                total += self.__getConcentrationFromFilament(location, fil)
        return total

    cpdef numpy.ndarray generateConcentrationMap2D(self, int iteration, float height):
        '''Returns 2D concentration map (in ppm) as a numpy.array of floats. Kind of slow.'''
        cdef numpy.ndarray[numpy.float_t, ndim=2] concentration_map = numpy.zeros((self.number_of_cells_x, self.number_of_cells_y), float)
        cdef int k = int( (height-self.env_min.z) / self.cell_size )

        cdef int[:] env = self.Environment #memoryview
        cdef Vector3 location
        for i in range(concentration_map.shape[0]) :
            for j in range(concentration_map.shape[1]) : 
                    if self.Environment[indexFrom3D(i,j,k, self.number_of_cells_x,self.number_of_cells_y)]:
                        location = Vector3(i + 0.5, j + 0.5, k + 0.5) * self.cell_size + self.env_min 
                        concentration_map[i,j] = self.__getConcentration(iteration, location, env)
        return concentration_map

    cpdef numpy.ndarray generateConcentrationMap3D(self, int iteration):
        '''Returns 3D concentration map (in ppm) as a numpy.array of floats. Quite slow.'''
        cdef numpy.ndarray[numpy.float_t, ndim=3] concentration_map = numpy.zeros((self.number_of_cells_x, self.number_of_cells_y, self.number_of_cells_z), float)
        for k in range(concentration_map.shape[2]) :
            concentration_map[:,:,k] = self.generateConcentrationMap2D(iteration, k*self.cell_size + self.env_min.z)
        return concentration_map



    cdef float __getConcentrationFromFilament(self, Vector3 location, Filament filament) :
        cdef float distance_cm = 100 * (<Vector3>(location-filament.position)).magnitude();

        cdef float num_moles_target_cm3 = (self.total_moles_in_filament /
            (math.sqrt(8*(math.pi**3)) * (filament.stdDev**3) )) * math.exp( -(distance_cm**2)/(2*(filament.stdDev**2)) );

        cdef float ppm = num_moles_target_cm3/self.num_moles_all_gases_in_cm3 * 1000000; #parts of target gas per million

        return ppm;

    cdef __readFile(self, int iteration):
        '''Safe to call without checking if iteration is the current one.'''
        if iteration == self.currentIteration :
            return

        self.currentIteration = iteration
        self.filaments = {}
        
        cdef object data
        with open(self.simulationFolder+"/iteration_"+str(iteration), "rb") as file:
            data = zlib.decompress(file.read())
        
        cdef int index = struct.calcsize("=i")
        self.env_min = Vector3.fromTuple( struct.unpack_from("=ddd", data, index) )
        index += struct.calcsize("=ddd")
        self.env_max = Vector3.fromTuple( struct.unpack_from("=ddd", data, index) )
        index += struct.calcsize("=ddd")

        cdef tuple num_cells = struct.unpack_from("=iii", data, index)
        self.number_of_cells_x = num_cells[0]
        self.number_of_cells_y = num_cells[1]
        self.number_of_cells_z = num_cells[2]
        index += struct.calcsize("=iii")
        
        self.cell_size = struct.unpack_from("=d", data, index)[0]
        index += struct.calcsize("=d")
        
        index += 5 * struct.calcsize("=d") #discard these bytes

        self.gas_type = struct.unpack_from("=i", data, index)[0]
        index += struct.calcsize("=i")

        self.total_moles_in_filament = struct.unpack_from("=d", data, index)[0]
        index += struct.calcsize("=d")
        self.num_moles_all_gases_in_cm3 = struct.unpack_from("=d", data, index)[0]
        index += struct.calcsize("=d")

        cdef int wind_iteration = struct.unpack_from("=i", data, index)[0]
        index += struct.calcsize("=i")
        
        cdef filament_ind
        while data.__sizeof__()>index+struct.calcsize("=idddd"):
            filament_ind = struct.unpack_from("=i", data, index)[0]
            index += struct.calcsize("=i")
            self.filaments[int(filament_ind)] = Filament(struct.unpack_from("=dddd", data, index))
            index += struct.calcsize("=dddd")
            
        self.__loadWind(wind_iteration)

    cdef __loadWind(self, int wind_iteration):
        cdef object data
        cdef int index = 0
        cdef object format
        if wind_iteration is not self.currentWind:
            self.currentWind = wind_iteration
            with open(self.simulationFolder+"/wind/wind_iteration_"+str(wind_iteration), "rb") as file:
                data = file.read() 

            format = "="+ self.number_of_cells_x * self.number_of_cells_y *self.number_of_cells_z* "d"
            self.U = numpy.array( struct.unpack_from(format, data, index) )
            index += struct.calcsize(format)

            self.V = numpy.array( struct.unpack_from(format, data, index) )
            index += struct.calcsize(format)

            self.W = numpy.array( struct.unpack_from(format, data, index) )

    cdef __loadOccupancyFile(self):
        cdef object lines
        with open(self.occupancyFile, "r") as file:
            lines = file.read().splitlines()
        
        self.Environment = numpy.zeros(self.number_of_cells_x * self.number_of_cells_y * self.number_of_cells_z, numpy.intc)

        cdef int idx_x = 0
        cdef int idx_y = 0
        cdef int idx_z = 0

        cdef object line
        for lineIndex in range(4, len(lines)):
            line = lines[lineIndex].replace(" ","")
            if line == ";":
                idx_x = 0
                idx_y = 0
                idx_z+=1
                continue

            for ch in line:
                self.Environment[indexFrom3D(idx_x, idx_y, idx_z, self.number_of_cells_x,self.number_of_cells_y)] =int( not(bool(int(ch))) )
                idx_y+=1

            idx_x += 1
            idx_y = 0
                
    cdef bint __checkPath(self, Vector3 start, Vector3 end, int[:] env) :
        cdef Vector3 startIndex = (start-self.env_min) / self.cell_size
        cdef Vector3 endIndex = (end-self.env_min) / self.cell_size
        
        cdef Vector3 delta = end - start
        cdef int steps = <int>(delta.magnitude()/self.cell_size)
        if steps == 0 : 
            return True
        delta = delta / steps
        cdef Vector3 current = start.copy()
        cdef Vector3 currentIndex
        for i in range(steps):
            current = current + delta
            currentIndex = (current-self.env_min) / self.cell_size
            if not(env[indexFrom3D(currentIndex.x, currentIndex.y, currentIndex.z, self.number_of_cells_x,self.number_of_cells_y)]) :
                return False

        return True
    


