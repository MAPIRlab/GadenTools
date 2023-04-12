from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

extensions = Extension("*", ["gadentools/*.pyx"],
    include_dirs=[numpy.get_include()]
)

setup(
    name='gadentools',
    ext_modules = cythonize([extensions]),
    version='0.1.0',
    description='Python API to directly access GADEN simulation results, without depending on ROS.',
    author='Pepe Ojeda',
    license='MIT',
)

