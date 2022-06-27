from setuptools import find_packages, setup
from Cython.Build import cythonize

setup(
    name='gaden',
    ext_modules = cythonize(["src/Simulation.pyx", "src/Utils.pyx"]),
    version='0.1.0',
    description='Python API to directly access GADEN simulation results, without depending on ROS.',
    author='Pepe Ojeda',
    license='MIT',
)

