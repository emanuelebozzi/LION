# setup.py

from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("waveform_stacking.pyx")
)
