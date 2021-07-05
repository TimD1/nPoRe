from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [ 
        Extension("aln", ["aln.pyx"]),
        Extension("cig", ["cig.pyx"]) 
]

setup(
        name="realign",
        ext_modules = cythonize(extensions, language_level='3str')
)
