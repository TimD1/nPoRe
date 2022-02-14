from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [ 
        Extension("aln", ["src/aln.pyx"]),
        Extension("cig", ["src/cig.pyx"]),
        Extension("bam", ["src/bam.pyx"]) 
]

setup(
        name="nPoRe",
        ext_modules = cythonize(extensions, language_level='3str')
)
