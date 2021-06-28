from setuptools import setup, Extension
from Cython.Build import cythonize

# https://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line
from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()
directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

extensions = [
        Extension("aln", ["aln.pyx"],
            define_macros=[('CYTHON_TRACE', '1')]
        )
]

setup(
        name="aln",
        ext_modules = cythonize(extensions, language_level='3str')
)
