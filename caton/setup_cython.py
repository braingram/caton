from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    cmdclass = {'build_ext': build_ext},
	include_dirs = [numpy.get_include(),'.'],
    ext_modules = [Extension("CEM_extensions", ["CEM_extensions.pyx"])]
)
