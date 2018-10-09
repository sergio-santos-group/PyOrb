from distutils.core import setup, Extension
import numpy.distutils.misc_util


import os
# os.environ["CC"] = "/usr/local/gfortran/bin/gcc"
# python setup.py build_ext  --inplace

c_ext = Extension("_cube",
                  ["_cube.c", "cube.c"],
                  #extra_compile_args=['-O3','-fopenmp'],
                  #extra_compile_args=['-O3','-fopenmp'],
                  #extra_link_args=['-lgomp']
                  )

setup(
    ext_modules=[c_ext],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
