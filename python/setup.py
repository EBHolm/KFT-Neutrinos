from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    name='kftneutrinos',
    version='23.0.0',
    description='Python interface to KFT Neutrinos',
    url='',
    ext_modules = cythonize(
        Extension('kftneutrinos',
            sources=['../source/perturbation_first.cpp', '../source/perturbation_second_kft.cpp',
                     '../source/quadrature.cpp', '../source/cosmology.cpp', 'kftneutrinos.pyx'],
            include_dirs=['../source/'],
            language='c++',
            extra_compile_args=['-std=c++20', '-O3', '-ffast-math'],
            extra_link_args=['-std=c++20']
        ),
        language_level=3,
        annotate=True
    )
)