import numpy
from setuptools import setup
from Cython.Build import cythonize, build_ext

setup(
    name='cascy_ded',
    version='0.0.0',
    author='Benjamin Spreng',
    author_email='sprengjamin@gmail.com',
    packages=['cascy_ded'],
    license='LICENSE.txt',
    description='Casimir interaction involving cylinders',
    long_description=open('README.md').read(),
    install_requires=[],
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(
        ['cascy_ded/scattering_amplitude.pyx', 'cascy_ded/bessel.pyx'],
        annotate=True),
    include_dirs=[numpy.get_include()],
    options = {'build_ext':{'inplace':True, 'force':True}},
)