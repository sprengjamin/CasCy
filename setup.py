import numpy
from setuptools import setup
from Cython.Build import cythonize, build_ext

setup(
    name='CasCy',
    version='0.0.0',
    author='Benjamin Spreng',
    author_email='sprengjamin@gmail.com',
    packages=['CasCy'],
    license='LICENSE.txt',
    description='Casimir interaction involving cylinders',
    long_description=open('README.md').read(),
    install_requires=[],
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(
        ['CasCy/scattering_amplitude.pyx', 'CasCy/bessel.pyx'],
        annotate=True),
    include_dirs=[numpy.get_include()],
    compiler_directives={'embedsignature': True},
    options = {'build_ext':{'inplace':True, 'force':True}},
)