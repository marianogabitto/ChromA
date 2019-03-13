from setuptools import setup, Extension
import setuptools
import os

#from distutils.command.install import install
#from distutils.command.build_ext import build_ext

from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info

#############################################
# HOW TO COMPILE THE CPP CORE
# export EIGENPATH=/Users/mgabitto/Desktop/Projects/PeakCalling/PeakCalling/Chrom/util/eigen/
# g++ FwdBwdRowMajor.cpp -o libfwdbwd.so --shared -fPIC -DNDEBUG -O3 -I$EIGENPATH
#############################################

extensions = [Extension("libfwdbwdcpp",
			["ChromA/util/FwdBwdRowMajor.cpp"],
			include_dirs=["ChromA/util/eigen"])]

setup(
    name='ChromA',
    version='0.0.3',
    packages=setuptools.find_packages(),
    # note that we need to explicitly list the .so file so it gets copied
    package_data={'': ['test/*', 'data/*', 'data/promoters/*', 'data/blacklisted/*',
                       'util/libfwdbwdcpp.so', 'util/*.so']},
    url='',
    license='',
    author='Mariano Gabitto',
    author_email='mgabitto@simonsfoundation.org',
    description='Chromatin Annotation Tool',
    scripts=['bin/ChromA'],
    install_requires=[
        'matplotlib',
        'seaborn',
        'numpy',
        'pysam',
        'ray',
        'scipy',
        'setproctitle',
        'psutil',
        'nose'
    ],
    ext_modules=extensions,
    test_suite='nose.collector',
)
