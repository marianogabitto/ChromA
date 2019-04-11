from setuptools import setup, Extension
import setuptools
import os

from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info

extensions = [Extension("libfwdbwdcpp",
                        ["ChromA/util/FwdBwdRowMajor.cpp"],
                        include_dirs=["ChromA/util/eigen"])]

setup(
      name='ChromA',
      version='0.0.9',
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
