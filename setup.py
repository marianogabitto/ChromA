from setuptools import setup, Extension
import setuptools

extensions = [Extension("libfwdbwdcpp",
                        ["ChromA/util/FwdBwd.cpp"],
                        include_dirs=["ChromA/util/eigen"])]

setup(
      name='ChromA',
      version='3.0',
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
                        'numpy>=1.16.2',
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
