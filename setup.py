from setuptools import setup, Extension
import setuptools
import os

src_dir=os.path.dirname(os.path.abspath(__file__))

#cpp_args = ['--shared','-fPIC']
cpp_args = ['-O3']
core_hmm = Extension('ChromA_HMM',
                     define_macros=[('MAJOR_VERSION', '1'),
                                    ('MINOR_VERSION', '0')],
                     include_dirs=['/usr/local/include',src_dir+'/ChromA/util/eigen'],
                     libraries=[],
                     library_dirs=['/usr/local/lib'],
                     sources=[src_dir+'/ChromA/util/FwdBwdRowMajor.cpp'],
                     extra_compile_args=cpp_args
                     )
# export EIGENPATH=/Users/mgabitto/Desktop/Projects/PeakCalling/PeakCalling/Chrom/util/eigen/
# g++ FwdBwdRowMajor.cpp -o libfwdbwd.so --shared -fPIC -DNDEBUG -O3 -I$EIGENPATH


setup(
    name='ChromA',
    version='0.0.2',
    #packages=['ChromA'],
    packages=setuptools.find_packages(),
    package_data={'': ['data/*','data/promoters/*','data/blacklisted/*']},
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
        'psutil'
    ],
    ext_modules = [core_hmm]
)
