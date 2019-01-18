from setuptools import setup, Extension


cpp_args = ['-Wall']
core_hmm = Extension('HMM',
                     define_macros=[('MAJOR_VERSION', '1'),
                                    ('MINOR_VERSION', '0')],
                     include_dirs=['/usr/local/include'],
                     libraries=['tcl83'],
                     library_dirs=['/usr/local/lib'],
                     sources=['/util/FwdBwdRowMajor.cpp'],
                     extra_compile_args=cpp_args
                     )
# export EIGENPATH=/Users/mgabitto/Desktop/Projects/PeakCalling/PeakCalling/Chrom/util/eigen/
# g++ FwdBwdRowMajor.cpp -o libfwdbwd.so --shared -fPIC -DNDEBUG -O3 -I$EIGENPATH


setup(
    name='ChromA',
    version='0.0.2',
    packages=['', 'data', 'util', 'classes'],
    package_dir={'': 'Chrom_light2'},
    package_data={'': ['data/*']},
    url='',
    license='',
    author='Mariano Gabitto',
    author_email='mgabitto@simonsfoundation.org',
    description='Chromatin Annotation Tool',
    install_requires=['matplotlib', 'seaborn', 'numpy', 'pysam', 'ray', 'scipy']
    # ext_modules = [core_hmm]
)
