from setuptools import setup, Extension
import setuptools
import os

from distutils.command.install import install


class CustomInstall(install):
    # Compile the C++ code
    def run(self):
        src_dir=os.path.dirname(os.path.abspath(__file__))
        src_file=src_dir + '/ChromA/util/FwdBwdRowMajor.cpp'
        so_file=src_dir + '/ChromA/util/libfwdbwdcpp.so'
        cmd='g++ {} -o {} --shared -fPIC -DNDEBUG -O3 -I{}'.format(src_file,so_file,src_dir + '/ChromA/util/eigen')
        print('=============================================================================')
        print(cmd)
        os.system(cmd)

        install.run(self)
# export EIGENPATH=/Users/mgabitto/Desktop/Projects/PeakCalling/PeakCalling/Chrom/util/eigen/
# g++ FwdBwdRowMajor.cpp -o libfwdbwd.so --shared -fPIC -DNDEBUG -O3 -I$EIGENPATH


setup(
    name='ChromA',
    version='0.0.2',
    #packages=['ChromA'],
    packages=setuptools.find_packages(),
    # note that we need to explicitly list the .so file so it gets copied
    package_data={'': ['data/*','data/promoters/*','data/blacklisted/*','util/*.so']}, 
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
    cmdclass={'install': CustomInstall}
    # ext_modules = [core_hmm]
)
