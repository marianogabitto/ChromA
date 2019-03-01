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


def _post_install(_):
    src_dir = os.path.dirname(os.path.abspath(__file__))
    src_file = src_dir + '/ChromA/util/FwdBwdRowMajor.cpp'
    so_file = src_dir + '/ChromA/util/libfwdbwdcpp.so'
    cmd = 'g++ {} -o {} --shared -fPIC -DNDEBUG -O3 -I{}'.format(src_file, so_file, src_dir + '/ChromA/util/eigen')
    print('=============================================================================')
    print(cmd)
    print('=============================================================================')
    os.system(cmd)
    os.system("cp {} {}".format(so_file, src_dir + '/ChromA/data/libfwdbwdcpp.so'))
    print(os.system("ls {}".format(src_dir + '/ChromA/util/')))

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        self.execute(_post_install, " ", msg="Compiling C++ Core")

class CustomDevelopCommand(develop):
    def run(self):
        develop.run(self)
        self.execute(_post_install, " ", msg="Compiling C++ Core")

class CustomEggInfoCommand(egg_info):
    def run(self):
        egg_info.run(self)
        self.execute(_post_install, " ", msg="Compiling C++ Core")

class CustomInstall(install):
    # Compile the C++ code
    def run(self):
        install.run(self)
        src_dir=os.path.dirname(os.path.abspath(__file__))
        src_file=src_dir + '/ChromA/util/FwdBwdRowMajor.cpp'
        so_file=src_dir + '/ChromA/util/libfwdbwdcpp.so'
        cmd='g++ {} -o {} --shared -fPIC -DNDEBUG -O3 -I{}'.format(src_file,so_file,src_dir + '/ChromA/util/eigen')
        print('=============================================================================')
        print(cmd)
        os.system(cmd)
        os.system("pwd")
        print(os.system("ls"))
        # install.run(self)
        """
        self.execute(_post_install, " ", msg="Compiling C++ Core")
        """

setup(
    name='ChromA',
    version='0.0.5',
    packages=setuptools.find_packages(),
    # note that we need to explicitly list the .so file so it gets copied
    package_data={'': ['test/*', 'data/*', 'data/promoters/*', 'data/blacklisted/*',
                       'util/*.so']},
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
    test_suite='nose.collector',
    cmdclass={
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
        'egg_info': CustomEggInfoCommand
    }
)
