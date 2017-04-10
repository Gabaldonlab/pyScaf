"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

from setuptools import find_packages

NAME="pyScaf"

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name=NAME,
      version='0.12a4',
      description='',
      long_description=long_description,
      author='Leszek Pryszcz',
      author_email='l.p.pryszcz+distutils@gmail.com',
      url='https://github.com/lpryszcz/%s'%NAME,
      license='GPLv3',

      install_requires=["FastaIndex", ],
      #tests_require=["FastaIndex", ],
      #packages=find_packages(),
      
      py_modules=[NAME, ],

      entry_points = {"console_scripts":
                      ["%s = %s:main" % (NAME, NAME), 
                       ]
                     }, 
      
      classifiers=[
          # How mature is this project? Common values are  3 - Alpha 4 - Beta 5 - Production/Stable
          'Development Status :: 4 - Beta',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
      ],

      keywords='assembly scaffolding paired-end long-reads synteny',
      
     )
