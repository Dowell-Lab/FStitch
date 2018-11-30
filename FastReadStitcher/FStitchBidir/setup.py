#!/usr/bin/env python

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
#with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()

setup(name='FStitchBidir',
      version='0.1',
      description='Annotates bidirections using FStitch segment data.',
      url='http://github.com/storborg/funniest',
      author='Margaret Gruca',
      author_email='margaret.gruca@colorado.edu',
      license='MIT',
      python_requires='>=3.6',
      classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3'
      ],
      
keywords='bioinformatics genomics chromatin ATAC-seq motif transcription_factor',

    package_dir={'bidir' : 'bidir'},
    packages=['FStitchBidir'],

    install_requires=[
        'argparse',
        'datetime',
        'numpy',
        'pybedtools',
        'chartify',
        'pandas',
    ],

    entry_points={
        'console_scripts': [
            'bidir=FStitch:bidir',
        ],
    },

    project_urls={
        'Bug Reports': 'https://github.com/Dowell-Lab/FStitch',
        'Source': 'https://github.com/Dowell-Lab/FStitch',
        'Human RefSeq Gene Annotation Download': 'http://genome.ucsc.edu/cgi-bin/hgTables',
    },
)