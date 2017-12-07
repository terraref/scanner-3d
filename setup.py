#!/usr/bin/env python

from setuptools import setup, find_packages

from codecs import open
from os import path

here=path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description=f.read()

print find_packages()

setup(
        name='scanner_3d',
        version='1.0.0',

        description='useful 3d-scanner modules that contain extractor functions',
        long_description=long_description,

        entry_points={
            'console_scripts': [
                'terra_ply2las.py=ply2las.terra_ply2las:main',
            ],
        },

        install_requires=[
            'pika>=0.10.0',
            'requests>=2.11.0',
        ],

        dependency_links=['https://opensource.ncsa.illinois.edu/bitbucket/rest/archive/latest/projects/CATS/repos/pyclowder2/archive?format=zip'],

        packages=find_packages(),

        # basic package metadata
        url='https://github.com/terraref/scanner_3d',
        author='Max Burnette',
        author_email='mburnet2@illinois.edu',

        license='NCSA',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: NCSA License',

            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
        ],
        keywords='terraref 3d scanner extractor'

)
