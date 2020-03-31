from setuptools import setup

import os

base_path = os.path.dirname(__file__)

setup(
    name='usdaars-primerlib',
    version='1.0',
    description='A primer library for calculating primers for Nested Loop PCR and STARP.',
    author='Kaleb Burnham',
    author_email='kaleb.burnham@usda.gov',
    packages=['nestedloop', 'starp'],
    package_dir={"": "src"},
    install_requires=[
        'Biopython',
        'bs4',
        'regex',
        'lxml',
   ]
)