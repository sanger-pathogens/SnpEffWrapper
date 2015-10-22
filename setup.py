import os
from setuptools import setup, find_packages
import multiprocessing

setup(name='annotateVCF',
      version='0.0.0',
      scripts=[
      ],
      test_suite='nose.collector',
      tests_require=[
        'nose',
        'mock'
      ],
      install_requires=[
        'PyYAML',
      ],
      include_package_data=True,
      package_data={
        'data': 'annotateVCF/data/*',
        'test_data': 'annotateVCF/tests/data/*'
      },
      packages=find_packages(),
      zip_safe=False
)
