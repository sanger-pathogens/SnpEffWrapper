import os
from setuptools import setup, find_packages
import multiprocessing

setup(name='annotateVCF',
      version='0.1.0',
      scripts=[
        'scripts/annotateVCF'
      ],
      install_requires=[
        'Jinja2',
        'PyVCF',
        'PyYAML'
      ],
      include_package_data=True,
      package_data={
        'data': 'annotateVCF/data/*',
        'test_data': 'annotateVCF/tests/data/*'
      },
      packages=find_packages(),
      zip_safe=False
)
