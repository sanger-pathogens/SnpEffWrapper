import os
from setuptools import setup, find_packages
import multiprocessing

setup(name='snpEffWrapper',
      version='0.2.2',
      scripts=[
        'scripts/snpEffBuildAndRun'
      ],
      install_requires=[
        'Jinja2',
        'PyVCF',
        'PyYAML'
      ],
      include_package_data=True,
      package_data={
        'data': 'snpEffWrapper/data/*',
        'test_data': 'snpEffWrapper/tests/data/*'
      },
      packages=find_packages(),
      zip_safe=False
)
