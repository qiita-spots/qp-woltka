#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from setuptools import setup


__version__ = "2024.09"


classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""


with open('README.rst') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='qp-woltka',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Qiita Plugin: Woltka',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/biocore/qiita',
      test_suite='nose.collector',
      packages=['qp_woltka'],
      package_data={
        'qp_woltka': [
            'support_files/*',
            'databases/*']},
      scripts=['scripts/configure_woltka', 'scripts/start_woltka',
               'scripts/finish_woltka', 'scripts/woltka_merge'],
      extras_require={'test': ["nose >= 0.10.1", "pep8"]},
      install_requires=['click >= 3.3', 'future', 'pandas >= 0.15', 'duckdb',
                        'h5py >= 2.3.1', 'biom-format', 'lxml',
                        'polars-lts-cpu', 'html5lib',
                        'woltka @ git+https://github.com/qiyunzhu/'
                        'woltka.git#egg=woltka',
                        'pysyndna @ git+https://github.com/AmandaBirmingham/'
                        'pysyndna.git#egg=pysyndna',
                        'mxdx @ git+https://github.com/wasade/'
                        'mxdx.git#egg=mxdx',
                        'micov @ git+https://github.com/biocore/'
                        'micov.git#egg=micov'],
      classifiers=classifiers
      )
