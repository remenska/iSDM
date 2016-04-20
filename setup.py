from __future__ import print_function

import six

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = open('README.md').read()

requirements = open('requirements.txt').read().splitlines()

if six.PY2:
    requirements.append('configparser')

setup(
    name='iSDM',
    version='0.0.1',
    url='http://github.com/remenska/iSDM/',
    description='interacting Species Distribution Modeling framework.',
    author='Daniela Remenska',
    author_email='d.remenska@esciencecenter.nl',
    packages = ['iSDM'],
    license='Apache Software License',
    tests_require=['nosetests'],
    platforms='any',
    long_description=readme,
    install_requires=requirements,
    classifiers = [
    'Programming Language :: Python',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    ],
    test_suite='tests'

    )