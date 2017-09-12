#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    'Click>=6.0',
    'xmltodict',
    'interop'
]

setup_requirements = [
    'pytest-runner',
]

test_requirements = [
    'pytest',
]

setup(
    name='projman_filler',
    version='0.1.0',
    description="App for filling ProjMan",
    long_description=readme,
    author="Johan Dahlberg",
    author_email='johan.dahlberg@medsci.uu.se',
    url='https://github.com/johandahlberg/projman_filler',
    packages=find_packages(include=['projman_filler']),
    entry_points={
        'console_scripts': [
            'projman_filler=projman_filler.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="GNU General Public License v3",
    zip_safe=False,
    keywords='projman_filler',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
