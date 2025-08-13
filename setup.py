#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

from projman_filler import __version__ as projman_filler_version

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'Click~=8.1.8',
    'interop~=1.4.0',
    'SQLAlchemy~=2.0.40',
    'pymssql~=2.3.4',
    'pandas~=2.2.3',
    'numpy~=2.2.4',
    # checkQC must be installed manually
    # git clone https://github.com/Molmed/checkQC.git
    # cd checkQC
    # sed -i 's/python_requires=.*$/python_requires=">=3.10",/' setup.py
    # sed -i 's/^interop~=.*$/interop==1.4.0/' requirements/prod
    # pip install --no-deps -r requirements/dev .

]

setup_requirements = [
    'pytest-runner',
]

test_requirements = [
    'pytest',
]

setup(
    name='projman_filler',
    version=projman_filler_version,
    description="App for filling ProjMan",
    long_description=readme,
    author="Johan Dahlberg",
    author_email='johan.dahlberg@medsci.uu.se',
    url='https://github.com/johandahlberg/projman_filler',
    packages=find_packages(exclude=["tests/"]),
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
        'Programming Language :: Python :: 3.13.1',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
