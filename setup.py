#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
from projman_filler import __version__ as projman_filler_version
import subprocess
import re
import logging

# 1. Configure the logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

with open('README.md') as readme_file:
    readme = readme_file.read()


def get_latest_version(package_name):
    try:
        # Run the pip index versions command
        result = subprocess.run(
            ["python", "-m", "pip", "index", "versions", "--pre", package_name],
            capture_output=True, text=True, check=True
        )
        
        # The first line usually contains "package_name (X.Y.Z)"
        # Or parse the "Available versions:" line
        output = result.stdout
        
        # Regex to find the first version block, e.g., "package_name (1.2.3)"
        match = re.search(rf"{package_name}\s*\(([^)]+)\)", output)
        if match:
            return match.group(1)
            
        # Fallback: parse the list of available versions
        versions_match = re.search(r"Available versions:\s*([^\n]+)", output)
        if versions_match:
            # Grab the first version in the comma-separated list
            versions = [v.strip() for v in versions_match.group(1).split(",")]
            return versions[0]
            
    except (subprocess.CalledProcessError, FileNotFoundError, IndexError) as e:
        # Fallback version if the pip command fails or is unavailable
        logger.exception(f"Failed to fetch latest version for package '{package_name}' via pip.")
        return '4.1.1rc1' 
    
latest_checkqc_version = get_latest_version("checkQC")
logger.info(f"Using latest checkQC version = {latest_checkqc_version}")

requirements = [
    'Click~=8.1.8',
    'interop~=1.4.0',
    'SQLAlchemy~=2.0.40',
    'pymssql~=2.3.4',
    'pandas~=2.2.3',
    'numpy~=2.2.4',
    'checkQC@git+https://github.com/nkongenelly/checkQC.git@DATAOPS-1796_fix_bclconvert_sequencing_metrics#egg=checkQC',

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
