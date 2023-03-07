#!/usr/bin/env python3

import setuptools
from kmerator import info

setuptools.setup(
    name = 'kmerator',
    version = info.VERSION,
    author = info.AUTHOR,
    author_email = info.AUTHOR_EMAIL,
    description = info.SHORTDESC,
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    url="https://github.com/Transipedia/kmerator",
    packages = setuptools.find_packages(),
    classifiers = [
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
    ],
    entry_points = {
        'console_scripts': [
            'kmerator = kmerator.kmerator:main',
        ],
    },
    include_package_data = True,
    install_requires=['bs4', 'lxml', 'requests'],
    python_requires = ">=3.7",
    licence = "GPLv3"
)
