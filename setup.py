#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['toil[all]==5.4.0', ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Ronak Shah",
    author_email='rons.shah@gmail.com',
    classifiers=[
        'Development Status :: 2 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
    description="Workflow that creates all the bam files for the MSK-ACCESS fastq file",
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + "Workflow that creates all the bam files for the MSK-ACCESS fastq file",
    include_package_data=True,
    keywords='fastq_to_bam',
    name='fastq_to_bam',
    packages=find_packages(include=['nucleo']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/msk-access/nucleo',
    version='3.0.2',
    zip_safe=False,
)
