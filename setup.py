import os
import re

from setuptools import find_packages, setup


# Recommendations from https://packaging.python.org/
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def read(*parts):
    with open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='e3nn',
    version=1.1,
    description='Hidden Markov state model analysis for vesicle simulations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ka-zuumi/vesicle-HMM-analysis',
    packages=find_packages(exclude=["tests.*", "tests"]),
    install_requires=[
        'numpy',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pre-commit',
        ],
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
    ],
    python_requires='>=3.7',
    license="MIT",
    license_files="LICENSE",
)
