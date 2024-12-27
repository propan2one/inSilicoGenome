# InSilicoGenome 

[![Python](https://img.shields.io/badge/python-3.9-blue)]()
[![codecov](https://codecov.io/gh/propan2one/insilicogenome/branch/main/graph/badge.svg)](https://codecov.io/gh/propan2one/insilicogenome)
[![Documentation Status](https://readthedocs.org/projects/insilicogenome/badge/?version=latest)](https://insilicogenome.readthedocs.io/en/latest/?badge=latest)


## Installation

### Using pip

```bash
$ pip install -i https://test.pypi.org/simple/ insilicogenome
```

### Using github clone conda envs with poetry


```bash
# 1) Clone the repo 
git clone git@github.com:propan2one/inSilicoGenome.git
 
# 2) create a conda env where all tools work together
#    in my case, for ease of use, I'm going to use biopython
conda create -y -p ~/envs/insilicogenome \
    --channel conda-forge python=3.11.11 Poetry
conda activate ~/envs/insilicogenome
poetry install
```

## Features

- TODO

## Dependencies

- TODO

## Usage

- TODO

## Documentation

The official documentation is hosted on Read the Docs: https://insilicogenome.readthedocs.io/en/latest/

## Contributors

We welcome and recognize all contributions. You can see a list of current contributors in the [contributors tab](https://github.com/propan2one/insilicogenome/graphs/contributors).

### Credits

This package was created with Cookiecutter and the UBC-MDS/cookiecutter-ubc-mds project template, modified from the [pyOpenSci/cookiecutter-pyopensci](https://github.com/pyOpenSci/cookiecutter-pyopensci) project template and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage).
