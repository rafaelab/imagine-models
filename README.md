# Imagine Models
A model library usable for Galactic inference engines.
The library is (mostly) written C++, but can be accessed via python.
Some (non-essential) elements are written in pure python, which presently cannot be accessed in the c-library.

## Installation

### Python

#### Requirements

Requires:

    [Python 3](https://www.python.org/) (>3.6)


    [NumPy](https://numpy.org/) (>1.22)

Optional:



#### Installation procedure

Easiest via

    pip3 install --user git+https://github.com/IMAGINE-Consortium/imagine-models.git

If you install the library in a virtual environment, remove the --user tag.


If you want to specify a branch, you can do so by adding @branch-name to the above command.

If you want to add your own model you need to clone the repository via

    git clone -b branch_name --recursive https://github.com/IMAGINE-Consortium/imagine-models.git

The `--recursive` flag makes sure that also the pybind module is cloned (to the `/extern` folder).

The package can then be installed with

    pip3 install folder/where/setup/py/is/


### C++

TBD
