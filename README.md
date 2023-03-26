# Imagine Models
A model library usable for Galactic inference engines. 
The library is (mostly) written C++, but can be accessed via both C++ and Python.
Some (non-essential) elements are written in pure Python, which presently cannot be accessed in the C++ version of the library.


**TODO: Add model reference at bottom and link it here** 

## Installation (Python)

#### Requirements

**TODO: check fftw3 dependence, make it optional** 

Requires:

- [Python 3](https://www.python.org/) (>=3.6)
- [C++](https://www.python.org/) (>=11)

Python Libraries:

- [NumPy](https://numpy.org/) (>=1.22)

Optional (Developers):

- [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html)

Note that pybind is included in the current repository, so in the current state you should not need to install it.
This may change in the future.  


#### Installation procedure

Easiest via

    pip3 install --user git+https://github.com/IMAGINE-Consortium/imagine-models.git

If you install the library in a virtual environment, remove the --user tag


If you want to specify a branch, you can do so by adding @branch-name to the above command.

Alternatively, (e.g. if you want to add your own model) you need to clone the repository via

    git clone -b branch_name --recursive https://github.com/IMAGINE-Consortium/imagine-models.git

The `--recursive` flag makes sure that also the pybind module is cloned (to the `/extern` folder).

The package can then be installed with

    pip3 install folder/where/setup/py/is/


## Installation (C++)

From 


Note that we build within the 

## Examples

#### Including the pipeline

Example scripts demonstrating how to include both the python and C++ version are located in the ./examples folder. 
The python
 
#### Adding new models

Defining your own models is easy in both Python and C++, if you want to 

**TODO: model templates** 