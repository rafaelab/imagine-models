# IMAGINE Model Library

A model library usable for Galactic inference engines. 
The library is (mostly) written C++, but can be accessed via both C++ and Python.
Some (non-essential) elements are written in pure Python, which presently cannot be accessed in the C++ version of the library.
A full list of implemented models can be found [here](#list-of-models).

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

To install the C++ library (Debian only until now)
 
```
cd ./c_library  #important
mkdir build
cd build 
cmake ..
sudo make install 
```

Other systems have not been tested yet.

Note that we build within the `c_library` directory to avoid with the build of the Python package, which is 

## Examples

### Including the pipeline

Example scripts demonstrating how to include both the python and C++ version are located in the ./demo folder. 
In the python case, we also include a Jupyter notebook. 
 
### Adding new models **TBD** 

Defining your own models is easy in both Python and C++. 
For that you can start with the model templates provided in the `\templates` folder. 
Binding a C++ model to Python is a bit more involved. For the simplest case, a template exists as well. 


## List of Models

| MODEL NAME  | PYTHON      | C++         | reference   | notes       |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| **Regular models** |           |             |             |             |
| Uniform | U+2713 | U+2713 |             |             |
| Helix | U+2713 | U+2713 |             |             |
| Axissymetric spiral | U+2713 | U+2717 |              |             |
| Jaffe | U+2713 | U+2713 |             |             |
| Jansson Farrar | U+2713 | U+2713 |             |             |
| **Random models** |             |             |             |             |
| Jansson Farrar | U+2713 | U+2713 |             |             |