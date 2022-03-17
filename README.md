# imagine-models
A model library usable for Galactic inference engines. 
The version on this branch is focused on wrapping hammurabi x models. 
To this date, only the regular JF12 model  has been wrapped.

## Installation

First clone the branch via

    git clone -b wrapping_hammurabi_models --recursive https://github.com/IMAGINE-Consortium/imagine-models.git

The `--recursive` flag makes sure that also the pybind module is cloned (to the `/extern` folder) 

The package can then be installed with 

    pip3 install folder/where/setup/py/is/

Note that this mode of installation is still preliminary and might change in the future.
Also, directly installing via git (by e.g. `pip install git+...` ) is not yet supported.

