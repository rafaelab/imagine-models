# IMAGINE Model Library

A model library usable for Galactic inference engines. 
The library is (mostly) written C++, but can be accessed via both C++ and Python.
Some (non-essential) elements are written in pure Python, which presently cannot be accessed in the C++ version of the library.
A full list of implemented models can be found [here](#list-of-models).

<span style="color:red;font-weight:700;font-size:18px">
    Please note that this library is work in progress and in early development. 
    By far not everything has been tested, and interfaces may still change.
    Don't hesitate to use this software, but beware of the risks.
</span>

## Requirements

Requires:

- [Python 3](https://www.python.org/) (>=3.6)
- [C++](https://www.python.org/) (>=11)

Python Libraries:

- [NumPy](https://numpy.org/) (>=1.22)

Optional:

- [FFTW3](http://fftw.org/) (>3.3, necessary for random models)
- [autodiff](https://autodiff.github.io/) (>1.0, necessary for building gradients w.r.t. to model parameters)

Optional (Developers):

- [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html)

Note that pybind is included in the current repository, so in the current state you should not need to install it.
This may change in the future.  


## Installation 

The installation procedures are preliminary only 


### Installation (C++)

To install the C++ library (tested only under Debian until now):
 
```
mkdir build
cd build 
cmake ..
sudo make install 
```


### Installation (Python)

Easiest via

    pip3 install --user git+https://github.com/IMAGINE-Consortium/imagine-models.git

If you install the library in a virtual environment, remove the `--user` tag.

If you want to specify a branch, you can do so by adding @branch-name to the above command.

Alternatively, (e.g. if you want to add your own model) you need to clone the repository via

    git clone -b branch_name --recursive https://github.com/IMAGINE-Consortium/imagine-models.git

The `--recursive` flag makes sure that also the pybind module is cloned (to the `/extern` folder).

The package can then be installed with

    python3 -m pip install folder/where/setup/py/is/

The installer will automatically figure out which of the optional dependencies you have installed. 

If you want to disable those, you can do so by defining the `USE_AUTODIFF=OFF` and `USE_FFTW=OFF` environment variables BEFORE you run pip. 



## Examples

### Including the library

Example scripts demonstrating how to include both the python and C++ version are located in the ./demo folder. 
In the python case, we also include Jupyter notebooks. 

### Adding new models **TBD** 

Defining your own models is easy in both Python and C++. 
For that you can start with the model templates provided in the `\templates` folder. 
Binding a C++ model to Python is a bit more involved. For the simplest case, a template exists as well. 


## List of Models

### Magnetic/Vector Fields

| MODEL NAME          | PYTHON      | C++         | reference    | notes               | original implementation |
| -----------         | ----------- | ----------- | -----------  | -----------         | -----------             | 
| **Regular models**  |             |             |              |                     |                         |
| Uniform             | &#x2714;    | &#x2714;    |              | used for unit tests |                         |
| Helix               | &#x2714;    | &#x2714;    |              |                     |                         |
| Axissymetric spiral | &#x2714;    | &#x2718;    | Pelgrims, V. |                     |  Pelgrims, V.           |
| Archimedean spiral | &#x2714;    | &#x2714;    |  |                     |   [CRPRopa](https://github.com/CRPropa/CRPropa3/tree/master)        |
| Local Bubble | &#x2714;    | &#x2718;    | [Pelgrims et al.](https://www.aanda.org/articles/aa/full_html/2020/04/aa37157-19/aa37157-19.html) |   Only defined on the shell                   |  Pelgrims, V.           |
| Jaffe               | &#x2714;    | &#x2714;    | [Jaffe et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010MNRAS.401.1013J/abstract)        |             | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)  |
| Sun2008               | &#x2714;    | &#x2714;    | [Sun et al. (2008)](https://www.aanda.org/articles/aa/abs/2008/02/aa8671-07/aa8671-07.html)        |             | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  |
| Han2018               | &#x2714;    | &#x2714;    | [Han et al. (2018)](https://iopscience.iop.org/article/10.3847/1538-4365/aa9c45)        |             |  |
| Pshirkov | &#x2714;    | &#x2714;    | [Pshirkov et al.](https://iopscience.iop.org/article/10.1088/0004-637X/738/2/192)  |                     |   [CRPRopa](https://github.com/CRPropa/CRPropa3/tree/master)        |
| HMR              | &#x2714;    | &#x2714;    | [Harari et al. (1999)](https://arxiv.org/abs/astro-ph/9906309)        |             | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) and [Kachelrieß (2007)](https://arxiv.org/pdf/astro-ph/0510444.pdf )  |
| TT              | &#x2714;    | &#x2714;    | [Tinyakov and Tkachev (2017)](https://arxiv.org/abs/astro-ph/0111305)        |             | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) and [Kachelrieß (2007)](https://arxiv.org/pdf/astro-ph/0510444.pdf )  |
| TF              | &#x2714;    | &#x2714;    | [Terral and Ferriere (2017)](https://arxiv.org/abs/1611.10222)        |             | [CRPRopa](https://github.com/CRPropa/CRPropa3/tree/master)        | |
| Fauvet              | &#x2714;    | &#x2714;    | [Fauvet et al. (2012)](https://arxiv.org/abs/1201.5742)        |             | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  |
| Stanev            | &#x2714;    | &#x2714;    | [Stanev (1996)](https://arxiv.org/abs/astro-ph/9607086)        |             | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  |
| WMAP              | &#x2714;    | &#x2714;    | [Page et al. (2006)](https://arxiv.org/pdf/astro-ph/0603450.pdf)        |             | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  |
| Jansson Farrar      | &#x2714;    | &#x2714;    | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...757...14J/abstract)        |             | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)   |
| **Random models**   |             |             |              |                     |                         |
| Jansson Farrar      | &#x2714;    |&#x2714;     | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...761L..11J/abstract)        | depends on the JF12 regular model   |[Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)    |
| Ensslin Steininger  | &#x2714;    |&#x2714;     |              |                      | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)             |

### Thermal electron/Scalar Fields

| MODEL NAME  | PYTHON      | C++         | reference   | notes       |  original implementation |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| **Regular models**        |             |             |             |             |             |
| YMW16       | &#x2714;    | &#x2714;    | [Yao et al. (2016)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract)     |             |       [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)         |
| **Random models** |             |             |             |             |             |
| GaussianScalar | &#x2714;|&#x2714; |             |      used for random number unit testing (TBD!)       |             |
| LogNormalScalar | &#x2714;|&#x2714; |             |             |             |
