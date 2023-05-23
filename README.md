Code for computing the relic neutrino density with kinetic field theory.

## Installing the Python wrapper

You need to have the `cython` package for Python installed. How to do this depends on your installation; if you are using Anaconda, do `conda install cython`.

Go to the directory `python/` and execute the command `python setup.py install`.

You should now be able to import the package `kftneutrinos` in Python. 

## The `kftneutrinos` package

The Python package `kftneutrinos` gives access to the following functions. For more details on each computation method, see the paper `XXXX.XXXX`.

* `py_first_order`: Computes the neutrino clustering factor to first order in kinetic field theory perturbation theory. For an example of how to call the function, see `/notebooks/first_order.ipynb`.
