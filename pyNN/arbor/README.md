# Introduction

This is the Arbor backend for PyNN.

# Wrapping Strategy

PyNN has it's own way of describing and building a network of neurons that needs
to be mapped to Arbor's way of doing this. A model description in Arbor is 
called a *recipe*. This recipe is a class that contains all information regarding
the cell types and their interconnections. The underlying C++ recipe class is 
mapped to an equivalent Python class using *PyBind11* (see 
[`arbor/python/recipe.hpp`](https://github.com/bcumming/arbor/blob/feature/pybind11/python/recipe.hpp). An example of an implemented recipe can found [here](https://github.com/bcumming/arbor/blob/feature/pybind11/python/example/ring.py).

To map a PyNN network description in terms of Populations and Projections to
an Arbor recipe, we can use a singleton `Recipe` class similar to how PyNN's
NEURON backend uses [`_Initializer` or `_State`](https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/simulator.py) that registers and initializes PyNN Populations.

# Installation

Temporary installation instructions for unmerged Arbor Python interface:

- Clone and build the unmerged Python feature branch of Arbor:

```sh
git clone https://github.com/bcumming/arbor.git
cd arbor
git checkout feature/pybind11

# make sure Python headers of desired version can be found
conda create -n arb-dev python=<version>
source activate arb-dev
conda_home="$(conda info --base)"
export PATH=$PATH:$conda_home
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$conda_home/lib

# Build arbor with Python interface
mkdir build && cd build
cmake .. -DARB_WITH_PYTHON=ON
make -j 4
```

- Put the `pyarb` module on your Python path

```sh
cd arbor/build/lib
ln -s pyarb.cpython-<version>.so pyarb # symlink pyarb Python module
# use any method that puts pyarb module on Python path
export PYTHONPATH=$PYTHONPATH:`pwd`
```