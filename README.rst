.. image:: https://travis-ci.org/fabiommendes/euclid.svg?branch=master
    :target: https://travis-ci.org/fabiommendes/euclid


`euclid` is a lightweight library that implements linear algebra operations
in low dimensions. These objects were create to be used in a game engine, but
may be useful elsewhere. Under the hood, it wraps GLM and Eigen, but does not
try to reproduce the API of neither of those libraries.

The library is optimized using C-extensions written in Cython and it is one of
the fastest Python library for linear algebra and shape manipulation in low
dimension.

Euclid supports:

* Vectors and directions: 2d


Installation
============

You just need pip and Python 3.6+

    $ python3 -m pip install euclid-math

