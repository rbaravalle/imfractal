imfractal
=======

An open source (BSD) library to compute (multi)fractal dimensions of images written in Python. Different implementations of the multifractal formalism are present, e.g., the Sandbox method, the MFS and the Singularity multifractal spectrum. The support for the use of the pyopencl library is in process of development.

## Required dependencies for imfractal:

* Python 2.7
* Numpy 1.1 or higher       (http://numpy.scipy.org/)   
* Scipy 0.7 or higher       (http://www.scipy.org/)
* PIL 1.1.7 or higher       (http://www.pythonware.com/products/pil/)
* scikit-learn 0.10.0 or higher (http://scikit-learn.org/)

## Installation

* In Debian/Ubuntu, just install the required Python packages:


    apt-get install python-numpy python-scipy python-sklearn

## Usage

The usage of imfractal is shown using test.py

## Multifractal Spectrums implemented (so far)

1. Singularity Multifractal spectrum (Haussdorf dimentions through Hölder exponents)
2. Sandbox Multifractal spectrum
3. Multifractal spectrum (MFS)
