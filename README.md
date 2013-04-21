imfractal
=======

An open source (BSD) library to compute (multi)fractal dimensions of images written in Python. Different implementations of the multifractal formalism are present, e.g., the Sandbox method, the MFS and the Singularity multifractal spectrum. The support for the use of the pyopencl library is in process of development.

## Required dependencies for imfractal:

* Python 2.7
* Numpy 1.1 or higher       (http://numpy.scipy.org/)   
* Scipy 0.7 or higher       (http://www.scipy.org/)
* PIL 1.1.7 or higher       (http://www.pythonware.com/products/pil/)

Optionaly, the classification test requires:

* scikit-learn 0.10.0 or higher (http://scikit-learn.org/)

## Installation

* In Debian/Ubuntu, just install the required Python packages:


    apt-get install python-numpy python-scipy python-sklearn

## Example

The most interesting example is a classifier to distiguish [bread images](https://github.com/rbaravalle/imfractal/tree/master/images/test/bread) from [others types of images](https://github.com/rbaravalle/imfractal/tree/master/images/test/nonbread) using imfractal. The code is available [here](https://github.com/rbaravalle/imfractal/blob/master/tests/test_classifier.py).

## Multifractal Spectrums implemented (so far)

1. Singularity Multifractal spectrum (Haussdorf dimentions through Hölder exponents)
2. Sandbox Multifractal spectrum
3. Multifractal spectrum (MFS)
