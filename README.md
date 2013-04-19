--------------------------------------------------------------------
Copyright (c) 2013 Rodrigo Baravalle
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

imfractal
=======

An open source library to compute (multi)fractal dimensions of images written in Python. Different implementations of the multifractal formalism are present, e.g., the Sandbox method, the MFS and the Singularity multifractal spectrum. The support for the use of the pyopencl library is in process of development.

## Required dependencies for imfractal:

* Python 2.7
* Numpy 1.1 or higher       (http://numpy.scipy.org/)   
* Scipy 0.7 or higher       (http://www.scipy.org/)
* PIL 1.1.7 or higher       (http://www.pythonware.com/products/pil/)

## Installation

* In Debian/Ubuntu, just install the required Python packages:


    apt-get install python-numpy python-scipy python-image 

## Usage

The usage of imfractal is shown using test.py

## Multifractal Spectrums implemented (so far)

1. Singularity Multifractal spectrum (Haussdorf dimentions through Hölder exponents)
2. Sandbox Multifractal spectrum
3. Multifractal spectrum (MFS)
