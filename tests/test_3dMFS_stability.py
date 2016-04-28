"""
Copyright (c) 2016 Rodrigo Baravalle
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
"""
from imfractal import *

from pylab import *

import sys
import os
import numpy as np

sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'imfractal', 'imfractal',  "Algorithm"))

import qs3D


def do_test(_path, N):

    print "PATH: " + _path

    patients = ["32"]
    scans = ["1"]
    # amount of volumes of interest
    vois = 1
    dims = 10

    # BioAsset bone's multifractal spectra database
    mfss = np.zeros([len(patients),len(scans),vois,2*dims+1])

    aux = CSandbox3D(dims)

    params = {
        "zero": 1,
        "one": 0.75,
        "two": 3.7,
        "three": 1,
        "four": 15,
        "five": 0,
        "mask_filename": '',
        "seven": "no",
        "eight": 'S',
        "nine": 'M',
        "threshold": 200
    }

    fmask = _path + "BA" + patients[0] + "_120_" + scans[0] + "Mask.mat"

    params["five"] = 0
    params["mask_filename"] = fmask

    aux.setDef(40, 1.02, True, params)

    slice_filename = _path + "BA" + patients[0] + "_120_" + scans[0] + "Slices.mat"

    print fmask
    print slice_filename

    if N > 1 :
        print "Repeating", N, " times"

        mfss = np.array(aux.getFDs(slice_filename)).astype(np.double)
        for i in range(1, int(N)):
            print "Computing ", i, " th time"
            mfss = np.vstack((mfss, aux.getFDs(slice_filename)))

        print "MFSS SHAPE: ", mfss.shape

        # show variations
        print "Standard Deviation by dimension: "
        for j in range(0, mfss.shape[1]):
            print np.std(mfss[:, j]), " j: ", j


    
