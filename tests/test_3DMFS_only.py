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
from imfractal import MFS_3D as MFS
import matplotlib.pyplot as plt

import time
import numpy as np
import math

def checkerboard(w, h, c0, c1, blocksize):
        tile = np.array([[c0,c1],[c1,c0]]).repeat(blocksize, axis=0).repeat(blocksize, axis=1)
        grid = np.tile(tile,(int(math.ceil((h+0.0)/(2*blocksize))),int(math.ceil((w+0.0)/(2*blocksize)))))
        return grid[:h,:w]

def checker_3d(size, width):
    res = np.zeros((size, size, size))

    res[:,:,size-1] = checkerboard(size, size, 0, 1, width)

    return res

def do_test():

    # construct checkerboard
    checker = checker_3d(50, 25)

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
        "threshold": 200,
        "total_pixels":6000,
        "adaptive" : False, # adaptive threshold (only for not holder)
        "laplacian": False,
        "gradient" : False
    }

    i = MFS()
    i.setDef(1,20,3,'', '', params)
    print "Computing Lipschitz-Holder 3D Multifractal Spectrum... (monofractal checkerboard)"
    t =  time.clock()
    fds = i.getFDs(checker)
    t =  time.clock()-t
    print "Time 3D MFS: ", t
    print fds

    plt.title('Monofractal')
    plt.ylim((0.0, 3.2))
    plt.plot(fds, 'x', label = 'Mono')
    
    multif = np.load('exps/data/img3d.npy')
    print "Computing Lipschitz-Holder 3D Multifractal Spectrum... (multifractal)"
    t =  time.clock()
    fds = i.getFDs(multif)
    t =  time.clock()-t
    print "Time 3D MFS: ", t
    print fds

    plt.plot(fds,'-', label = 'Multi')
    plt.legend()
    plt.show()
