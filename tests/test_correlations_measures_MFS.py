"""
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
"""
from imfractal import *

import Image
import time
import matplotlib.pyplot as plt
from pylab import *

import os
import scipy

def createFolders(dirs):

    for f in dirs:
        if not os.path.isdir(f): 
            os.mkdir (f) 

def do_test():

    dims = 10

    data = np.load('mfs_and_standard_params.npy')

    mfs = data[:, : 2 * dims + 1]
    measures = data[:, 2 * dims + 2:]


    fsize = 18
    e = 0.1
    s = 10.0
    x1 = 1
    x2 = 20
    x = np.arange(2*dims+1)+1


    print "Measures.shape : ", measures.shape
    print "mfss.shape : ", mfs.shape

    # one correlation for each measure and each mfs dimention
    correls = np.zeros((mfs.shape[1], measures.shape[1]))

    print "MFS0", mfs[:, 13]
    print "measures",measures[:, 13]

    plt.plot(mfs[:, 13])
    plt.plot(measures[:, 13])
    plt.show()

    #print "CORR", scipy.stats.stats.spearmanr(mfs[:, 15], measures[:, 15])[0]



    for d in range(mfs.shape[1]):
        for m in range(measures.shape[1]):
            correls[d, m] = scipy.stats.stats.spearmanr(mfs[:, d], measures[:, m])[0]

    print correls
    print "Higher correlations: ", np.min(correls), np.max(correls)

    plt.plot(correls)
    plt.show()

    np.save("correls_measures_mfs.npy", correls)


