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

def compute_linear_model(mfs, measures):
    from sklearn.linear_model import Ridge
    from sklearn import linear_model

    # try different ones
    clf = Ridge(alpha = 1.0)
    #clf = RidgeCV(alphas=[0.1, 1.0, 10.0])
    #clf = linear_model.LinearRegression()

    # explain fexp using BMD + the MFS data
    fexp = measures[:, measures.shape[1]-1]

    bmd = measures[:, 0]
    bmd = bmd.reshape((bmd.shape[0], 1))

    print "BMD: ", bmd
    print "FEXP: ", fexp
    print "MFS; ", mfs

    #PCA
    #from sklearn.decomposition import PCA
    #pca = PCA(n_components=12)
    #pca.fit(mfs)
    #mfs_pca = pca.transform(mfs)

    X = np.hstack((bmd, mfs))
    clf.fit(X, fexp)

    # Results
    print "Coefs:", clf.coef_
    print "Score (R^2):", clf.score(X, fexp)

def do_test():

    dims = 10

    MFS_HOLDER = True
    APPLY_LAPLACIAN = False
    APPLY_GRADIENT = True

    TRANSFORMED_INPUT_STR = ''
    MFS_STR = ''

    data_path = "exps/data/"

    if MFS_HOLDER:
        MFS_STR = '_holder'

    if APPLY_LAPLACIAN:
        TRANSFORMED_INPUT_STR = '_laplacian'
    else:
        if APPLY_GRADIENT:
            TRANSFORMED_INPUT_STR = '_gradient'

    BASE_NAME = 'mfs' + MFS_STR + TRANSFORMED_INPUT_STR + '_BioAsset'

    data = np.load(data_path + BASE_NAME + '_and_standard_params_adaptive.npy')

    mfs_last_d = 20
    if MFS_HOLDER:
        mfs_last_d = 19

    mfs = data[:, : mfs_last_d]
    measures = data[:, mfs_last_d+1 :]


    fsize = 18
    e = 0.1
    s = 10.0
    x1 = 1
    x2 = 20
    x = np.arange(mfs_last_d)+1


    print "Measures.shape : ", measures.shape
    print "mfss.shape : ", mfs.shape

    # one correlation for each measure and each mfs dimension
    correls = np.zeros((mfs.shape[1], measures.shape[1]))

    which_d = 0
    print "MFS0", mfs[:, which_d]
    print "measures",measures[:, which_d]

    #plt.plot(mfs[:, which_d])
    #plt.plot(measures[:, which_d])
    #plt.show()

    #print "CORR", scipy.stats.stats.spearmanr(mfs[:, 15], measures[:, 15])[0]

    for d in range(mfs.shape[1]):
        for m in range(measures.shape[1]):
            correls[d, m] = scipy.stats.stats.spearmanr(mfs[:, d], measures[:, m])[0]

    print correls
    print "Higher correlations: ", np.min(correls), np.max(correls)

    #plt.plot(correls)
    #plt.show()

    np.save(data_path + "correls_measures_adaptive_mfs" + BASE_NAME + ".npy", correls)

    # compute linear models (BMD + fractal measures)
    compute_linear_model(mfs, measures)