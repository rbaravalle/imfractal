import numpy as np
from imfractal import *
from numpy import recfromcsv
import scipy
import math

# for Fexp
#fexp_names = np.load(data_path + 'bioAsset_meta.npy')

# Adaptive Metadata and mfs
measures = recfromcsv(data_path + 'default_BioAsset_Adaptive.csv', delimiter=',')
mfs = recfromcsv('exps/data/mfs_holder_BioAsset.csv', delimiter=',')
mfs_normalized = recfromcsv('exps/data/BioAssetAdaptiveThresh/mfs_holder_BioAsset.csv', delimiter=',')
mfs_sandbox_adaptive = np.load('exps/data/BioAssetAdaptiveThresh/mfs_Sandbox_BioAsset_adaptive_0.75.npy')

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

    #print "BMD: ", bmd
    #print "FEXP: ", fexp
    #print "MFS; ", mfs

    #PCA
    #from sklearn.decomposition import PCA
    #pca = PCA(n_components=12)
    #pca.fit(mfs)
    #mfs_pca = pca.transform(mfs)

    X = np.hstack((bmd, mfs))
    clf.fit(X, fexp)

    # Results
    #print "Coefs:", clf.coef_
    print "Score (R^2):", clf.score(X, fexp)

def compute_correlations(measures_matrix, mfs, mfs_pos_start_data,
                                        mfs_pos_end_data):
    print mfs.shape

    result = []

    # convert from ugly format to matrix
    mfs_matrix = []

    for i in range(mfs.shape[0]):

        mfs_i = tuple(mfs[i])
        mfs_i = mfs_i[mfs_pos_start_data : mfs_pos_end_data + 1]


        if len(mfs_matrix) == 0:
            mfs_matrix = mfs_i
        else:
            mfs_matrix = np.vstack((mfs_matrix, mfs_i))


    correls = np.zeros((mfs_matrix.shape[1], measures_matrix.shape[1]))
    # compute correlations
    for d in range(mfs_matrix.shape[1]):
        for m in range(measures_matrix.shape[1]):
            correls[d, m] = scipy.stats.stats.spearmanr(mfs_matrix[:, d],
                                                        measures_matrix[:, m])[0]

    #DEBUG CORRELATIONS:
    #import matplotlib.pyplot as plt
    #plt.plot(correls)
    #plt.show()
    print "Higher correlations: ", np.min(correls), np.max(correls)

def compute_subset(measures_matrix, mfs,
                   mfs_pos_start_data, mfs_pos_end_data):
    mfs_subset = []
    for i in range(len(measures)):
        if not (math.isnan(measures_matrix[i][pos_fexp])):

            mfs_i = tuple(mfs[i])
            mfs_i = mfs_i[mfs_pos_start_data: mfs_pos_end_data]
            if len(mfs_subset) == 0:
                mfs_subset = np.array(mfs_i)
            else:
                mfs_subset = np.vstack((mfs_subset, mfs_i))

    return mfs_subset


measures_pos_start_data = 1
measures_pos_end_data = 18

measures_matrix = []

for i in range(measures.shape[0]):
    measures_i = tuple(measures[i])
    measures_i = measures_i[measures_pos_start_data: measures_pos_end_data + 1]

    if len(measures_matrix) == 0:
        measures_matrix = np.array(measures_i)
    else:
        measures_matrix = np.vstack((measures_matrix, measures_i))

print "MEASURES MATRIX: ", measures_matrix.shape

################################################################

mfs_pos_start_data = 1
mfs_pos_end_data = 20



pos_fexp = 17 #check


#print "Correlations with raw MFS..."
#compute_correlations(measures, mfs)
print "Correlations with normalized MFS..."
compute_correlations(measures_matrix, mfs_normalized, mfs_pos_start_data,
                                        mfs_pos_end_data)


# obtain subsets of 17 scans for Fexp
mfs_subset = np.array([])
measures_subset = np.array([])

i = 0

# subset of measures
for i in range(len(measures)):
    if not(math.isnan(measures_matrix[i][pos_fexp])):
        if len(measures_subset) == 0:
            measures_subset = np.array(measures_matrix[i])
        else:
            measures_subset = np.vstack((measures_subset, measures_matrix[i]))

mfs_subset = compute_subset(measures_matrix, mfs_normalized, mfs_pos_start_data,
                            mfs_pos_end_data)

print "MFS_SUBSET: ", mfs_subset.shape
print "MEASURES_SUBSET: ", measures_subset.shape

compute_linear_model(mfs_subset, measures_subset)


mfs_pos_start_data = 0
mfs_pos_end_data = 21

print "Correlations with Sandbox MFS Adaptive..."
compute_correlations(measures_matrix, mfs_sandbox_adaptive, mfs_pos_start_data,
                                        mfs_pos_end_data)


mfs_sandbox_subset = compute_subset(measures_matrix, mfs_sandbox_adaptive,
                                    mfs_pos_start_data, mfs_pos_end_data)

print "MFS_SUBSET: ", mfs_subset.shape
print "MEASURES_SUBSET: ", measures_subset.shape

compute_linear_model(mfs_sandbox_subset, measures_subset)



