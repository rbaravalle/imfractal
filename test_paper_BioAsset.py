import numpy as np
from imfractal import *
from numpy import recfromcsv
import scipy
import math
import pandas
from pandas import DataFrame, Series
import statsmodels.formula.api as sm
import statsmodels
from statsmodels.tools.eval_measures import aicc, bic, hqic

np.seterr(divide='ignore', invalid='ignore')

# for Fexp
#fexp_names = np.load(data_path + 'bioAsset_meta.npy')

# Adaptive Metadata and mfs
measures = recfromcsv('exps/data/BioAssetAdaptiveThresh/default_BioAsset_Adaptive.csv', delimiter=',')
mfs = np.load('exps/data/mfs_holder_BioAsset_raw.npy')
mfs_normalized = recfromcsv('exps/data/BioAssetAdaptiveThresh/mfs_holder_BioAsset.csv', delimiter=',')
mfs_sandbox_adaptive = np.load('exps/data/BioAssetAdaptiveThresh/mfs_Sandbox_BioAsset_adaptive_0.75.npy')
# should be similar to mfs_sandbox_adaptive?
mfs_sandbox_absolute_normalized = np.load('exps/data/mfs_Sandbox_BioAsset_normalized.npy')
mfs_local = np.load('exps/data/mfs_holder_local_BioAsset.npy')
mfs_local_pyramid = np.load('exps/data/mfs_holder_local_BioAsset_pyramid.npy')
mfs_gradient_pyramid = np.load('exps/data/mfs_holder_pyramid_gradient_BioAsset.npy')
mfs_slices_x = np.load('exps/data/mfs_holder_slices_x_BioAsset.npy')
mfs_slices_y = np.load('exps/data/mfs_holder_slices_y_BioAsset.npy')
mfs_slices_z = np.load('exps/data/mfs_holder_slices_z_BioAsset.npy')
mfs_pure_pyramid = np.load('exps/data/mfs_pure_pyramid_BioAsset.npy')
mfs_pure_pyramid_gradient = np.load('exps/data/mfs_pure_pyramid_gradient_BioAsset.npy')
mfs_sigmoid = np.load('exps/data/mfs_holder_sigmoid_BioAsset.npy')
mfs_holder_10 = np.load('exps/data/mfs_holder_10_BioAsset.npy')
mfs_holder_5 = np.load('exps/data/mfs_holder_5_BioAsset.npy')
stats_mfs_holder = np.load('exps/data/stats_mfs_holder_BioAsset.npy')


pos_fexp = 17 #check

def normalize(vector):
    if np.std(vector) != 0:
        return (vector - np.mean(vector))/ np.std(vector)
    else:
        print "std equals 0"
        return vector

def compute_best_aicc(X, fexp):
    X2 = statsmodels.tools.tools.add_constant(X)

    # one dimension

    best_aicc = 100000
    best_aicc2 = 100000
    best_aicc3 = 100000
    best_i = -1
    best_i_j = [-1, -1]
    best_i_j_k = [-1, -1, -1]
    best_r2_ij = -1
    best_r2 = -1
    best_r2_ijk = -1


    for i in range(1, X.shape[1]):
        Xi = X2[:, [0, 1, i]]

        model = sm.OLS(fexp, Xi)
        res = model.fit()

        aic = aicc(res.llf, res.nobs, res.params.shape[0])

        if aic < best_aicc :
            best_aicc = aic
            best_i = i
            best_r2 = res.rsquared_adj

    for i in range(1, X.shape[1]):
        for j in range(i+1, X.shape[1]):
            Xij = X2[:, [0, 1, i, j]]


            model = sm.OLS(fexp, Xij)
            res = model.fit()

            aic = aicc(res.llf, res.nobs, res.params.shape[0])

            if aic < best_aicc2:
                best_aicc2 = aic
                best_i_j = [i, j]
                best_r2_ij = res.rsquared_adj


    for i in range(1, X.shape[1]):
        for j in range(i+1, X.shape[1]):
            for k in range(j + 1, X.shape[1]):
                Xijk = X2[:, [0, 1, i, j, k]]


                model = sm.OLS(fexp, Xijk)
                res = model.fit()

                aic = aicc(res.llf, res.nobs, res.params.shape[0])

                if aic < best_aicc3:
                    best_aicc3 = aic
                    best_i_j_k = [i, j, k]
                    best_r2_ijk = res.rsquared_adj

    return best_aicc, best_i, best_r2,\
           best_aicc2, best_i_j, best_r2_ij,\
           best_aicc3, best_i_j_k, best_r2_ijk

def compute_linear_model(mfs, measures, output_file="standarized.csv"):
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
    # print "Coefs:", clf.coef_
    #print ""

    #print "Score (R^2):", clf.score(X, fexp)

    #cols = ['bmd']
    #for i in range(mfs.shape[1]):
        #cols.append('mfs_' + str(i))

    #### using statsmodel
    #df = DataFrame(np.hstack((X, np.array([fexp]).T)), columns=cols)
    #df = DataFrame(X, columns=cols)

    # BMD ALONE


    Xbmd = X[:, [0]]
    X2 = statsmodels.tools.tools.add_constant(Xbmd)


    model = sm.OLS(fexp, X2)
    res = model.fit()

    aic = aicc(res.llf, res.nobs, res.params.shape[0])
    print "BMD AICc, dimension, R2: " , aic, ' bmd ', res.rsquared_adj

    res =  compute_best_aicc(X, fexp)
    print "AICc, dimension, R2: ", res[0], ' bmd + ', res[1], res[2]
    print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[0])/2.0)
    print "AICc, dimensions, R2: ", res[3],' bmd + ',  res[4], res[5]
    print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[3]) / 2.0)
    print "AICc, dimensions, R2: ", res[6],' bmd + ',  res[7], res[8]
    print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[6]) / 2.0)

    return
    X_normalized = X
    for i in range(X.shape[1]):
        X_normalized[:,i] = normalize(X_normalized[:,i])

    clf.fit(X_normalized, fexp)


    print ""

    print "Normalized Variables - Score (R^2):", clf.score(X_normalized, fexp)

    X_2 = statsmodels.tools.tools.add_constant(X_normalized)
    model_2 = sm.OLS(fexp, X_2)
    res_2 = model_2.fit()

    print "R^2 adj", res.rsquared_adj
    print "AICc: ", aicc(res_2.llf, res_2.nobs, res_2.params.shape[0])
    #print "AIC: ", res.aic

    np.savetxt(data_path + output_file, X_normalized, delimiter=",")

def compute_correlations(measures_matrix, mfs, mfs_pos_start_data,
                                        mfs_pos_end_data, transform_mfs = True):
    print "///////////////////////"

    result = []

    # convert from ugly format to matrix
    mfs_matrix = []

    if transform_mfs:
        for i in range(mfs.shape[0]):
            mfs_i = tuple(mfs[i])
            mfs_i = mfs_i[mfs_pos_start_data : mfs_pos_end_data + 1]


            if len(mfs_matrix) == 0:
               mfs_matrix = mfs_i
            else:
               mfs_matrix = np.vstack((mfs_matrix, mfs_i))

    else: mfs_matrix = mfs



    correls = np.zeros((mfs_matrix.shape[1], measures_matrix.shape[1]))
    # compute correlations
    for d in range(mfs_matrix.shape[1]):
        for m in range(measures_matrix.shape[1]):
            corr = scipy.stats.stats.spearmanr(mfs_matrix[:, d],
                                                        measures_matrix[:, m])[0]
            if not(math.isnan(corr)):
                correls[d, m] = corr
            else:
                correls[d, m] = 0

    #DEBUG CORRELATIONS:
    #import matplotlib.pyplot as plt
    #plt.plot(correls)
    #plt.show()
    print "Higher correlations: ", np.min(correls), np.max(correls)

    mfs_matrix_normalized = mfs_matrix.copy()
    measures_matrix_normalized = measures_matrix.copy()
    for i in range(mfs_matrix_normalized.shape[1]):
        mfs_matrix_normalized[:, i] = normalize(mfs_matrix_normalized[:, i])

    for i in range(measures_matrix_normalized.shape[1]):
        measures_matrix_normalized[:, i] = normalize(measures_matrix_normalized[:, i])

    correls_normalized = np.zeros((mfs_matrix_normalized.shape[1], measures_matrix_normalized.shape[1]))
    # compute correlations
    for d in range(mfs_matrix_normalized.shape[1]):
        for m in range(measures_matrix_normalized.shape[1]):
            corr = scipy.stats.stats.spearmanr(mfs_matrix_normalized[:, d],
                                               measures_matrix_normalized[:, m])[0]
            if not (math.isnan(corr)):
                correls_normalized[d, m] = corr
            else:
                correls_normalized[d, m] = 0

    print "Higher normalized correlations: ", np.min(correls_normalized), np.max(correls_normalized)

def compute_subset(measures_matrix, mfs,
                   mfs_pos_start_data, mfs_pos_end_data):
    mfs_subset = []
    for i in range(len(measures)):
        if not (math.isnan(measures_matrix[i][pos_fexp])):
            mfs_i = tuple(mfs[i])
            mfs_i = mfs_i[mfs_pos_start_data: mfs_pos_end_data + 1]
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

################################################################

measures_subset = np.array([])

# subset of measures
for i in range(len(measures)):
    if not(math.isnan(measures_matrix[i][pos_fexp])):
        if len(measures_subset) == 0:
            measures_subset = np.array(measures_matrix[i])
        else:
            measures_subset = np.vstack((measures_subset, measures_matrix[i]))


##############################################


mfs_pos_start_data = 1
mfs_pos_end_data = 20
print "Correlations with normalized MFS..."
compute_correlations(measures_matrix, mfs_normalized, mfs_pos_start_data,
                                        mfs_pos_end_data, True)


mfs_subset = compute_subset(measures_matrix, mfs_normalized, mfs_pos_start_data,
                            mfs_pos_end_data)


compute_linear_model(mfs_subset, measures_subset)
print ""
###########################################
#mfs_pos_start_data = 0
#mfs_pos_end_data = 21

#print "Correlations with Sandbox MFS Adaptive..."
#compute_correlations(measures_matrix, mfs_sandbox_adaptive, mfs_pos_start_data,
#                                        mfs_pos_end_data, True)


#mfs_sandbox_subset = compute_subset(measures_matrix, mfs_sandbox_adaptive,
#                                    mfs_pos_start_data, mfs_pos_end_data)

#compute_linear_model(mfs_sandbox_subset, measures_subset)
#print ""

###########################################
#mfs_pos_start_data = 0
#mfs_pos_end_data = 21

#print "Correlations with Sandbox_Absolute MFS_normalized ..."
#compute_correlations(measures_matrix, mfs_sandbox_absolute_normalized, mfs_pos_start_data,
#                                        mfs_pos_end_data, True)


#mfs_sandbox_absolute_subset = compute_subset(measures_matrix, mfs_sandbox_absolute_normalized,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_sandbox_absolute_subset, measures_subset)
#print ""


###############################################

mfs_pos_start_data = 0
mfs_pos_end_data = 20
print "Correlations with MFS..."
compute_correlations(measures_matrix, mfs, mfs_pos_start_data,
                                        mfs_pos_end_data)

# obtain subsets of 17 scans for Fexp

mfs_subset = compute_subset(measures_matrix, mfs,
                                    mfs_pos_start_data, mfs_pos_end_data)


compute_linear_model(mfs_subset, measures_subset, 'mfs_normalized')

print ""

###############################################
#mfs_pos_start_data = 0
#mfs_pos_end_data = 6*5
#print "Correlations with Pyramid MFS (local or global) ..."
#compute_correlations(measures_matrix, mfs_local_pyramid, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_local_pyramid,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""



###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 6
#print "Correlations with Local MFS..."
#compute_correlations(measures_matrix, mfs_local, mfs_pos_start_data,
#                                        mfs_pos_end_data)

# obtain subsets of 17 scans for Fexp

#mfs_subset = compute_subset(measures_matrix, mfs_local,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""

###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 6
#print "Correlations with Pyramid Gradient MFS..."
#compute_correlations(measures_matrix, mfs_gradient_pyramid, mfs_pos_start_data,
#                                        mfs_pos_end_data)

# obtain subsets of 17 scans for Fexp

#mfs_subset = compute_subset(measures_matrix, mfs_gradient_pyramid,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""

###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 100
#print "Correlations with Slices X 2.5D MFS..."
#compute_correlations(measures_matrix, mfs_slices_x, mfs_pos_start_data,
#                                        mfs_pos_end_data)

# obtain subsets of 17 scans for Fexp

#mfs_subset = compute_subset(measures_matrix, mfs_slices_x,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""

###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 100
#print "Correlations with Slices Y 2.5D MFS..."
#compute_correlations(measures_matrix, mfs_slices_y, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_slices_y,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""


###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 100
#print "Correlations with Slices Z 2.5D MFS..."
#compute_correlations(measures_matrix, mfs_slices_z, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_slices_z,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""

###############################################

mfs_pos_start_data = 0
mfs_pos_end_data = 100
print "Correlations with Pure Pyramid MFS..."
compute_correlations(measures_matrix, mfs_pure_pyramid, mfs_pos_start_data,
                                        mfs_pos_end_data)

mfs_subset = compute_subset(measures_matrix, mfs_pure_pyramid,
                                    mfs_pos_start_data, mfs_pos_end_data)

np.savetxt('pyramid.csv', mfs_subset, delimiter=",")


compute_linear_model(mfs_subset, measures_subset)

print ""

###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 100
#print "Correlations with Pure Pyramid Gradient MFS..."
#compute_correlations(measures_matrix, mfs_pure_pyramid_gradient, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_pure_pyramid_gradient,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""

###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 20
#print "Correlations with Sigmoid MFS..."
#compute_correlations(measures_matrix, mfs_sigmoid, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_sigmoid,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""


###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 10
#print "Correlations with MFS Holder (10 Dimentions)..."
#compute_correlations(measures_matrix, mfs_holder_10, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_holder_10,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""


###############################################

#mfs_pos_start_data = 0
#mfs_pos_end_data = 5
#print "Correlations with MFS Holder (5 Dimentions)..."
#compute_correlations(measures_matrix, mfs_holder_5, mfs_pos_start_data,
#                                        mfs_pos_end_data)

#mfs_subset = compute_subset(measures_matrix, mfs_holder_5,
#                                    mfs_pos_start_data, mfs_pos_end_data)


#compute_linear_model(mfs_subset, measures_subset)

#print ""

###############################################

mfs_pos_start_data = 0
mfs_pos_end_data = 5
print "Correlations with Stats MFS Holder..."
compute_correlations(measures_matrix, stats_mfs_holder, mfs_pos_start_data,
                                        mfs_pos_end_data)

mfs_subset = compute_subset(measures_matrix, stats_mfs_holder,
                                    mfs_pos_start_data, mfs_pos_end_data)


compute_linear_model(mfs_subset, measures_subset)

print ""