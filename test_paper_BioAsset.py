import numpy as np
from imfractal import *
from numpy import recfromcsv
import scipy
import math
import pandas
from pandas import DataFrame, Series
import statsmodels.formula.api as sm
import statsmodels
from statsmodels.tools.eval_measures import aicc, bic, hqic, rmse

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
stats_mfs_pyramid = np.load('exps/data/mfs_stats_pyramid.npy')
stats_mfs_pyramid_gradient = np.load('exps/data/mfs_stats_pyramid_gradient.npy')
stats_mfs_slices_x = np.load('exps/data/stats_mfs_slices_x.npy')
stats_mfs_slices_y = np.load('exps/data/stats_mfs_slices_y.npy')
stats_mfs_slices_z = np.load('exps/data/stats_mfs_slices_z.npy')


pos_fexp = 17 #check

def normalize(vector):
    if np.std(vector) != 0:
        return (vector - np.mean(vector))/ np.std(vector)
    else:
        print "std equals 0"
        return vector

def compute_robust_r2(y, X, total_model):
    # leave-one-out cross validation
    # (robust R^2 = cross-validated R^2)

    # Predictive Error Sum of Squares (PRESS)
    # Total Sum of Squares (TSS)
    # http://www.moleculardescriptors.eu/tutorials/T5_moleculardescriptors_models.pdf

    PRESS = 0.0
    TSS = 0.0

    # average experimental response
    y_line = np.mean(y)

    for i in range(X.shape[0]):
        arr = np.arange(X.shape[0])
        mask = np.ones(arr.shape, dtype=bool)
        mask[i] = 0

        model = sm.OLS(y[mask], X[mask])
        model = model.fit()

        p = model.predict(X[i])
        t_p = total_model.predict(X[i])

        PRESS += (y[i] - p) * (y[i] - p)
        TSS   += (y[i] - y_line) * (y[i] - y_line)

    return 1.0 - PRESS / TSS

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
    best_rmse = 100000
    best_rmse2 = 100000
    best_rmse3 = 100000
    best_rob_r2 = 0.0
    best_rob_r2_ij = 0.0
    best_rob_r2_ijk = 0.0


    for i in range(2, X2.shape[1]):
        Xi = X2[:, [0, 1, i]]

        model = sm.OLS(fexp, Xi)
        res = model.fit()

        aic = aicc(res.llf, res.nobs, res.params.shape[0])

        if aic < best_aicc :
            best_aicc = aic
            best_i = i-2
            best_r2 = res.rsquared_adj
            best_rmse = np.sqrt(res.mse_resid)
            best_rob_r2 = compute_robust_r2(fexp, Xi, res)

    for i in range(2, X2.shape[1]):
        for j in range(i+1, X2.shape[1]):
            Xij = X2[:, [0, 1, i, j]]


            model = sm.OLS(fexp, Xij)
            res = model.fit()

            aic = aicc(res.llf, res.nobs, res.params.shape[0])

            if aic < best_aicc2:
                best_aicc2 = aic
                best_i_j = [i-2, j-2]
                best_r2_ij = res.rsquared_adj
                best_rmse2 = np.sqrt(res.mse_resid)
                best_rob_r2_ij = compute_robust_r2(fexp, Xij, res)


    for i in range(2, X2.shape[1]):
        for j in range(i+1, X2.shape[1]):
            for k in range(j + 1, X2.shape[1]):
                Xijk = X2[:, [0, 1, i, j, k]]


                model = sm.OLS(fexp, Xijk)
                res = model.fit()

                aic = aicc(res.llf, res.nobs, res.params.shape[0])

                if aic < best_aicc3:
                    best_aicc3 = aic
                    best_i_j_k = [i-2, j-2, k-2]
                    best_r2_ijk = res.rsquared_adj
                    best_rmse3 = np.sqrt(res.mse_resid)
                    best_rob_r2_ijk = compute_robust_r2(fexp, Xijk, res)

    return best_aicc, best_i, best_r2, best_rob_r2, best_rmse,\
           best_aicc2, best_i_j, best_r2_ij, best_rob_r2_ij, best_rmse2,\
           best_aicc3, best_i_j_k, best_r2_ijk, best_rob_r2_ijk, best_rmse3

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


    import statsmodels.robust.robust_linear_model
    Xbmd = X[:, [0]]
    X2 = statsmodels.tools.tools.add_constant(Xbmd)

    #huber_t = sm.RLM(fexp, X2, M=statsmodels.robust.norms.HuberT())
    #m = huber_t.fit()
    #print m.rsquared
    #print m.summary()
    #exit()


    model = sm.OLS(fexp, X2)
    res = model.fit()

    model_bmd = sm.OLS(fexp, Xbmd)
    model_bmd = model_bmd.fit()

    aic = aicc(res.llf, res.nobs, res.params.shape[0])
    r2 = res.rsquared_adj
    rmsee = np.sqrt(res.mse_resid)
    rob_r2 = compute_robust_r2(fexp, Xbmd, model_bmd)
    #print "BMD AICc, dimension, R2: " , aic, ' bmd ', r2

    res = compute_best_aicc(X, fexp)
    #print "AICc, dimension, R2: ", res[0], ' bmd + ', res[1], res[2]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[0])/2.0)
    #print "AICc, dimensions, R2: ", res[3],' bmd + ',  res[4], res[5]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[3]) / 2.0)
    #print "AICc, dimensions, R2: ", res[6],' bmd + ',  res[7], res[8]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[6]) / 2.0)

    return aic, r2, rob_r2, rmsee, res

    #X_normalized = X
    #for i in range(X.shape[1]):
    #    X_normalized[:,i] = normalize(X_normalized[:,i])

    #res_n = compute_best_aicc(X_normalized, fexp)
    #print "AICc, dimension, R2: ", res_n[0], ' bmd + ', res_n[1], res_n[2]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res_n[0]) / 2.0)
    #print "AICc, dimensions, R2: ", res_n[3], ' bmd + ', res_n[4], res_n[5]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res_n[3]) / 2.0)
    #print "AICc, dimensions, R2: ", res_n[6], ' bmd + ', res_n[7], res_n[8]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res_n[6]) / 2.0)

    #print ""

    #print "Normalized Variables - Score (R^2):", clf.score(X_normalized, fexp)

    #X_2 = statsmodels.tools.tools.add_constant(X_normalized)
    #model_2 = sm.OLS(fexp, X_2)
    #res_2 = model_2.fit()

    np.savetxt(data_path + output_file, X_normalized, delimiter=",")


def compute_correlations(measures_matrix, mfs, mfs_pos_start_data,
                                        mfs_pos_end_data, transform_mfs = True):
    print "///////////////////////"

    result = []

    # convert from ugly format to matrix
    mfs_matrix = []

    #print mfs.shape

    #if transform_mfs:
    if len(mfs.shape) == 1:
        for i in range(mfs.shape[0]):
            mfs_i = tuple(mfs[i])
            mfs_i = mfs_i[mfs_pos_start_data : mfs_pos_end_data + 1]


            if len(mfs_matrix) == 0:
               mfs_matrix = mfs_i
            else:
               mfs_matrix = np.vstack((mfs_matrix, mfs_i))

    else: mfs_matrix = mfs


    #print mfs_matrix.shape
    #print measures_matrix.shape

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
    #print "Higher correlations: ", np.min(correls), np.max(correls)

    return np.min(correls), np.max(correls)

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

str_method = [
    "Normalized MFS",
    "Stats Pyramid MFS",
    "Stats Pyramid Gradient MFS",
    "Stats MFS Slices X",
    "Stats MFS Slices Y",
    "Stats MFS Slices Z",
    "MFS",
    #"Pyramid MFS (local or global) ...",
    #"Sandbox MFS Adaptive",
    #"Sandbox_Absolute MFS_normalized ...",
    #"Local MFS",
    #"Sigmoid MFS",
    #"10-MFS",
    #"5-MFS",
    "Stats MFS Holder",
    "Pyramid MFS",
    #"Pyramid Gradient MFS",
    "Slices X MFS",
    "Slices Y MFS",
    "Slices Z MFS",
    "Pure Pyramid Gradient MFS"
]
method_array = [
    mfs_normalized,
    stats_mfs_pyramid,
    stats_mfs_pyramid_gradient,
    stats_mfs_slices_x,
    stats_mfs_slices_y,
    stats_mfs_slices_z,
    mfs,
    #mfs_local_pyramid,
    #mfs_sandbox_adaptive,
    #mfs_sandbox_absolute_normalized,
    #mfs_local,
    #mfs_sigmoid,
    #mfs_holder_10,
    #mfs_holder_5,
    stats_mfs_holder,
    mfs_pure_pyramid,
    #mfs_gradient_pyramid,
    mfs_slices_x,
    mfs_slices_y,
    mfs_slices_z,
    mfs_pure_pyramid_gradient
]


num_array_begin = [
    1,
    0,
    0,
    0,
    0,
    0,
    0,
    #0,
    #0,
    #0,
    #0,
    #0,
    #0,
    #0,
    0,
    0,
    0,
    0,
    0,
    0,
    #0
]

num_array_end = [
    20,
    50,
    50,
    50,
    50,
    50,
    20,
    #30,
    #21,
    #21,
    #6,
    #20,
    #10,
    #5,
    10,
    100,
    100,
    100,
    100,
    100,
    #100
]

for i in range(len(method_array)):

    c1, c2 = compute_correlations(measures_matrix, method_array[i], num_array_begin[i],
                                            num_array_end[i], True)

    print str_method[i],  " #", num_array_end[i]

    if(len(method_array[i]) > 17):
        mfs_subset = compute_subset(measures_matrix, method_array[i], num_array_begin[i],
                                num_array_end[i])
    else:
        mfs_subset = method_array[i]

    #np.savetxt('pyramid.csv', mfs_subset, delimiter=",")

    aic, r2, rob_r2, rmsee, res = compute_linear_model(mfs_subset, measures_subset)

    aic_s = "%.4f" % aic
    rmsee = "%.4f" % rmsee
    r2 = "%.4f" % r2
    rob_r2 = "%.4f" % rob_r2

    aicc1 = "%.4f" % res[0]
    aicc2 = "%.4f" % res[5]
    aicc3 = "%.4f" % res[10]
    r2_1 = "%.4f" % res[2]
    r2_2 = "%.4f" % res[7]
    r2_3 = "%.4f" % res[12]
    dims_1 = res[1]
    dims_2 = res[6]
    dims_3 = res[11]
    p1 = "%.7f" % np.exp((min(aic, res[0])-max(aic, res[0])) / 2.0)
    p2 = "%.7f" % np.exp((min(aic, res[5])-max(aic, res[5])) / 2.0)
    p3 = "%.7f" % np.exp((min(aic, res[10])-max(aic, res[10])) / 2.0)
    rmse1 = "%.5f" % res[4]
    rmse2 = "%.5f" % res[9]
    rmse3 = "%.5f" % res[14]

    rob_r2_1 = "%.4f" % res[3]
    rob_r2_2 = "%.4f" % res[8]
    rob_r2_3 = "%.4f" % res[13]



    print "AICC     -  Adj R^2 -   Rob. R^2  -   DIMS            - p-value   -   RMSE"
    print aic_s, ' |',  r2,   '  |',  rob_r2,   '     | bmd   ', "           | --------- | ", rmsee
    print aicc1, ' |',  r2_1, '  |',  rob_r2_1, '     | bmd + ', dims_1, "         |", p1, "| ", rmse1
    print aicc2, ' |',  r2_2, '  |',  rob_r2_2, '     | bmd + ', dims_2, "    |", p2, "| ", rmse2
    print aicc3, ' |',  r2_3, '  |',  rob_r2_3, '     | bmd + ', dims_3, "|", p3, "| ", rmse3



