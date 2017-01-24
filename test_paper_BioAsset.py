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
import matplotlib.pyplot as plt

from matplotlib.mlab import griddata
import itertools
import numpy as np
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D
from pylab import *


np.seterr(divide='ignore', invalid='ignore')

def ribbon1(data):
    x=data[:,0]
    fig=plt.figure()
    ax=fig.gca(projection='3d')

    for i in range(1,data.shape[1]-1,1):
        y=data[:,i]
        z=data[:,i+1]
        xi=np.linspace(min(x),max(x))
        yi=np.linspace(min(y),max(y))
        X,Y=np.meshgrid(xi,yi)
        Z=griddata(x,y,z,xi,yi)
        ax.plot_surface(X,Y,Z,rstride=50,cstride=1,cmap='RdYlBu')
        ax.set_zlim3d(np.min(Z),np.max(Z))

    ax.set_zlim3d(0.0,3.0)
    ax.set_title('Gradient')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_yticks([])
    ax.set_zlabel('')
    plt.show()

def ribbon2(spectra, alpha, title):
    import plotly.plotly as py
    import plotly.graph_objs as go

    import numpy as np

    traces = []
    y_raw = spectra[:, 0]
    sample_size = spectra.shape[1]
    for i in range(0, sample_size):
        print i
        z_raw = spectra[:, i]
        x = []
        y = []
        z = []
        ci = int(255/sample_size*i) # ci = "color index"
        for j in range(0, len(z_raw)):
            z.append([z_raw[j], z_raw[j]])
            y.append([alpha[j], alpha[j]])
            x.append([i, i+1])
        traces.append(dict(
            xsrc="xsrc",
            z=z,
            x=x,
            y=y,
            colorscale=[ [i, 'rgb(%d,%d,255)'%(ci, ci)] for i in np.arange(0,1.1,0.1) ],
            showscale=False,
            type='surface',
        ))

    #traces = traces[::-1]

    fig = go.Figure( data=traces, layout={'title':'Pyramid Gradient MFS'})
    py.iplot(fig, filename='ribbon-plot-python')


def ribbon3(spectra, alpha):
    matplotlib.rcParams.update({'font.size':10})
    fig=figure()
    ax=fig.gca(projection='3d')
    for i in range(0,spectra.shape[1]):
        y=spectra[:,i]
        x=sorted(range(1,len(y)+1)*2)
        a=[i,i+1]*len(y)
        b=list(itertools.chain(*zip(y,y)))
        xi=np.linspace(min(x),max(x))
        yi=np.linspace(min(a),max(a))
        X,Y=np.meshgrid(xi/(len(x)*0.5),yi)
        Z=griddata(x,a,b,xi,yi)
        ax.plot_surface(X,Y,Z, rstride=50, cstride=1, cmap='Spectral')
        ax.set_zlim3d(np.min(Z),np.max(Z))

    ax.grid(False)
    ax.w_xaxis.pane.set_visible(False)
    ax.w_yaxis.pane.set_visible(False)
    ax.w_zaxis.pane.set_color('gainsboro')
    ax.set_title('Pyramid Gradient MFS')
    ax.set_xlim3d(0,1)
    ax.set_xticks(alpha)
    ax.set_xticklabels(alpha)
    ax.set_xlabel(r'$\alpha$')
    #ax.set_yticks([0.5,1.5,2.5,3.5,4.5])
    #ax.set_yticklabels(['1','2','3','4','5'])
    ax.set_ylabel('Resolution')
    ax.set_zlim3d(0,3)
    ax.set_zlabel(r'$f(\alpha)$')
    show()

def csvToNumpy(X):

    x = tuple(X[0])
    x = np.array(x[1:])

    X_res = x

    for i in range(1, len(X)):
        x = tuple(X[i])
        x = np.array(x[1:])

        X_res = np.vstack((X_res, x))

    return X_res

# for Fexp
#fexp_names = np.load(data_path + 'bioAsset_meta.npy')

# Adaptive Metadata and mfs
measures = recfromcsv('exps/data/BioAssetAdaptiveThresh/default_BioAsset_Adaptive.csv', delimiter=',')
measures_npy = csvToNumpy(measures)
mfs = np.load('exps/data/mfs_holder_BioAsset_raw.npy')
mfs_normalized = recfromcsv('exps/data/BioAssetAdaptiveThresh/mfs_holder_BioAsset.csv', delimiter=',')
mfs_normalized = csvToNumpy(mfs_normalized)
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

#mfs_gradient = recfromcsv('exps/data/mfs_holder_gradient_BioAsset.csv', delimiter=',')
#mfs_gradient = csvToNumpy(mfs_gradient)
mfs_gradient = np.load('exps/data/mfs_gradient.npy')

mfs_laplacian = np.load('exps/data/mfs_laplacian.npy')

stats_mfs_holder = np.load('exps/data/stats_mfs_holder_BioAsset.npy')
stats_mfs_holder2 = np.load('exps/data/mfs_stats_2.npy')


stats_mfs_pyramid = np.load('exps/data/mfs_stats_pyramid.npy')
stats_mfs_pyramid_gradient = np.load('exps/data/mfs_stats_pyramid_gradient.npy')
stats_mfs_slices_x = np.load('exps/data/stats_mfs_slices_x.npy')
stats_mfs_slices_y = np.load('exps/data/stats_mfs_slices_y.npy')
stats_mfs_slices_z = np.load('exps/data/stats_mfs_slices_z.npy')

#stats_mfs_pyramid_gradient = np.hstack((stats_mfs_pyramid_gradient[:, 5:45])) #, #stats_mfs_pyramid_gradient[:, 36:47]))
#stats_mfs_pyramid_gradient = stats_mfs_pyramid_gradient[:, [1, 17,46]] #, #stats_mfs_pyramid_gradient[:, 36:47]))
#print stats_mfs_pyramid_gradient
#mfs = np.hstack([mfs[:, 6:10], mfs[:, 17:20]])

#mfs_slices_x[i] = np.hstack([mfs_slices_x[:, 6:10], mfs_slices_x[:, 17:20]])

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

    sum_rmse = 0.0

    for i in range(X.shape[0]):
        arr = np.arange(X.shape[0])
        mask = np.ones(arr.shape, dtype=bool)
        mask[i] = 0

        model = sm.OLS(y[mask], X[mask])
        model = model.fit()

        p = model.predict(X[i])

        PRESS += (y[i] - p) * (y[i] - p)
        TSS   += (y[i] - y_line) * (y[i] - y_line)

        sum_rmse += model.mse_resid # RMSE^2

    rob_r2 = 1.0 - PRESS / TSS
    rob_rmse = np.sqrt(sum_rmse / X.shape[0])

    return rob_r2, rob_rmse

def compute_best_aicc(X, fexp):

    # X: [BMD MFS]
    # X2 : [CTE BMD MFS]

    #X2 = statsmodels.tools.tools.add_constant(X)
    #X2 = np.hstack((np.ones(X.shape[0]), X))
    X2 = np.append(np.ones((X.shape[0], 1)), X, axis = 1)

    #print np.ones((X.shape[0], 1)).shape
    #print X.shape

    if X2.shape == X.shape:
        print "Error in add_constant!"
        exit()

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

    best_rob_rmse = 100000
    best_rob_rmse2 = 100000
    best_rob_rmse3 = 100000

    best_rob_r2 = 0.0
    best_rob_r2_ij = 0.0
    best_rob_r2_ijk = 0.0


    for i in range(2, X2.shape[1]):

        # 0: constant, 1: BMD

        Xi = X2[:, [0, 1, i]]

        first_c = np.any(np.array(X2[:,0]).astype(np.int32) ==
                          np.ones(X2.shape[0]).astype(np.int32))
        if not(np.any(first_c)):
            print "NOT!!"
            continue
        #print ""
        #print i

        model = sm.OLS(fexp, Xi)
        res = model.fit()

        rob_r2, rob_rmse = compute_robust_r2(fexp, Xi, res)

        #if aic < best_aicc :
        if rob_r2 > best_rob_r2:
            #best_aicc = aic
            best_i = i-2
            best_r2 = res.rsquared_adj
            best_rmse = np.sqrt(res.mse_resid)
            best_rob_r2 = rob_r2
            best_rob_rmse = rob_rmse
            best_aicc = aicc(res.llf, res.nobs, res.params.shape[0])
            #best_rob_r2, best_rob_rmse = compute_robust_r2(fexp, Xi, res)

    for i in range(2, X2.shape[1]):
        for j in range(i+1, X2.shape[1]):
            Xij = X2[:, [0, 1, i, j]]

            first_c = np.any(np.array(X2[:, 0]).astype(np.int32) ==
                             np.ones(X2.shape[0]).astype(np.int32))
            if not (np.any(first_c)):
                print "NOT!!"
                continue

            model = sm.OLS(fexp, Xij)
            res = model.fit()

            rob_r2_ij, rob_rmse2 = compute_robust_r2(fexp, Xij, res)

            #if aic2 < best_aicc2:
            if rob_r2_ij > best_rob_r2_ij:
                #best_aicc2 = aic2
                best_i_j = [i-2, j-2]
                best_r2_ij = res.rsquared_adj
                best_rmse2 = np.sqrt(res.mse_resid)
                best_rob_r2_ij = rob_r2_ij
                best_rob_rmse2 = rob_rmse2
                best_aicc2 = aicc(res.llf, res.nobs, res.params.shape[0])

                #best_rob_r2_ij, best_rob_rmse2 = compute_robust_r2(fexp, Xij, res)


    for i in range(2, X2.shape[1]):
        for j in range(i+1, X2.shape[1]):
            for k in range(j + 1, X2.shape[1]):
                Xijk = X2[:, [0, 1, i, j, k]]

                first_c = np.any(np.array(X2[:, 0]).astype(np.int32) ==
                                 np.ones(X2.shape[0]).astype(np.int32))
                if not (np.any(first_c)):
                    print "NOT!!"
                    continue


                model = sm.OLS(fexp, Xijk)
                res = model.fit()

                rob_r2_ijk, rob_rmse3 = compute_robust_r2(fexp, Xijk, res)

                if rob_r2_ijk > best_rob_r2_ijk:
                #if aic3 < best_aicc3:
                    #best_aicc3 = aic3
                    best_i_j_k = [i-2, j-2, k-2]
                    best_r2_ijk = res.rsquared_adj
                    best_rmse3 = np.sqrt(res.mse_resid)
                    #best_rob_r2_ijk, best_rob_rmse3 = compute_robust_r2(fexp, Xijk, res)
                    best_rob_r2_ijk = rob_r2_ijk
                    best_rob_rmse3 = rob_rmse3
                    best_aicc3 = aicc(res.llf, res.nobs, res.params.shape[0])

    return best_aicc, best_i, best_r2, best_rob_r2, best_rmse, best_rob_rmse,\
           best_aicc2, best_i_j, best_r2_ij, best_rob_r2_ij, best_rmse2, best_rob_rmse2,\
           best_aicc3, best_i_j_k, best_r2_ijk, best_rob_r2_ijk, best_rmse3, best_rob_rmse3

def compute_linear_model(mfs, measures, output_file="standarized.csv"):
    #from sklearn.linear_model import Ridge
    from sklearn import linear_model

    # try different ones
    #clf = Ridge(alpha = 1.0)
    #clf = RidgeCV(alphas=[0.1, 1.0, 10.0])
    #clf = linear_model.LinearRegression()

    # explain fexp using BMD + the MFS data

    bmd = measures[:, 0]
    fexp = measures[:, measures.shape[1]-1]


    #print "BMD: ", bmd
    #print "FEXP: ", fexp
    #print "MFS; ", mfs

    #PCA
    #from sklearn.decomposition import PCA
    #pca = PCA(n_components=8)
    #pca.fit(mfs)
    #mfs = pca.transform(mfs)

    X = np.hstack((bmd, mfs))
    #clf.fit(X, fexp)
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


    #import statsmodels.robust.robust_linear_model
    #Xbmd = X[:, [0]]
    X2 = statsmodels.tools.tools.add_constant(bmd)

    #huber_t = sm.RLM(fexp, X2, M=statsmodels.robust.norms.HuberT())
    #m = huber_t.fit()
    #print m.rsquared
    #print m.summary()
    #exit()


    model = sm.OLS(fexp, X2)
    res = model.fit()

    aic = aicc(res.llf, res.nobs, res.params.shape[0])
    r2 = res.rsquared_adj
    rmsee = np.sqrt(res.mse_resid)
    rob_r2, rob_rmse = compute_robust_r2(fexp, X2, res)
    #print "BMD AICc, dimension, R2: " , aic, ' bmd ', r2

    res = compute_best_aicc(X, fexp)
    #print "AICc, dimension, R2: ", res[0], ' bmd + ', res[1], res[2]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[0])/2.0)
    #print "AICc, dimensions, R2: ", res[3],' bmd + ',  res[4], res[5]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[3]) / 2.0)
    #print "AICc, dimensions, R2: ", res[6],' bmd + ',  res[7], res[8]
    #print "AICc p-value (significance): ", 1.0 / np.exp((aic - res[6]) / 2.0)

    return aic, r2, rob_r2, rmsee, rob_rmse, res

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


def compute_correlations(X, Y):#, mfs_pos_start_data,
                                        #mfs_pos_end_data, transform_mfs = True):
    #result = []

    # convert from ugly format to matrix
    #mfs_matrix = mfs

    #print mfs.shape

    #if transform_mfs:
    #if len(mfs.shape) == 1:
    #    for i in range(mfs.shape[0]):
    #        mfs_i = tuple(mfs[i])
    #        mfs_i = mfs_i[mfs_pos_start_data : mfs_pos_end_data + 1]


    #        if len(mfs_matrix) == 0:
    #           mfs_matrix = mfs_i
    #        else:
    #           mfs_matrix = np.vstack((mfs_matrix, mfs_i))

    #else: mfs_matrix = mfs


    #print mfs_matrix.shape
    #print measures_matrix.shape

    correls = np.zeros((X.shape[1], Y.shape[1]))
    # compute correlations
    for x in range(X.shape[1]):
        for y in range(Y.shape[1]):
            corr = scipy.stats.stats.spearmanr(X[:, x],
                                               Y[:, y])[0]

            if not(math.isnan(corr)):
                correls[x, y] = corr
            else:
                #print "Nan in correlations", x, y
                correls[x, y] = 0

    #DEBUG CORRELATIONS:
    #import matplotlib.pyplot as plt
    #plt.plot(correls)
    #plt.show()
    #print "Higher correlations: ", np.min(correls), np.max(correls)

    return np.min(correls), np.max(correls), correls

    #mfs_matrix_normalized = mfs_matrix.copy()
    #measures_matrix_normalized = measures_matrix.copy()
    #for i in range(mfs_matrix_normalized.shape[1]):
    #    mfs_matrix_normalized[:, i] = normalize(mfs_matrix_normalized[:, i])

    #for i in range(measures_matrix_normalized.shape[1]):
    #    measures_matrix_normalized[:, i] = normalize(measures_matrix_normalized[:, i])

    #correls_normalized = np.zeros((mfs_matrix_normalized.shape[1], measures_matrix_normalized.shape[1]))
    # compute correlations
    #for d in range(mfs_matrix_normalized.shape[1]):
    #    for m in range(measures_matrix_normalized.shape[1]):
    #        corr = scipy.stats.stats.spearmanr(mfs_matrix_normalized[:, d],
    #                                           measures_matrix_normalized[:, m])[0]
    #        if not (math.isnan(corr)):
    #            correls_normalized[d, m] = corr
    #        else:
    #            correls_normalized[d, m] = 0

    #print "Higher normalized correlations: ", np.min(correls_normalized), np.max(correls_normalized)

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

#################################################################
# Paper figures

# Figure : MFS Sandbox of volume 1_2
fsize = 15
mfs1_2 = mfs_sandbox_absolute_normalized[1] #mfs_sandbox_adaptive[1]
x1 = -10
x2 = 10
x = np.arange(x1, x2+1)
plt.ylim((2.2, 3.5))
plt.xlim(x1,x2)
plt.ylabel('$D_{q}$',fontsize=fsize)
plt.xlabel('Generalised dimension',fontsize=fsize)

print len(x), len(mfs1_2)

plt.plot(x, mfs1_2, '*-', linewidth=2.0)
plt.show()

# Figure: MFS Holder of volume 1_2

xt = np.arange(1,20,2)
alpha = [   6.,   19.,   32.,   45.,   58.,   71.,   84.,   97.,  110.,  123. , 136. , 149.,
162.,  175.,  188.,  201. , 214. , 227. , 240.  ,253.]
alpha = map(lambda i:"%.2f" % float(i/255.), alpha)

alpha_orig = alpha
alpha = alpha[0:len(alpha):2]

plt.xticks(xt,alpha) # translate


fsize = 15
mfs1_2 = mfs[1]
#x1 = -10
#x2 = 10
#x = np.arange(x1, x2+1)
plt.ylim((0.0, 3.0))
#plt.xlim(x1,x2)
plt.ylabel(r'$f(\alpha)$', fontsize=fsize)
plt.xlabel(r'$\alpha$', fontsize=fsize)

x = np.arange(len(mfs1_2))
plt.xticks(xt,alpha) # translate
plt.plot(map(lambda i: i+1, x), mfs1_2, '*-', linewidth=2.0)
plt.show()


# Figure: MFS Pyramid of volume 1_2
fsize = 15
mfs1_2 = mfs_pure_pyramid[1]

# reverse (bug)

a = mfs1_2[:20]
b = mfs1_2[20:40]
c = mfs1_2[40:60]
d = mfs1_2[60:80]
e = mfs1_2[80:]

plt.xticks(xt,alpha) # translate
mfs1_2 = np.hstack((e,d,c,b,a))

#x1 = -10
#x2 = 10
#x = np.arange(x1, x2+1)
plt.ylim((0.0, 3.0))
#plt.xlim(x1,x2)
plt.ylabel(r'$f(\alpha)$', fontsize=fsize)
#plt.xlabel(r'$\alpha$', fontsize=fsize)

#x = np.arange(len(mfs1_2))
#plt.plot(mfs1_2, '*-', linewidth=2.0)
plt.plot(e, '*-', linewidth=2.0, label='0')
plt.plot(d, '*-', linewidth=2.0, label='1')
plt.plot(c, '*-', linewidth=2.0, label='2')
plt.plot(b, '*-', linewidth=2.0, label='3')
plt.plot(a, '*-', linewidth=2.0, label='4')
plt.xlabel(r'$\alpha$', fontsize=fsize)
plt.legend()
plt.show()

# Figure: Gradient MFS Pyramid of volume 1_2
fsize = 15
mfs1_2 = mfs_pure_pyramid_gradient[1]
plt.xticks(xt,alpha) # translate

a = mfs1_2[:20]
b = mfs1_2[20:40]
c = mfs1_2[40:60]
d = mfs1_2[60:80]
e = mfs1_2[80:]

mfs1_2 = np.hstack((e,d,c,b,a))
plt.plot(e, '*-', linewidth=2.0, label='0')
plt.plot(d, '*-', linewidth=2.0, label='1')
plt.plot(c, '*-', linewidth=2.0, label='2')
plt.plot(b, '*-', linewidth=2.0, label='3')
plt.plot(a, '*-', linewidth=2.0, label='4')
plt.xlabel(r'$\alpha$', fontsize=fsize)
#x1 = -10
#x2 = 10
#x = np.arange(x1, x2+1)
plt.ylim((0.0, 3.0))
#plt.xlim(x1,x2)
plt.ylabel(r'$f(\alpha)$', fontsize=fsize)
plt.legend()
#plt.xlabel(r'$\alpha$', fontsize=fsize)

#x = np.arange(len(mfs1_2))
#plt.plot(mfs1_2, '*-', linewidth=2.0)
plt.show()

#ribbon = np.vstack((a, b, c, d, e))
#print "Shape ribbon:" , ribbon.shape
#ribbon2(ribbon.T, alpha_orig, 'Pyramid Gradient MFS')
#exit()

####### 3D SCATTER BMD vs FEXP vs SK_{0}
from mpl_toolkits.mplot3d import Axes3D

indexes = np.load('exps/data/valid_fexp_indexes.npy')
bmd17 = np.load('exps/data/bmd17.npy')
fexp = np.load('exps/data/fexp.npy')
sk0 = np.load('exps/data/sk0.npy')
sk0_17 = sk0[indexes]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = ax.scatter(sk0_17, bmd17, fexp, s=(4*fexp**2)**2, c="r")
ax.set_xlabel(r'$SK_{0}$')
ax.set_ylabel('BMD')
ax.set_zlabel(r'$F_{exp}$')

plt.show()

colors=['b', 'c', 'y', 'm', 'r']

ax = plt.subplot(111, projection='3d')

ax.plot(bmd17, sk0_17, fexp, color=colors[0], label='Fexp')
#ax.plot(random(10), random(10), random(10), 'o', color=colors[0], label='LoLo')
#ax.plot(random(10), random(10), random(10), 'o', color=colors[1], label='Lo')
#ax.plot(random(10), random(10), random(10), 'o', color=colors[2], label='Average')
#ax.plot(random(10), random(10), random(10), 'o', color=colors[3], label='Hi')
#ax.plot(random(10), random(10), random(10), 'o', color=colors[4], label='HiHi')
#ax.plot(random(10), random(10), random(10), 'x', color=colors[4], label='High Outlier')

plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8, bbox_to_anchor=(0, 0))

plt.show()

# 2D Scatters
plt.xlabel('BMD')
plt.ylabel(r'F$^{exp}$')
plt.scatter(bmd17, fexp)
plt.show()

#
plt.xlabel(r'SK$_{0}$')
plt.ylabel(r'F$^{exp}$')
plt.scatter(sk0_17, fexp)
plt.show()

#
plt.xlabel(r'SK$_{0}$')
plt.ylabel('BMD')
plt.scatter(sk0_17, bmd17)
plt.show()





#################################################################

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
    #"MFS",
    #"Standard Measures",
    #"Normalized MFS",
    "Gradient MFS",
    "Stats Pyramid Gradient MFS",
    #"Normalized MFS",
    #"Stats Pyramid MFS",
    #"Pyramid MFS (local or global) ...",
    #"Sandbox MFS Adaptive",
    #"Sandbox_Absolute MFS_normalized ...",
    #"Local MFS",
    #"Sigmoid MFS",
    #"10-MFS",
    #"5-MFS",
    "Stats MFS Holder",
    #"Stats MFS Holder 2",
    #"Pyramid MFS",
    "Pyramid Gradient MFS",
    #"Slices X MFS",
    #"Slices Y MFS",
    #"Slices Z MFS",
    "Pure Pyramid Gradient MFS",
    #"Stats MFS Slices X",
    #"Stats MFS Slices Y",
    #"Stats MFS Slices Z",
]
method_array = [
    #mfs,
    #measures_npy[:, :-2],
    #mfs_normalized,
    mfs_gradient,
    stats_mfs_pyramid_gradient,
    #mfs_normalized,
    #stats_mfs_pyramid,
    #mfs_local_pyramid,
    #mfs_sandbox_adaptive,
    #mfs_sandbox_absolute_normalized,
    #mfs_local,
    #mfs_sigmoid,
    #mfs_holder_10,
    #mfs_holder_5,
    stats_mfs_holder,
    #stats_mfs_holder2,
    #mfs_pure_pyramid,
    mfs_gradient_pyramid,
    #mfs_slices_x,
    #mfs_slices_y,
    #mfs_slices_z,
    mfs_pure_pyramid_gradient,
    #stats_mfs_slices_x,
    #stats_mfs_slices_y,
    #stats_mfs_slices_z,
]

#print "Auto-correlations:"
#print "MFS: "
c1, c2, c = compute_correlations(mfs, mfs)

mask1 = []
sums = []
for i in range(len(c)):
    #print i, c[i]
    s = np.median(np.abs(c[i]))
    sums = np.append(sums, s)
    arr = np.arange(c.shape[1])
    mask = np.ones(arr.shape, dtype = bool)
    mask[i] = 0
    th = 0.89
    if not(np.any(c[i][mask] > th) or np.any(c[i][mask] < -th)):
        #print c[i][mask]
        #print i
        mask1 = np.append(mask1, i)

mask1 = np.append(mask1, 1)
mask1 = np.append(mask1, 17)

mask1 = np.sort(mask1)

#print sums

#plt.plot(sums.T)
#plt.show()
#exit()

mask1 = mask1.astype(np.int32)
#print mask1

sk0 = stats_mfs_pyramid_gradient[:, 36:37]
np.save('exps/data/sk0.npy', sk0)

exit()

if(False):

    print "Standard Measures intra-correlations:"
    c1, c2, c = compute_correlations(measures_matrix, measures_matrix)


    np.set_printoptions(suppress=True)
    print c

    print "Multifractal Skewness - Standard Measures : correlations:"
    skew_levels = stats_mfs_pyramid_gradient[:, [6,36,46]]
    c1, c2, c = compute_correlations(measures_matrix, skew_levels)

    np.set_printoptions(suppress=True)
    print c

    print "Multifractal Skewness intra-correlations:"
    skew_levels = stats_mfs_pyramid_gradient[:, [6,36,46]]
    c1, c2, c = compute_correlations(skew_levels, skew_levels)

    np.set_printoptions(suppress=True)
    print c

    fexp = np.array(measures_matrix[:, measures_matrix.shape[1] - 1]).reshape(measures_matrix.shape[0],1)
    rest = np.hstack((skew_levels, measures_matrix))

    rest_subset = compute_subset(measures_matrix, rest, 0, rest.shape[1])
    fexp_subset = compute_subset(measures_matrix, fexp, 0, fexp.shape[1])

    print "FEXP", fexp
    print fexp.shape

    indexes = []
    for i in range(len(fexp)):
        if not(np.isnan(fexp[i])): indexes.append(i)

    indexes = np.array(indexes).astype(np.uint32)

    print indexes

    print "Fexp against all correlations:"
    skew_levels = stats_mfs_pyramid_gradient[:, [6,36,46]]
    #c1, c2, c = compute_correlations(fexp_subset, rest_subset)

    np.set_printoptions(suppress=True)
#print c

#exit()
#print "Std Measures: ", compute_correlations(measures_matrix, measures_matrix)

for i in range(len(method_array)):

    #c1, c2, _ = compute_correlations(measures_matrix, method_array[i])

    #print str_method[i], " Min, Max correlations with #", measures_matrix.shape, "std measures: ", c1, c2



    if False:
        #mask1 = np.array([0,6,9,19]) #np.hstack((np.array(np.arange(1, 7)), np.array([16])))
        if method_array[i].shape[1] == 20:
            method_array[i] = method_array[i][:, mask1]

        if method_array[i].shape[1] == 100:

            #method_array[i] = method_array[i][:,60:80]

            if True:
                mask = []
                for j in range(5):
                    mask = np.hstack(( mask , (20*j)+  mask1)).astype(np.int32)

                #print mask
                #method_array[i] = np.hstack((method_array[i][:, 0:20], method_array[i][:, 80:100]))
                method_array[i] = method_array[i][:, mask]


    if(len(method_array[i]) > 17):
        mfs_subset = compute_subset(measures_matrix, method_array[i], 0, method_array[i].shape[1])
    else:
        mfs_subset = method_array[i]

    #plt.plot(mfs_subset.T)
    #plt.show()
    #exit()

    print str_method[i],  " #", method_array[i].shape[1]


    #np.savetxt('pyramid.csv', mfs_subset, delimiter=",")

    aic, r2, rob_r2, rmsee, rob_rmse, res = compute_linear_model(mfs_subset, measures_subset)

    aic_s = "%.4f" % aic
    rmsee = "%.4f" % rmsee
    r2 = "%.4f" % r2
    rob_r2 = "%.4f" % rob_r2
    rob_rmse = "%.4f" % rob_rmse

    aicc1 = "%.4f" % res[0]
    aicc2 = "%.4f" % res[6]
    aicc3 = "%.4f" % res[12]
    r2_1 = "%.4f" % res[2]
    r2_2 = "%.4f" % res[8]
    r2_3 = "%.4f" % res[14]
    dims_1 = res[1]
    dims_2 = res[7]
    dims_3 = res[13]
    p1 = "%.7f" % np.exp((min(aic, res[0])-max(aic, res[0])) / 2.0)
    p2 = "%.7f" % np.exp((min(aic, res[6])-max(aic, res[6])) / 2.0)
    p3 = "%.7f" % np.exp((min(aic, res[12])-max(aic, res[12])) / 2.0)
    rmse1 = "%.5f" % res[4]
    rmse2 = "%.5f" % res[10]
    rmse3 = "%.5f" % res[16]

    rob_rmse1 = "%.5f" % res[5]
    rob_rmse2 = "%.5f" % res[11]
    rob_rmse3 = "%.5f" % res[17]

    rob_r2_1 = "%.4f" % res[3]
    rob_r2_2 = "%.4f" % res[9]
    rob_r2_3 = "%.4f" % res[15]



    print "AICC     -  Adj R^2 - Rob. R^2 -   DIMS            - p-value   -   RMSE  - Rob. RMSE "
    print aic_s, ' |',  r2,   '  |',  rob_r2,   '  | bmd   ', "           | --------- | ", rmsee,  " | ", rob_rmse
    print aicc1, ' |',  r2_1, '  |',  rob_r2_1, '  | bmd + ', dims_1, "         |", p1, "| ", rmse1,  "| ", rob_rmse1
    print aicc2, ' |',  r2_2, '  |',  rob_r2_2, '  | bmd + ', dims_2, "    |", p2, "| ", rmse2,  "| ", rob_rmse2
    print aicc3, ' |',  r2_3, '  |',  rob_r2_3, '  | bmd + ', dims_3, "|", p3, "| ", rmse3,  "| ", rob_rmse3



