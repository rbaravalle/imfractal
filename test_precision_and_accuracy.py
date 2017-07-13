import numpy as np
from numpy import recfromcsv
import scipy.io as sio
from imfractal import *
from scipy import stats

def openMatlab(typ, filename, threshold = 100, adaptive = False):

    arr = np.array(sio.loadmat(filename)[typ])
    if typ == "False":
        if adaptive:
            threshold = self.determine_threshold(arr)

        arr = arr > threshold

        a_v = arr.cumsum()

        print "Amount of white pixels: ", a_v[len(a_v) - 1]

    # debug - to see the spongious structure
    # plt.imshow((arr[:,:,50]), cmap=plt.gray())
    # plt.show()

    return arr

#def prec_and_acc(data, f):
    #return f(data)

def get_ltp(data, xct):
    medians = np.zeros((data.shape[0], data.shape[1]))

    # median for every VOI
    for i in range(data.shape[0]):
        for ii in range(data.shape[1]):
            # tilde x_i
            medians[i,ii] = np.median(data[i,ii])

    print xct, xct.shape
    slope, intercept, r_value, p_value, std_err = stats.linregress(medians.flatten(), xct.flatten())
    print "SLOPE: " , slope

    ltp = 0
    for i in range(xct.shape[0]):
        for ii in range(xct.shape[1]):
            ltp += (xct[i,ii] - (slope * (medians[i,ii]) + intercept))**2

    print xct

    print ltp, ((xct.shape[0]*xct.shape[1] - 2) * (np.max(xct) - np.min(xct))**2)
    ltp = ltp/((xct.shape[0]*xct.shape[1] - 2) * (np.max(xct) - np.min(xct))**2)
    ltp = np.sqrt(ltp)

    return ltp

def get_stp(data):
    # data.shape = 4,5,10

    stp = 0

    means = np.zeros((data.shape[0], data.shape[1]))
    medians = np.zeros((data.shape[0], data.shape[1]))

    # mean and median for every VOI
    for i in range(data.shape[0]):
        for ii in range(data.shape[1]):
            means[i,ii] = np.mean(data[i,ii])
            medians[i,ii] = np.median(data[i,ii])

    res = 0
    for i in range(data.shape[0]):
        for ii in range(data.shape[1]):
            for j in range(data.shape[2]):
                res += (data[i,ii,j] - means[i,ii])**2

    res = res/(data.shape[0]*data.shape[1]*(data.shape[2]-1) * (np.max(medians) - np.min(medians))**2)
    res = np.sqrt(res)

    return res



def stp():
    sk0 = np.load('sk0.npy')
    sk1 = np.load('sk1.npy')
    sk2 = np.load('sk2.npy')
    sk3 = np.load('sk3.npy')
    sk4 = np.load('sk4.npy')
    kt0 = np.load('kt0.npy')
    kt1 = np.load('kt1.npy')
    kt2 = np.load('kt2.npy')
    kt3 = np.load('kt3.npy')
    kt4 = np.load('kt4.npy')
    bmd_npy = np.load("bmd.npy")
    bvtv_npy = np.load("bvtv.npy")
    bsbv_npy = np.load("bsbv.npy")
    tbN_npy = np.load("tbN.npy")
    tbth_npy = np.load("tbth.npy")
    tbsp_npy = np.load("tbsp.npy")
    mil_npy = np.load("mil.npy")
    tv_npy = np.load("tv.npy")
    tmc_npy = np.load("tmc.npy")
    tmd_npy = np.load("tmd.npy")
    bmc_npy = np.load("bmc.npy")


    stp_sk0 = get_stp(sk0)
    stp_sk1 = get_stp(sk1)
    stp_sk2 = get_stp(sk2)
    stp_sk3 = get_stp(sk3)
    stp_sk4 = get_stp(sk4)
    stp_kt0 = get_stp(kt0)
    stp_kt1 = get_stp(kt1)
    stp_kt2 = get_stp(kt2)
    stp_kt3 = get_stp(kt3)
    stp_kt4 = get_stp(kt4)

    print "STP - SK0 : ", stp_sk0
    print "STP - SK1 : ", stp_sk1
    print "STP - SK2 : ", stp_sk2
    print "STP - SK3 : ", stp_sk3
    print "STP - SK4 : ", stp_sk4
    print "STP - KT0 : ", stp_kt0
    print "STP - KT1 : ", stp_kt1
    print "STP - KT2 : ", stp_kt2
    print "STP - KT3 : ", stp_kt3
    print "STP - KT4 : ", stp_kt4

    print "STP - BMD: ", get_stp(bmd_npy)
    print "STP - BMC: ", get_stp(bmc_npy)
    print "STP - TMD: ", get_stp(tmd_npy)
    print "STP - TMC: ", get_stp(tmc_npy)
    print "STP - TV: ", get_stp(tv_npy)
    print "STP - BV/TV: ", get_stp(bvtv_npy)
    print "STP - BS/BV: ", get_stp(bsbv_npy)
    print "STP - Tb.N: ", get_stp(tbN_npy)
    print "STP - Tb.Th: ", get_stp(tbth_npy)
    print "STP - Tb.Sp: ", get_stp(tbsp_npy)
    print "STP - MIL: ", get_stp(mil_npy)



def ltp():

    sk0 = np.load('sk0.npy')
    sk1 = np.load('sk1.npy')
    sk2 = np.load('sk2.npy')
    sk3 = np.load('sk3.npy')
    sk4 = np.load('sk4.npy')
    kt0 = np.load('kt0.npy')
    kt1 = np.load('kt1.npy')
    kt2 = np.load('kt2.npy')
    kt3 = np.load('kt3.npy')
    kt4 = np.load('kt4.npy')
    bmd = np.load("bmd.npy")
    bvtv = np.load("bvtv.npy")
    bsbv = np.load("bsbv.npy")
    tbN = np.load("tbN.npy")
    tbth = np.load("tbth.npy")
    tbsp = np.load("tbsp.npy")
    mil = np.load("mil.npy")
    tv = np.load("tv.npy")
    tmc = np.load("tmc.npy")
    tmd = np.load("tmd.npy")
    bmc = np.load("bmc.npy")

    sk0_xct = np.load('sk0_xct.npy')
    sk1_xct = np.load('sk1_xct.npy')
    sk2_xct = np.load('sk2_xct.npy')
    sk3_xct = np.load('sk3_xct.npy')
    sk4_xct = np.load('sk4_xct.npy')
    kt0_xct = np.load('kt0_xct.npy')
    kt1_xct = np.load('kt1_xct.npy')
    kt2_xct = np.load('kt2_xct.npy')
    kt3_xct = np.load('kt3_xct.npy')
    kt4_xct = np.load('kt4_xct.npy')
    bmd_xct = np.load("bmd_xct.npy")
    bvtv_xct = np.load("bvtv_xct.npy")
    bsbv_xct = np.load("bsbv_xct.npy")
    tbN_xct = np.load("tbN_xct.npy")
    tbth_xct = np.load("tbth_xct.npy")
    tbsp_xct = np.load("tbsp_xct.npy")
    mil_xct = np.load("mil_xct.npy")
    tv_xct = np.load("tv_xct.npy")
    tmc_xct = np.load("tmc_xct.npy")
    tmd_xct = np.load("tmd_xct.npy")
    bmc_xct = np.load("bmc_xct.npy")


    ltp_sk0 = get_ltp(sk0, sk0_xct)
    ltp_sk1 = get_ltp(sk1, sk1_xct)
    ltp_sk2 = get_ltp(sk2, sk2_xct)
    ltp_sk3 = get_ltp(sk3, sk3_xct)
    ltp_sk4 = get_ltp(sk4, sk4_xct)
    ltp_kt0 = get_ltp(kt3, kt0_xct)
    ltp_kt1 = get_ltp(kt3, kt1_xct)
    ltp_kt2 = get_ltp(kt3, kt2_xct)
    ltp_kt3 = get_ltp(kt3, kt3_xct)
    ltp_kt4 = get_ltp(kt3, kt4_xct)
    ltp_bmd = get_ltp(bmd, bmd_xct)
    ltp_bmc = get_ltp(bmc, bmc_xct)
    ltp_tmd = get_ltp(tmd, tmd_xct)
    ltp_tmc = get_ltp(tmc, tmc_xct)
    ltp_tv = get_ltp(tv, tv_xct)
    ltp_bvtv = get_ltp(bvtv, bvtv_xct)
    ltp_mil = get_ltp(mil, mil_xct)
    ltp_tbsp = get_ltp(tbsp, tbsp_xct)
    ltp_tbN = get_ltp(tbN, tbN_xct)
    ltp_tbth = get_ltp(tbth, tbth_xct)
    ltp_bsbv = get_ltp(bsbv, bsbv_xct)


    print "LTP - SK0 : ", ltp_sk0
    print "LTP - SK1 : ", ltp_sk1
    print "LTP - SK2 : ", ltp_sk0
    print "LTP - SK3 : ", ltp_sk1
    print "LTP - SK4 : ", ltp_sk4
    print "LTP - KT0 : ", ltp_kt0
    print "LTP - KT1 : ", ltp_kt1
    print "LTP - KT2 : ", ltp_kt2
    print "LTP - KT3 : ", ltp_kt3
    print "LTP - KT4 : ", ltp_kt4
    print "LTP - BMD : ", ltp_bmd
    print "LTP - BMC : ", ltp_bmc
    print "LTP - TMD : ", ltp_tmd
    print "LTP - TMC : ", ltp_tmc
    print "LTP - TV : ", ltp_tv
    print "LTP - BV/TV : ", ltp_bvtv
    print "LTP - Tb.TH : ", ltp_tbth
    print "LTP - Tb.SP : ", ltp_tbsp
    print "LTP - Tb.N : ", ltp_tbN
    print "LTP - MIL : ", ltp_mil
    print "LTP - BS/BV : ", ltp_bsbv

def precision_and_accuracy():
    mat_dirs = '/home/rodrigo/members.imaglabs.org/felix.thomsen/VertebraPhantom/normalized(MSC)/mats/'

    patients = ["5c", "6b", "8b", "8c", "V12"]
    scans_ltp = ["XCT"]
    scans_stp = ["O1", "O2", "O3", "M1", "M2", "O1", "O2", "O3", "M1", "M2"]
    # amount of volumes of interest
    vois = 4
    sk0 = np.zeros((vois, len(patients), len(scans_stp)))
    sk1 = np.zeros((vois, len(patients), len(scans_stp)))
    sk2 = np.zeros((vois, len(patients), len(scans_stp)))
    sk3 = np.zeros((vois, len(patients), len(scans_stp)))
    sk4 = np.zeros((vois, len(patients), len(scans_stp)))
    kt0 = np.zeros((vois, len(patients), len(scans_stp)))
    kt1 = np.zeros((vois, len(patients), len(scans_stp)))
    kt2 = np.zeros((vois, len(patients), len(scans_stp)))
    kt3 = np.zeros((vois, len(patients), len(scans_stp)))
    kt4 = np.zeros((vois, len(patients), len(scans_stp)))

    #functions = [['MEAN', np.mean], ['MAX', np.max]]

    for i in range(len(patients)):
        for j in range(len(scans_stp)):
            if j < 5:
                filename_slices = mat_dirs + patients[i] + '_' + scans_stp[j] + '_120Slices.mat'
                filename_masks = mat_dirs + patients[i] + '_' + scans_stp[j] + '_120Mask.mat'
                print filename_slices
                print filename_masks

                slices = openMatlab('S', filename_slices)
                masks = openMatlab('M', filename_masks)

                for k in range(1,vois+1):
                    print "k: ", k
                    # voi k
                    voi_k = slices * (masks == k)
                    #np.save(patient + scan + '_voi_k_' + str(k) , voi_k)

                    aux = Stats_MFS_Pyramid_3D()
                    params={'gradient':True}
                    aux.setDef(1, 20, 3, filename_slices, filename_masks, params)
                    stats = aux.getFDs(voi_k)

                    # get those for the paper
                    sk0[(k-1), i, j] = stats[0]
                    kt0[(k-1), i, j] = stats[1]
                    sk1[(k-1), i, j] = stats[2]
                    kt1[(k-1), i, j] = stats[3]
                    sk2[(k-1), i, j] = stats[4]
                    kt2[(k-1), i, j] = stats[5]
                    sk3[(k-1), i, j] = stats[6]
                    kt3[(k-1), i, j] = stats[7]
                    sk4[(k-1), i, j] = stats[8]
                    kt4[(k-1), i, j] = stats[9]

                    np.save("sk0.npy", sk0)
                    np.save("sk1.npy", sk1)
                    np.save("sk2.npy", sk2)
                    np.save("sk3.npy", sk3)
                    np.save("sk4.npy", sk4)
                    np.save("kt0.npy", kt0)
                    np.save("kt1.npy", kt1)
                    np.save("kt2.npy", kt2)
                    np.save("kt3.npy", kt3)
                    np.save("kt4.npy", kt4)

            else:
                filename_slices = mat_dirs + patients[i] + '_' + scans_stp[j] + '_120_140Slices.mat'
                filename_masks = mat_dirs + patients[i] + '_' + scans_stp[j] + '_120_140Mask.mat'
                print filename_slices
                print filename_masks

                slices = openMatlab('S', filename_slices)
                masks = openMatlab('M', filename_masks)

                for k in range(1,vois+1):
                    print "k: ", k
                    # voi k
                    voi_k = slices * (masks == k)

                    aux = Stats_MFS_Pyramid_3D()
                    params={'gradient':True}
                    aux.setDef(1, 20, 3, filename_slices, filename_masks, params)
                    stats = aux.getFDs(voi_k)

                    # get those for the paper
                    sk0[(k-1), i, j] = stats[0]
                    kt0[(k-1), i, j] = stats[1]
                    sk1[(k-1), i, j] = stats[2]
                    kt1[(k-1), i, j] = stats[3]
                    sk2[(k-1), i, j] = stats[4]
                    kt2[(k-1), i, j] = stats[5]
                    sk3[(k-1), i, j] = stats[6]
                    kt3[(k-1), i, j] = stats[7]
                    sk4[(k-1), i, j] = stats[8]
                    kt4[(k-1), i, j] = stats[9]

                    np.save("sk0.npy", sk0)
                    np.save("sk1.npy", sk1)
                    np.save("sk2.npy", sk2)
                    np.save("sk3.npy", sk3)
                    np.save("sk4.npy", sk4)
                    np.save("kt0.npy", kt0)
                    np.save("kt1.npy", kt1)
                    np.save("kt2.npy", kt2)
                    np.save("kt3.npy", kt3)
                    np.save("kt4.npy", kt4)


    stp()


#precision_and_accuracy()
#stp()

def compute_xct():
    mat_dirs = '/home/rodrigo/members.imaglabs.org/felix.thomsen/VertebraPhantom/normalized(MSC)/mats/'

    patients = ["5c", "6b", "8b", "8c", "V12"]
    scans_ltp = ["XCT"]
    vois = 4

    sk0_xct = np.zeros((vois, len(patients)))
    sk1_xct = np.zeros((vois, len(patients)))
    sk2_xct = np.zeros((vois, len(patients)))
    sk3_xct = np.zeros((vois, len(patients)))
    sk4_xct = np.zeros((vois, len(patients)))
    kt0_xct = np.zeros((vois, len(patients)))
    kt1_xct = np.zeros((vois, len(patients)))
    kt2_xct = np.zeros((vois, len(patients)))
    kt3_xct = np.zeros((vois, len(patients)))
    kt4_xct = np.zeros((vois, len(patients)))

    for i in range(len(patients)):
        filename_slices = mat_dirs + patients[i] + '_XCT' + 'Slices.mat'
        filename_masks = mat_dirs + patients[i] + '_XCT' + 'Mask.mat'
        print filename_slices
        print filename_masks

        slices = openMatlab('S', filename_slices)
        masks = openMatlab('M', filename_masks)

        for k in range(1,vois+1):
            print "xct: k: ", k
            voi_k = slices * (masks == k)

            aux = Stats_MFS_Pyramid_3D()
            params={'gradient':True}
            aux.setDef(1, 20, 3, filename_slices, filename_masks, params)
            stats = aux.getFDs(voi_k)

            # get those for the paper
            sk0_xct[(k-1), i] = stats[0]
            kt0_xct[(k-1), i] = stats[1]
            sk1_xct[(k-1), i] = stats[2]
            kt1_xct[(k-1), i] = stats[3]
            sk2_xct[(k-1), i] = stats[4]
            kt2_xct[(k-1), i] = stats[5]
            sk3_xct[(k-1), i] = stats[6]
            kt3_xct[(k-1), i] = stats[7]
            sk4_xct[(k-1), i] = stats[8]
            kt4_xct[(k-1), i] = stats[9]

            np.save("sk0_xct.npy", sk0_xct)
            np.save("sk1_xct.npy", sk1_xct)
            np.save("sk2_xct.npy", sk2_xct)
            np.save("sk3_xct.npy", sk3_xct)
            np.save("sk4_xct.npy", sk4_xct)
            np.save("kt0_xct.npy", kt0_xct)
            np.save("kt1_xct.npy", kt1_xct)
            np.save("kt2_xct.npy", kt2_xct)
            np.save("kt3_xct.npy", kt3_xct)
            np.save("kt4_xct.npy", kt4_xct)

#compute_xct()
#ltp()


def csvToNumpy(X):

    x = tuple(X[0])
    x = np.array(x[1:])

    X_res = x

    for i in range(1, len(X)):
        x = tuple(X[i])
        x = np.array(x[1:])

        X_res = np.vstack((X_res, x))

    return X_res



def repl(s):
    if ',' in s:
        a,b = s.split(',')
        return a+'.'+b
    else:
        return s

# transform a 200 array into a 4*5*10
def make_npy(data):
    data_npy = np.zeros((4, 5, 10))

    for i in range(len(data)):
        data_npy[i%4,int(i/40),int(i/4)%10] = repl(data[i])

    return data_npy

def make_npy_xct(data):
    data_npy = np.zeros((4, 5))

    for i in range(len(data)):
        data_npy[i%4,int(i/4)] = repl(data[i])

    return data_npy

def test_xct():
    measures = recfromcsv('exps/data/xct_only.csv', delimiter=';')
    measures_npy = csvToNumpy(measures)

    sk0_xct = np.zeros((4, 5))

    bmd_xct = measures_npy[:, 3]
    bmc_xct = measures_npy[:, 4]
    tmd_xct = measures_npy[:, 5]
    tmc_xct = measures_npy[:, 6]
    tv_xct = measures_npy[:, 7]
    bvtv_xct = measures_npy[:, 8]
    tbN_xct = measures_npy[:, 9]
    mil_xct =  measures_npy[:, 10]
    tbsp_xct = measures_npy[:, 11]
    tbth_xct = measures_npy[:, 12]
    bsbv_xct = measures_npy[:, 13]

    bmd_xct = make_npy_xct(bmd_xct)
    bmc_xct = make_npy_xct(bmc_xct)
    tmd_xct = make_npy_xct(tmd_xct)
    tmc_xct = make_npy_xct(tmc_xct)
    tv_xct = make_npy_xct(tv_xct)
    bvtv_xct = make_npy_xct(bvtv_xct)
    tbN_xct = make_npy_xct(tbN_xct)
    mil_xct = make_npy_xct(mil_xct)
    tbsp_xct = make_npy_xct(tbsp_xct)
    tbth_xct = make_npy_xct(tbth_xct)
    bsbv_xct = make_npy_xct(bsbv_xct)

    np.save("bmd_xct.npy", bmd_xct)
    np.save("bmc_xct.npy", bmc_xct)
    np.save("tmd_xct.npy", tmd_xct)
    np.save("tmc_xct.npy", tmc_xct)
    np.save("tv_xct.npy", tv_xct)
    np.save("bvtv_xct.npy", bvtv_xct)
    np.save("tbN_xct.npy", tbN_xct)
    np.save("mil_xct.npy", mil_xct)
    np.save("tbsp_xct.npy", tbsp_xct)
    np.save("tbth_xct.npy", tbth_xct)
    np.save("bsbv_xct.npy", bsbv_xct)




    #print measures_npy


def test():
    measures = recfromcsv('exps/data/Parameters_HRQCT_XCT.csv', delimiter=';')
    measures_npy = csvToNumpy(measures)

    #BMD;BMC;TMD;TMC;TV;BV/TV;Tb.N;MIL;Tb.Sp;Tb.Th;BS/BV
    bmd = measures_npy[:, 3]
    bmc = measures_npy[:, 4]
    tmd = measures_npy[:, 5]
    tmc = measures_npy[:, 6]
    tv = measures_npy[:, 7]
    bvtv = measures_npy[:, 8]
    tbN = measures_npy[:, 9]
    mil =  measures_npy[:, 10]
    tbsp = measures_npy[:, 11]
    tbth = measures_npy[:, 12]
    bsbv = measures_npy[:, 13]

    bmd_npy = make_npy(bmd)
    bmc_npy = make_npy(bmc)
    tmd_npy = make_npy(tmd)
    tmc_npy = make_npy(tmc)
    tv_npy = make_npy(tv)
    bvtv_npy = make_npy(bvtv)
    tbN_npy = make_npy(tbN)
    mil_npy = make_npy(mil)
    tbsp_npy = make_npy(tbsp)
    tbth_npy = make_npy(tbth)
    bsbv_npy = make_npy(bsbv)

    np.save("bmd.npy", bmd_npy)
    np.save("bmc.npy", bmc_npy)
    np.save("tmd.npy", tmd_npy)
    np.save("tmc.npy", tmc_npy)
    np.save("tv.npy", tv_npy)
    np.save("bvtv.npy", bvtv_npy)
    np.save("tbN.npy", tbN_npy)
    np.save("mil.npy", mil_npy)
    np.save("tbsp.npy", tbsp_npy)
    np.save("tbth.npy", tbth_npy)
    np.save("bsbv.npy", bsbv_npy)


#precision_and_accuracy()
#compute_xct()
#test_xct()
ltp()
