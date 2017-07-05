import numpy as np
import scipy.io as sio
from imfractal import *

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
    sk4 = np.load('sk4.npy')
    kt3 = np.load('kt3.npy')

    stp_sk0 = get_stp(sk0)
    stp_sk1 = get_stp(sk1)
    stp_sk4 = get_stp(sk4)
    stp_kt3 = get_stp(kt3)

    print "STP - SK0 : ", stp_sk0
    print "STP - SK1 : ", stp_sk1
    print "STP - SK4 : ", stp_sk4
    print "STP - KT3 : ", stp_kt3

def ltp():
    pass

def precision_and_accuracy():
    mat_dirs = '/home/rodrigo/members.imaglabs.org/felix.thomsen/VertebraPhantom/normalized(MSC)/mats/'

    patients = ["5c", "6b", "8b", "8c", "V12"]
    scans_ltp = ["01", "02", "03", "M1", "M2", "XCT"]
    scans_stp = ["O1", "O2", "O3", "M1", "M2", "O1", "O2", "O3", "M1", "M2"]
    # amount of volumes of interest
    vois = 4
    sk0 = np.zeros((vois, len(patients), len(scans_stp)))
    sk1 = np.zeros((vois, len(patients), len(scans_stp)))
    sk4 = np.zeros((vois, len(patients), len(scans_stp)))
    kt3 = np.zeros((vois, len(patients), len(scans_stp)))

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
                    sk1[(k-1), i, j] = stats[2]
                    sk4[(k-1), i, j] = stats[8]
                    kt3[(k-1), i, j] = stats[7]

                    np.save("sk0.npy", sk0)
                    np.save("sk1.npy", sk1)
                    np.save("sk4.npy", sk4)
                    np.save("kt3.npy", kt3)

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
                    sk1[(k-1), i, j] = stats[2]
                    sk4[(k-1), i, j] = stats[8]
                    kt3[(k-1), i, j] = stats[7]

                    np.save("sk0.npy", sk0)
                    np.save("sk1.npy", sk1)
                    np.save("sk4.npy", sk4)
                    np.save("kt3.npy", kt3)


    stp()


#precision_and_accuracy()
stp()
