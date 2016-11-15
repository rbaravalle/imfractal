import Image
import numpy as np
import os
from imfractal import *
import time
import matplotlib.pyplot as plt

def load_synthetic_volume(mask):
    im = np.asarray(Image.open('/home/rodrigo/result/porous.tga0').convert("L"))

    num_imgs = len(os.listdir('/home/rodrigo/result/'))

    arr = np.zeros((num_imgs, im.shape[0], im.shape[1])).astype(np.uint8)
    

    for i in range(num_imgs):
        arr[i,:,:] = np.asarray(Image.open('/home/rodrigo/result/porous.tga'+str(i)).convert("L"))



    #padding, mask_ should be bigger than mask
    print arr.shape
    print mask.shape
    mask_ = np.zeros(arr.shape)
    mask_[:mask.shape[0], :mask.shape[1], :mask.shape[2]] = mask

    arr = (255-arr) * (mask_ > 0)

    plt.imshow(Image.fromarray(arr[num_imgs/2]))
    plt.show()
    np.save('synth.npy', arr)

    return arr


def openMatlab(name, filename, threshold, adaptive = False):

    import scipy.io as sio
    arr = np.array(sio.loadmat(filename)[name]).astype(np.uint8)
    if name == "S":
        #if adaptive:
        #    threshold = self.determine_threshold(arr)

        arr = arr * (arr > threshold)

        a_v = arr.cumsum()

        print "Amount of white pixels: ", a_v[len(a_v) - 1]

    # debug - to see the spongious structure
    plt.imshow((arr[50,:,:]), cmap=plt.gray())
    plt.show()

    return arr

def getFDs(filename, fmask):
    threshold = 200 # self.params["threshold"]
    data = openMatlab("S", filename, threshold)
    data_mask = openMatlab("M", fmask, threshold)

    # Masking
    data = data * (data_mask > 0)

    np.save('bone.npy', data)

    return data, data_mask

def main():

    fname = '/home/rodrigo/Documentos/members.imaglabs.org/felix.thomsen/Rodrigo/BioAsset/mats/BA01_120_1Slices.mat'
    fmask = '/home/rodrigo/Documentos/members.imaglabs.org/felix.thomsen/Rodrigo/BioAsset/mats/BA01_120_1Mask.mat'

    data, mask = getFDs(fname, fmask)
    synth = load_synthetic_volume(mask)

    print "Data: ", data.shape
    print "Synthetic: ", synth.shape

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

    aux = CSandbox3D(21)
    aux.setDef(40, 1.02, True, params)

    plt.title('Comparison')
    plt.ylim((2.5, 3.5))

    print "Computing Sandbox 3D Multifractal Spectrum... (real)"
    t =  time.clock()
    fds = aux.getFDs('', data)
    t =  time.clock()-t
    print "Time 3D MFS: ", t
    print fds

    plt.plot(fds,'-', label = 'Real')

    print "Computing Sandbox 3D Multifractal Spectrum... (synthetic)"
    t =  time.clock()
    fds = aux.getFDs('', synth)
    t =  time.clock()-t
    print "Time 3D MFS: ", t
    print fds


    plt.plot(fds, 'x', label = 'Synth')

    plt.legend()
    plt.show()

main()
