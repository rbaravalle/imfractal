import Image
import numpy as np
import os
from imfractal import *

def load_synthetic_volume():
    im = np.asarray(Image.open('/home/rodrigo/result/porous.tga0').convert("L"))

    num_imgs = len(os.listdir('/home/rodrigo/result/'))

    arr = np.zeros((num_imgs, im.shape[0], im.shape[1])).astype(np.uint8)
    

    for i in range(num_imgs):
        arr[:,i,:] = np.asarray(Image.open('/home/rodrigo/result/porous.tga'+str(i)).convert("L"))

#    plt.imshow(Image.fromarray(arr[num_imgs/2]))
#    plt.show()
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
    # plt.imshow((arr[:,:,50]), cmap=plt.gray())
    # plt.show()

    return arr

def getFDs(filename, fmask):
    threshold = 200 # self.params["threshold"]
    data = openMatlab("S", filename, threshold)
    data_mask = openMatlab("M", fmask, threshold)

    # Masking
    data = data * (data_mask > 0)

    np.save('bone.npy', data)

    return data

def main():

    fname = '/home/rodrigo/Documentos/members.imaglabs.org/felix.thomsen/Rodrigo/BioAsset/mats/BA01_120_1Slices.mat'
    fmask = '/home/rodrigo/Documentos/members.imaglabs.org/felix.thomsen/Rodrigo/BioAsset/mats/BA01_120_1Mask.mat'

    data = getFDs(fname, fmask)
    synth = load_synthetic_volume()

    print "Data: ", data.shape
    print "Synthetic: ", synth.shape

main()
