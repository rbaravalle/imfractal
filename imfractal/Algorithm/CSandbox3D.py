"""
Copyright (c) 2015 Rodrigo Baravalle
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


from Algorithm import *
from random import randrange,randint,seed
from math import log
from scipy import ndimage
import Image
import numpy as np
import scipy
import scipy.stats
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import time, qs3D	
import dicom
import os

class CSandbox3D (Algorithm):

    """

    :3D sandbox multifractal spectrum in cython
    :version: 1.0
    :author: Rodrigo Baravalle
    """

    # how many multifractal dimentions should the algorithm return
    def __init__(self, c):
        self.cant_dims = c

    def setDef(self,x,y,p,params):
        self.total_pixels = params["total_pixels"]  # number of pixels for averaging
        self.v = x
        self.b = y
        self.param = p
        self.params = params

       
    def openData(self, filename):
        return qs3D.volume(self.params, 256, 256)

        # test (should be = 3 for every DF)
        #data = np.ones((256,256,256))
        #return data

    def determine_threshold(self, arr):
        # compute histogram of values
        bins = range(np.min(arr), np.max(arr)+1)

        h = np.histogram(arr, bins = bins)

        threshold = np.min(arr)

        # get x% of mass -> threshold
        assert(len(arr.shape) == 3)

        total_pixels = arr.shape[0] * arr.shape[1] * arr.shape[2]

        for i in range(len(bins)+1) :
            # compute sum of h(x) from x = 0 to x = i
            partial_sum_vector = np.cumsum(h[0][ : (i+1)])
            partial_sum = partial_sum_vector[len(partial_sum_vector) - 1]

            percentage = (float)(partial_sum) / (float)(total_pixels)

            if percentage > 0.75 :
                threshold = np.min(arr) + i
                break

        return threshold

    def openMatlab(self, name, filename, threshold, adaptive = False):

        import scipy.io as sio
        arr = np.array(sio.loadmat(filename)[name]).astype(np.int32)
        if name == "S":
            if adaptive:
                threshold = self.determine_threshold(arr)

            arr = arr > threshold

            a_v = arr.cumsum()

            print "Amount of white pixels: ", a_v[len(a_v) - 1]

        # debug - to see the spongious structure
        # plt.imshow((arr[:,:,50]), cmap=plt.gray())
        # plt.show()

        return arr

    # loads a dicom set of files into a 3d numpy array
    def readDicom(self,path):
        lstFilesDCM = []  # create an empty list
        for dirName, subdirList, fileList in os.walk(path):
            for filename in fileList:
                if ".dcm" in filename.lower():  # check whether the file's DICOM
                    lstFilesDCM.append(os.path.join(dirName,filename))

        # Get ref file
        RefDs = dicom.read_file(lstFilesDCM[0])

        # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
        ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))

        # Load spacing values (in mm)
        ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

        x = numpy.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
        y = numpy.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
        z = numpy.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])

        # The array is sized based on 'ConstPixelDims'
        ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

        print "Loading Dicom..."
        # loop through all the DICOM files
        for filenameDCM in lstFilesDCM:
            # read the file
            ds = dicom.read_file(filenameDCM)
            # store the raw image data
            ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array

        ArrayDicom= numpy.logical_and(ArrayDicom > 3000, ArrayDicom < 6000)
        plt.imshow((ArrayDicom[:,:,220]), cmap=plt.gray())
        plt.show()
        print "loaded!"


        return ArrayDicom


    # get multifractal dimensions
    def getFDs(self,filename, data=[]):

        if len(data) == 0:
            fmask = self.params["mask_filename"]
            threshold = self.params["threshold"]
            adaptive = self.params["adaptive"]
            data = self.openMatlab(self.params["eight"], filename, threshold, adaptive)
            data_mask = self.openMatlab(self.params["nine"], fmask, threshold)

            # Masking
            data = data * (data_mask > 0)

        # debug
        # print "MAX, MIN: ", np.max(data), np.min(data)

        Nx, Ny, Nz = data.shape

        self.P = 30
        P = self.P

        while Nx < 2*P or Ny < 2*P or Nz < 2*P:
            P /= 2
            self.P = P
            print "P too large. New P: ",  P

        L = float(Nx*Ny*Nz)

        t = time.clock()

        Nx, Ny, Nz = data.shape

        white_pixels = np.zeros(((Nx-1-2*P) * (Ny-1-2*P) * (Nz-1-2*P), 3))

        print "Selecting points.."

        h = 0
        for i in range(P, Nx-1-P):
            for j in range(P, Ny-1-P):
                for k in range(P, Nz-1-P):

                    # list with selected points (the points should be in the "structure")
                    # points shouldn't be close to the borders, in order for the windows to have the same size

                    if(data[i][j][k] > 0):
                        white_pixels[h][0] = i
                        white_pixels[h][1] = j
                        white_pixels[h][2] = k
                        h+= 1

        
        white_pixels = white_pixels[:h]
        print "Finished selecting points"
        if len(white_pixels) == 0 :
            print "EMPTY Volume!!!"
            return np.zeros(self.cant_dims, dtype=np.double)

        first_point = white_pixels[randint(0, len(white_pixels)-1)]
        points = np.array([first_point])

        if self.total_pixels >= len(white_pixels):
            # take all points
            points = white_pixels
            self.total_pixels = len(white_pixels)

        else:
            for i in range(self.total_pixels):
                point = randint(0, len(white_pixels) - 1)
                points = np.vstack((points, white_pixels[point]))

                white_pixels = np.delete(white_pixels, point, axis = 0)

        np.set_printoptions(precision=5)
        np.set_printoptions(suppress=True)

        # Summed Area Table
        intImg = data.cumsum(0).cumsum(1).cumsum(2)

        res = qs3D.aux(self.P, self.total_pixels, Nx, Ny, Nz,
                       points.astype(np.int32),
                       np.array(intImg).astype(np.int32),
                       self.cant_dims)

        return res

