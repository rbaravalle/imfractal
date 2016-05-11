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


from Algorithm import *
import numpy as np
from math import log10
import scipy.signal
import scipy.io as sio


class MFS_3D (Algorithm):

    """
    :3D implementation of MFS through holder exponents f(alpha)
    :version: 1.0
    :author: Rodrigo Baravalle
    """

    def __init__(self):
        pass

    def setDef(self, ind, f, ite, filename, file_mask):

        # parameters: ind -> determines how many levels are used when computing the density 
        #                          choose 1 for using  directly the image measurement im or
        #                          >= 6 for computing the density of im (quite stable for >=5)      
        #              f ----> determines the dimension of  MFS vector
        #              ite  ---> determines how many levels are used when computing MFS for each 

        self.ind_num = ind      # number of pixels for averaging
        self.f_num = f          # window
        self.ite_num = ite
        self.filename = filename
        self.file_mask = file_mask



    def gauss_kern(self,size_x, size_y, size_z):
        """ Returns a normalized 3D gauss kernel array for convolutions """
        m = np.float32(size_x)
        n = np.float32(size_y)
        o = np.float32(size_z)
        sigma = 2;     # ???
        if(size_x <= 3): sigma = 1.5;
        if(size_x == 5): sigma = 2.5;
        z, y, x = np.mgrid[-(m-1)/2:(m-1)/2+1, -(n-1)/2:(n-1)/2+1, -(o-1)/2:(o-1)/2+1]

        b = 2*(sigma**2)
        square = lambda i : i**2
        fm = lambda i: map(square, i)

        x2 = map(fm, x)
        y2 = map(fm, y)
        z2 = map(fm, z)

        g = np.sum([x2, y2, z2], axis=0).astype(np.float32)
        g = np.exp(g).astype(np.float32)
        return g / g.sum()

    def determine_threshold(self, arr):
        # compute histogram of values
        bins = range(np.min(arr), np.max(arr) + 1)

        h = np.histogram(arr, bins=bins)

        threshold = np.min(arr)

        # get x% of mass -> threshold
        assert (len(arr.shape) == 3)

        total_pixels = arr.shape[0] * arr.shape[1] * arr.shape[2]

        for i in range(len(bins) + 1):
            # compute sum of h(x) from x = 0 to x = i
            partial_sum_vector = np.cumsum(h[0][: (i + 1)])
            partial_sum = partial_sum_vector[len(partial_sum_vector) - 1]

            percentage = (float)(partial_sum) / (float)(total_pixels)

            if percentage > 0.75:
                threshold = np.min(arr) + i
                break

        return threshold

    #def openMatlab(self, name, filename):

        #arr = np.array(sio.loadmat(filename)[name]).astype(np.int32)
        #return arr

    def openMatlab(self, name, filename):

        import scipy.io as sio
        arr = np.array(sio.loadmat(filename)[name]).astype(np.int32)
        if name == "S":
            threshold = self.determine_threshold(arr)

            arr = arr > threshold

            a_v = arr.cumsum()

            print "Amount of white pixels: ", a_v[len(a_v) - 1]

        # debug - to see the spongious structure
        # plt.imshow((arr[:,:,50]), cmap=plt.gray())
        # plt.show()

        return arr


    def getFDs(self):
        """
        @param string filename : volume location
        @param string file_mask : mask volume location
        @return [float] : 3D multi fractal dimentions
        @author: Rodrigo Baravalle. Code ported from Matlab and extended to 3D
        """

        # data is a 3D grayscale volume
        data = self.openMatlab('S', self.filename)
        data_mask = self.openMatlab('M', self.file_mask)

        # Masking
        data = data * (data_mask > 0)

        #Using [0..255] to denote the intensity profile of the image
        grayscale_box = [0, 255]

        #Preprocessing: default intensity value of image ranges from 0 to 255
        if abs(data).max()< 1:
            data = data * grayscale_box[1]
        else:
            # put every value into [0, 255]
            data = (data - data.min()) * 255 / (data.max() - data.min())
        #######################

        #DEBUG
        print data.max(), data.min(), data.sum()

        ### Estimating density function of the volume
        ### by solving least squares for D in  the equation  
        ### log10(bw) = D*log10(c) + b 
        r = 1.0 / max(data.shape)
        c = np.dot(range(1, self.ind_num+1), r)


        c = map(lambda i: log10(i), c)

        bw = np.zeros((self.ind_num, data.shape[0], data.shape[1], data.shape[2])).astype(np.float32)

        bw[0] = data + 1

        # DEBUG
        print "BW: ", bw.shape

        k = 1
        if(self.ind_num > 1):
            bw[1] = scipy.signal.convolve(bw[0], self.gauss_kern(k+1, k+1, k+1), mode="full")[1:,1:]*((k+1)**2)

        for k in range(2,self.ind_num):
            temp = scipy.signal.convolve(bw[0], self.gauss_kern(k+1, k+1, k+1), mode="full")*((k+1)**2)
            if(k==4):
                bw[k] = temp[k - 1 - 1 : temp.shape[0] - (k / 2),
                             k - 1 - 1 : temp.shape[1] - (k / 2),
                             k - 1 - 1 : temp.shape[2] - (k / 2)]
            else:
                bw[k] = temp[k - 1 : temp.shape[0] - (1),
                             k - 1 : temp.shape[1] - (1),
                             k - 1 : temp.shape[2] - (1)]


        print bw.min(), bw.max()
        bw = np.log10(bw)
        n1 = c[0] * c[0]
        n2 = bw[0] * c[0]

        for k in range(1,self.ind_num):
            n1 = n1 + c[k]*c[k]
            n2 = n2 + bw[k]*c[k]

        sum3 = bw[0]
        for i in range(1,self.ind_num):
            sum3 = sum3 + bw[i]

        if(self.ind_num >1):
            D = (n2*self.ind_num-sum(c)*sum3)/(n1*self.ind_num -sum(c)*sum(c));

        if (self.ind_num > 1):
            max_D  = np.float32(4)
            min_D = np.float32(1)
            D = grayscale_box[1]*(D-min_D)/(max_D - min_D)+grayscale_box[0]
        else:
            D = data

        #Partition the density
        # throw away the boundary
        D = D[self.ind_num - 1 : D.shape[0] - self.ind_num + 1,
              self.ind_num - 1 : D.shape[1] - self.ind_num + 1,
              self.ind_num - 1 : D.shape[2] - self.ind_num + 1]

        IM = np.zeros(D.shape)
        gap = np.ceil((grayscale_box[1] - grayscale_box[0])/np.float32(self.f_num));
        center = np.zeros(self.f_num);
        for k in range(1,self.f_num+1):
            bin_min = (k-1) * gap;
            bin_max = k * gap - 1;
            center[k-1] = round((bin_min + bin_max) / 2);
            D = ((D <= bin_max) & (D >= bin_min)).choose(D, center[k-1])

        D = ((D >= bin_max)).choose(D,0)
        D = ((D < 0)).choose(D,0)
        IM = D


        # Constructing the filter for approximating log fitting
        r = max(IM.shape)
        c = np.zeros(self.ite_num)
        c[0] = 1;
        for k in range(1,self.ite_num):
            c[k] = c[k-1]/(k+1)
        c = c / sum(c);

        # Construct level sets
        Idx_IM = np.zeros(IM.shape);
        for k in range(0, self.f_num):
            IM = (IM == center[k]).choose(IM,k+1)

        Idx_IM = IM
        IM = np.zeros(IM.shape)


        #Estimate MFS by box-counting
        num = np.zeros(self.ite_num)
        MFS = np.zeros(self.f_num)
        for k in range(1, self.f_num+1):
            print k, self.f_num
            IM = np.zeros(IM.shape)
            IM = (Idx_IM == k).choose(Idx_IM, 255 + k)
            IM = (IM<255 + k).choose(IM, 0)
            IM = (IM > 0).choose(IM, 1)
            temp = max(IM.sum(), 1)
            num[0] = log10(temp)/log10(r);    
            for j in range(2, self.ite_num+1):
                mask = np.ones((j, j, j))
                bw = scipy.signal.convolve(IM, mask, mode = "full")[1:, 1:, 1:]
                ind_x = np.arange(0, IM.shape[0], j)
                ind_y = np.arange(0, IM.shape[1], j)
                ind_z = np.arange(0, IM.shape[2], j)
                bw = bw[np.ix_(ind_x, ind_y, ind_z)]
                idx = (bw > 0 ).sum()
                temp = max(idx, 1)
                num[j-1] = log10( temp ) / log10( r / j )

            MFS[k-1] = sum(c*num)

        return MFS

