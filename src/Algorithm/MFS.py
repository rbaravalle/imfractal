from Algorithm import *
import Image
import numpy as np
from math import exp, log10
import scipy.ndimage.filters as sf
import matplotlib
from matplotlib import pyplot as plt
import scipy.signal

class MFS (Algorithm):

    """


    :version: 1.0
    :author: Rodrigo Baravalle
    """

    def __init__(self):
        """
        @return  :
        @author
        """
        pass

    def setDef(self,ind,f,ite,filen):
        self.ind_num = ind      # number of pixels for averaging
        self.f_num = f             # window
        self.ite_num = ite 
        self.FILENAME = filen


    def gauss_kern(self,size, sizey):
        """ Returns a normalized 2D gauss kernel array for convolutions """
        m = np.float32(size)
        n = np.float32(sizey)
        sigma = 2;     # ???
        if(size <= 3): sigma = 1.5;
        if(size == 5): sigma = 2.5;
        y, x = np.mgrid[-(m-1)/2:(m-1)/2+1, -(n-1)/2:(n-1)/2+1]

        b = 2*(sigma**2)
        x2 = map(lambda i: map( lambda j: j**2,i), x)
        y2 = map(lambda i: map( lambda j: j**2,i), y)
        g = np.sum([x2,y2],axis=0).astype(np.float32)
        g = np.array(map(lambda i: map( lambda j: exp(-j/b),i), g)).astype(np.float32)
        return g / g.sum()


    def getFDs(self, filename):
        """
         

        @param string filename : image location
        @return [float] : multi fractal dimentions
        @author
        """

        #mfs Computes the MFS vector for the input measurement image  im 
        #
        # parameters: ind_num -> determines how many levels are used when computing the density 
        #                          choose 1 for using  directly the image measurement im or
        #                          >= 6 for computing the density of im (quite stable for >=5)      
        #              f_num----> determines the dimension of  MFS vector
        #              ite_num  ---> determines how many levels are used when computing MFS for each level set
        #                            (quite stable for >= 3)
        #
        #MFS = mfs(im) computes the MFS for input  im with default setting
        #                  
        #MFS = mfs(im,ind_num) computes the MFS with ind_num density levels
        #
        #MFS = mfs(im,ind_num, f_num) computes the MFS of dimension f_num for input im 
        #                             with ind_num density levels
        #
        #MFS = mfs(im, ind_num, f_num,ite_num) computes the MFS of dimension f_num for input measurement im
        #                                  using ite_num level iterations in the
        #                                  estimation of the fractal dimension and using ind_num level
        #                                  iterations in the density estimation.
        #
        #Author: Yong Xu, Hui Ji
        #Date: Apr 24, 2007
        #Code ported to python : Rodrigo Baravalle. December 2012

        #self.ind_num = 1
        #if(len(extra) == 1):
        #    self.ind_num = extra[0]  #density counting levels
        #    self.f_num = 26          #the dimension of MFS
        #    self.ite_num = 3         # iteration levels in estimating fractal dimension
        #if(len(extra) == 2):
        #    self.ind_num = extra[0]
        #    self.f_num = extra[1]
        #    self.ite_num = 3
        #if(len(extra) >= 3):
        #    self.ind_num = extra[0]
        #    self.f_num = extra[1]
        #    self.ite_num = extra[2]


        # Extra[3] == True means what we are passing is a filename
        # Extra[3] == False means what we are passing is an array
        #self.FILENAME = extra[3]
        if(self.FILENAME):
            im = Image.open(filename)
            # Preprocessing: if IM is a color image convert it to a gray image 
            im = im.convert("L")
            im = np.array(im.getdata()).reshape(im.size)


        #Using [0..255] to denote the intensity profile of the image
        grayscale_box =[0, 255];

        #Preprocessing: default intensity value of image ranges from 0 to 255
        if(abs(im).max()< 1):
            im = im * grayscale_box[1];

        #######################

        ### Estimating density function of the image
        ### by solving least squares for D in  the equation  
        ### log10(bw) = D*log10(c) + b 
        r = 1.0/max(im.shape)
        c = np.dot(range(1,self.ind_num+1),r)

        c = map(lambda i: log10(i), c)
        bw = np.zeros((self.ind_num,im.shape[0],im.shape[1])).astype(np.float32)

        bw[0] = im + 1

        k = 1
        if(self.ind_num > 1):
            bw[1] = scipy.signal.convolve2d(bw[0], self.gauss_kern(k+1,(k+1)),mode="full")[1:,1:]*((k+1)**2)

        for k in range(2,self.ind_num):
            temp = scipy.signal.convolve2d(bw[0], self.gauss_kern(k+1,(k+1)),mode="full")*((k+1)**2)
            if(k==4):
                bw[k] = temp[k-1-1:temp.shape[0]-(k/2),k-1-1:temp.shape[1]-(k/2)]            
            else:
                bw[k] = temp[k-1:temp.shape[0]-(1),k-1:temp.shape[1]-(1)]


        bw = np.log10(bw)
        n1 = c[0]*c[0]
        n2 = bw[0]*c[0]

        for k in range(1,self.ind_num):
            n1 = n1+c[k]*c[k]
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
            D = im

        #Partition the density
        # throw away the boundary
        D = D[self.ind_num-1:D.shape[0]-self.ind_num+1, self.ind_num-1:D.shape[1]-self.ind_num+1]

        IM = np.zeros(D.shape)
        gap = np.ceil((grayscale_box[1] - grayscale_box[0])/np.float32(self.f_num));
        center = np.zeros(self.f_num);
        for k in range(1,self.f_num+1):
            bin_min = (k-1) * gap;
            bin_max = k * gap - 1;
            center[k-1] = round((bin_min + bin_max) / 2);
            D = ((D <= bin_max) & (D >= bin_min)).choose(D,center[k-1])

        D = ((D >= bin_max)).choose(D,0)
        D = ((D < 0)).choose(D,0)
        IM = D

        #Constructing the filter for approximating log fitting
        r = max(IM.shape)
        c = np.zeros(self.ite_num)
        c[0] = 1;
        for k in range(1,self.ite_num):
            c[k] = c[k-1]/(k+1)
        c = c / sum(c);

        #Construct level sets
        Idx_IM = np.zeros(IM.shape);
        for k in range(0,self.f_num):
            IM = (IM == center[k]).choose(IM,k+1)

        Idx_IM = IM
        IM = np.zeros(IM.shape)

        #Estimate MFS by box-counting
        num = np.zeros(self.ite_num)
        MFS = np.zeros(self.f_num)
        for k in range(1,self.f_num+1):
            IM = np.zeros(IM.shape)
            IM = (Idx_IM==k).choose(Idx_IM,255+k)
            IM = (IM<255+k).choose(IM,0)
            IM = (IM>0).choose(IM,1)
            temp = max(IM.sum(),1)
            num[0] = log10(temp)/log10(r);    
            for j in range(2,self.ite_num+1):
                mask = np.ones((j,j))
                bw = scipy.signal.convolve2d(IM, mask,mode="full")[1:,1:]
                indx = np.arange(0,IM.shape[0],j)
                indy = np.arange(0,IM.shape[1],j)
                bw = bw[np.ix_(indx,indy)]
                idx = (bw>0).sum()
                temp = max(idx,1)
                num[j-1] = log10(temp)/log10(r/j)

            MFS[k-1] = sum(c*num)

        return MFS

