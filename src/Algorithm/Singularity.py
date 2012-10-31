from Algorithm import *
from random import randrange,randint
from math import log
from scipy import ndimage
from pylab import plot, title, show , legend
import matplotlib
from matplotlib import pyplot as plt
import time
import Image
import numpy as np
import sys
import os

class Singularity (Algorithm):

    """

    :singularity multifractal spectrum
    :version: 1.0
    :author: Rodrigo Baravalle
    """

    def setCuantas(self,x):
        self.cuantas = x

    def getFDs(self,filename):
        t = time.clock()
        cantSelected = 0

        a = Image.open(filename)
        Nx, Ny = a.size
        L = Nx*Ny

        points = []     # number of elements in the structure
        gray = a.convert('L') # rgb 2 gray

        #plt.imshow(gray, cmap=matplotlib.cm.gray)
        #plt.show()

        alphaIm = np.zeros((Nx,Ny), dtype=np.double ) # Nx rows x Ny columns
        #measure = np.zeros(4, dtype=np.double ) # Ny rows x 4 columns

        l = 4 # (maximum window size-1) / 2
        # from 1 to l
        temp = np.log((1.0,3.0,5.0,7.0));
        measure = np.zeros(4);

        for i in range(0,Nx-1):
            for j in range(0,Ny-1):
                measure[0] = max(gray.crop((max(i-1,0),max(j-1,0),min(i+1,Nx-1),min(j+1,Ny-1))).getdata()) + 1
                measure[1] = max(gray.crop((max(i-2,0),max(j-2,0),min(i+2,Nx-1),min(j+2,Ny-1))).getdata()) + 1
                measure[2] = max(gray.crop((max(i-3,0),max(j-3,0),min(i+3,Nx-1),min(j+3,Ny-1))).getdata()) + 1
                measure[3] = max(gray.crop((max(i-4,0),max(j-4,0),min(i+4,Nx-1),min(j+4,Ny-1))).getdata()) + 1
                alphaIm[j,i] = np.polyfit(temp,np.log(measure),1)[0]

        maxim = np.max(alphaIm)
        minim = np.min(alphaIm)

        # Alpha image
        #plt.imshow(alphaIm, cmap=matplotlib.cm.gray)
        #plt.show()

        paso = (maxim-minim)/self.cuantas
        clases = np.arange(minim,maxim,paso)
        print "Clases: ", clases, self.cuantas

        # Window
        cant = int(np.floor(np.log(Nx)))

        hs = np.hstack((alphaIm,alphaIm))
        alphaIm = np.vstack((hs,hs))

        # Multifractal dimentions
        falpha = np.zeros(self.cuantas)

        for c in range(self.cuantas):
            N = np.zeros(cant+1)
            for k in range(cant+1):
                sizeBlocks = 2*k+1
                numBlocks_x = int(np.ceil(Nx/sizeBlocks))
                numBlocks_y = int(np.ceil(Ny/sizeBlocks))

                flag = np.zeros((numBlocks_x,numBlocks_y))

                for i in range(1,numBlocks_x):
                    for j in range(1,numBlocks_y):
                        xi = (i-1)*sizeBlocks
                        xf = i*sizeBlocks-1
                        yi = (j-1)*sizeBlocks
                        yf = j*sizeBlocks-1
                        if(xf == xi): xf = xf+1
                        if(yf == yi): yf = yf+1
                        block = alphaIm[xi : xf, yi : yf]

                        f = 0;
                        s1 = len(block)
                        s2 = len(block[0])
                        if(c != self.cuantas-1):
                            for w in range(s1):
                                for t in range(s2):
                                    b = block[w,t]
                                    if (b >= clases[c] and b < clases[c+1]):
                                       f = 1
                                    if(f == 1):
                                        break
                                if(f == 1):
                                    break
                        else:
                            for w in range(s1):
                                for t in range(s2):
                                    b = block[w,t]
                                    if (b == clases[c]+1): # !!
                                       f = 1
                                    if(f == 1):
                                        break
                                if(f == 1):
                                    break
                        
                        flag[i-1,j-1] = f
                        N[k] = N[k] + f;

            falpha[c] = -np.polyfit(map(lambda i: np.log(i*2+1),range(cant+1)),np.log(map(lambda i: i+1,N)),1)[0]

        t =  time.clock()-t
        s = np.hstack((clases,falpha))
        print "Time: ", t
        return s

