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
from random import randrange,randint
from math import log
from scipy import ndimage
from pylab import plot, title, show , legend
import matplotlib
from matplotlib import pyplot as plt
import Image
import numpy as np
import sys
import os
from scipy.sparse import coo_matrix
import scipy.sparse.linalg as linsolve

class Singularity (Algorithm):

    """

    :singularity multifractal spectrum
    :version: 1.1
    :author: Rodrigo Baravalle
    """

    # how many multifractal dimentions should thw algorithm return
    def __init__(self, c):
        self.cuantas = c

    def getFDs(self,filename):
        a = Image.open(filename)
        Nx, Ny = a.size
        L = Nx*Ny

        points = []     # number of elements in the structure
        gray = a.convert('L') # rgb 2 gray

        alphaIm = np.zeros((Nx,Ny), dtype=np.double ) # Nx rows x Ny columns
        #measure = np.zeros(4, dtype=np.double ) # Ny rows x 4 columns

        l = 4 # (maximum window size-1) / 2
        temp = np.log((1.0,3.0,5.0,7.0))
        measure = np.zeros(l*Ny)

        b = np.vstack((temp,np.ones((1,l)))).T
        AA=coo_matrix(np.kron(np.identity(Nx), b))
        
        # Instead of doing polyfits, a sparse linear system is constructed and solved
        for i in range(Nx):
            for j in range(Ny):
                measure[0+j*l] = max(gray.crop((max(i-1,0),max(j-1,0),min(i+1,Nx-1),min(j+1,Ny-1))).getdata()) + 1
                measure[1+j*l] = max(gray.crop((max(i-2,0),max(j-2,0),min(i+2,Nx-1),min(j+2,Ny-1))).getdata()) + 1
                measure[2+j*l] = max(gray.crop((max(i-3,0),max(j-3,0),min(i+3,Nx-1),min(j+3,Ny-1))).getdata()) + 1
                measure[3+j*l] = max(gray.crop((max(i-4,0),max(j-4,0),min(i+4,Nx-1),min(j+4,Ny-1))).getdata()) + 1

            bb=np.log(measure)
            z = linsolve.lsqr(AA,bb)[0]
            z = z.reshape(2,Ny,order = 'F')
            alphaIm[i] = z[0]

        maxim = np.max(alphaIm)
        minim = np.min(alphaIm)
        alphaIm = alphaIm.T

        # Alpha image
        #plt.imshow(alphaIm, cmap=matplotlib.cm.gray)
        #plt.show()

        paso = (maxim-minim)/self.cuantas
        clases = np.arange(minim,maxim,paso)

        # Window
        cant = int(np.floor(np.log(Nx)))

        # concatenate the image A as [[A,A],[A,A]]
        hs = np.hstack((alphaIm,alphaIm))
        alphaIm = np.vstack((hs,hs))

        # Multifractal dimentions
        falpha = np.zeros(self.cuantas)

        for c in range(self.cuantas):
            N = np.zeros(cant+1)
            # window sizes
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
                            # f = 1 if any pixel in block is between clases[c] and clases[c+1]
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
                            # f = 1 if any pixel in block is equal to classes[c]+1
                            for w in range(s1):
                                for t in range(s2):
                                    b = block[w,t]
                                    if (b == clases[c]): # !!
                                       raise SystemError
                                       f = 1
                                    if(f == 1):
                                        break
                                if(f == 1):
                                    break
                        
                        flag[i-1,j-1] = f

                        # number of blocks with holder exponents for this class (c)
                        # and for this window size (k)
                        N[k] = N[k] + f;

            # Haussodorf (box) dimention of the alpha distribution
            falpha[c] = -np.polyfit(map(lambda i: np.log(i*2+1),range(cant+1)),np.log(map(lambda i: i+1,N)),1)[0]

        s = np.hstack((clases,falpha))
        return s

