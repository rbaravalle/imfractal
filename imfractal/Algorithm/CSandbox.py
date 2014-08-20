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
import time, qs

class CSandbox (Algorithm):

    """

    :sandbox multifractal spectrum in cython
    :version: 1.0
    :author: Rodrigo Baravalle
    """

    # how many multifractal dimentions should the algorithm return
    def __init__(self, c):
        self.cant = c

    def setDef(self,x,y,p):
        self.total = 1000#*3      # number of pixels for averaging
        self.v = x
        self.b = y
        self.param = p

    # returns the sum of (summed area) image pixels in the box between
    # (x1,y1) and (x2,y2)
    def mww(self,x1,y1,x2,y2,intImg):
        sum = intImg[x2][y2]
        if (x1>= 1 and y1 >= 1):
            sum = sum + intImg[x1-1][y1-1]
        if (x1 >= 1):
            sum = sum - intImg[x1-1][y2];
        if (y1 >= 1):
            sum = sum - intImg[x2][y1-1]
        return sum/((x2-x1+1)*(y2-y1+1));

    # constructs summed area table
    def sat(self,img,Nx,Ny):
        intImg = np.empty((Nx,Ny))
            
        intImg[0][0] = img[0][0]
        intImg[1:,0] = intImg[0:-1,0] + img[1:,0]           
        intImg[0,1:] = intImg[0,0:-1] + img[0,1:]
        
        for f in range(1,Nx):
            for g in range(1,Ny):
               intImg[f][g] = img[f][g]+intImg[f-1][g]+intImg[f][g-1]-intImg[f-1][g-1]

        return intImg

    # white's algorithm
    # local threshold schema
    def white(self,img,Nx,Ny):             
        im = np.zeros((Nx,Ny))
        intImg = self.sat(np.asarray(img).astype(np.int32),Nx,Ny)
        vent = int(self.v)
        for i in range(Nx):
            for j in range(Ny):
                if(self.mww(max(0,i-vent),max(0,j-vent),min(Nx-1,i+vent),min(Ny-1,j+vent),intImg) >= img[i,j]*self.b ): 
                    v = img[i,j]
                    if (v > 0):
                        im[i,j] = 255
        return im.T

    # get multifractal dimensions
    def getFDs(self,filename):
        cantSelected = 0

        a = Image.open(filename)
        Nx, Ny = a.size

        self.P = min(Nx,Ny)/6

        L = float(Nx*Ny)

        t = time.clock()
        points = []     # number of elements in the structure
        if(self.param):
            gray = a.convert('L') # rgb 2 gray

            gray = np.asarray(gray).T.astype(np.int32)
            t = time.clock()
            gray = self.white(gray,Nx,Ny).T # local thresholding algorithm
            print "Time white :", time.clock()-t
        else: 
            b = a.getdata()
            if(type(b[0]) is int): a=b
            else: a = np.array(map (lambda i: i[0], np.array(b))) # argh!
            gray = np.array(a).reshape(b.size[1],b.size[0])
        #plt.imshow(gray, cmap=matplotlib.cm.gray)
        #plt.show()

        Nx, Ny = gray.shape

        intImg = self.sat(np.array(gray).astype(np.int32),Nx,Ny)

        m0 = intImg[Nx-1][Ny-1]

        if(m0 == 0):
            print "Empty IMAGE structure!!!"
            return np.zeros(self.cant*2+1, dtype=np.double )

        if(m0 < self.total):
            print "Warning: structure with less points than expected"
            self.total = m0/2 # FIX ME
            

        x = randint(0,Nx-1)
        y = randint(0,Ny-1)
        while(gray[x][y] == 0):
            x = randint(0,Nx-1)
            y = randint(0,Ny-1)

        # list with selected points (the points should be in the "structure")
        # points shouldn't be close to the borders, in order for the windows to have the same size

        while cantSelected < self.total:
            while(([x,y] in points) or gray[x][y] == 0):
                x = randint(0,Nx-1)
                y = randint(0,Ny-1)
            # new point, add to list
            points.append([x,y])
            cantSelected = cantSelected+1

        np.set_printoptions(precision=5)
        np.set_printoptions(suppress=True)

        points = np.array(points).astype(np.int32)

        res = qs.aux(self.P,self.total,Nx,Ny,points,np.array(intImg).astype(np.int32),m0,self.cant)

        return res

