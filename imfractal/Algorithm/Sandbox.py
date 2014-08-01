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

class Sandbox (Algorithm):

    """

    :sandbox multifractal spectrum
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
    # img: original image
    # Nx: img size(x)
    # Ny: img size(y)
    # which: type of img: 
    #   'img': Image
    #   else:  array
    def sat(self,img,Nx,Ny,which):
        # summed area table
        intImg = np.zeros((Nx,Ny))
        if(which == 'img'):
            intImg[0][0] = img.getpixel((0,0))
            
            arrNx = range(1,Nx)
            arrNy = range(1,Ny)
            for h in arrNx:
                intImg[h][0] = intImg[h-1][0] + img.getpixel((h,0))
            
            for w in arrNy:
                intImg[0][w] = intImg[0][w-1] + img.getpixel((0,w))
            
            for f in arrNx:
                for g in arrNy:
                   intImg[f][g] = img.getpixel((f,g))+intImg[f-1][g]+intImg[f][g-1]-intImg[f-1][g-1]
        else:
            intImg[0][0] = img[0][0]
            
            arrNx = range(1,Nx)
            arrNy = range(1,Ny)
            for h in arrNx:
                intImg[h][0] = intImg[h-1][0] + img[h][0]
            
            for w in arrNy:
                intImg[0][w] = intImg[0][w-1] + img[0][w]
            
            for f in arrNx:
                for g in arrNy:
                   intImg[f][g] = img[f][g]+intImg[f-1][g]+intImg[f][g-1]-intImg[f-1][g-1]

        return intImg

    # white's algorithm
    # local threshold schema
    def white(self,img,Nx,Ny):
               
        im = np.zeros((Nx,Ny))
        
        intImg = self.sat(img,Nx,Ny,'img')
            
        arrNx = range(Nx)
        arrNy = range(Ny)

        vent = int(self.v)
        for i in arrNx:
            for j in arrNy:
                if(self.mww(max(0,i-vent),max(0,j-vent),min(Nx-1,i+vent),min(Ny-1,j+vent),intImg) >= img.getpixel((i,j))*self.b ): 
                    v = img.getpixel((i,j))
                    if (v > 0):
                        im[i,j] = 255#img.getpixel((i,j))

        # do an opening operation to remove small elements
        return im.T #ndimage.binary_opening(im.T, structure=np.ones((0,0))).astype(np.int)


    # sum of values in the region (x1,y1), (x2,y2) in intImg
    # intImg: summed area table
    def count(self,x1,y1,x2,y2,intImg):
        sum = intImg[x2][y2]
        if (x1>= 1 and y1 >= 1):
            sum = sum + intImg[x1-1][y1-1]
        if (x1 >= 1):
            sum = sum - intImg[x1-1][y2];
        if (y1 >= 1):
            sum = sum - intImg[x2][y1-1]
        return sum
                
    # get multifractal dimensions
    def getFDs(self,filename):
        startR = 1
        cantSelected = 0

        a = Image.open(filename)
        Nx, Ny = a.size

        self.P = min(Nx,Ny)/6

        L = float(Nx*Ny)

        points = []     # number of elements in the structure
        if(self.param):
            gray = a.convert('L') # rgb 2 gray

            gray = self.white(gray,Nx,Ny) # local thresholding algorithm
        else: 
            b = a.getdata()
            if(type(b[0]) is int): a=b
            else: a = np.array(map (lambda i: i[0], np.array(b))) # argh!
            gray = np.array(a,np.uint8).reshape(b.size[1],b.size[0])
        #plt.imshow(gray, cmap=matplotlib.cm.gray)
        #plt.show()

        Nx = gray.shape[0]
        Ny = gray.shape[1]

        intImg = self.sat(gray,Nx,Ny,'array')

        m0 = intImg[Nx-1][Ny-1]

        l = range(-self.cant,self.cant+1)
        res = np.zeros(len(l), dtype=np.double )

        if(m0 == 0):
            print "Empty IMAGE structure!!!"
            return res

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

        # R's
        stepR = 1
        rvalues = np.array(range(startR,self.P+startR,stepR))
        # ln (R/L)
        sizes = np.log(rvalues/float(Nx))

        h = 0
        for q in l:
            down = 1.0/(m0**(q-1))

            c = np.zeros((len(rvalues)), dtype=np.double )
            # ln< M(R)/M0 ** q-1 >
            if(q!=1):
                ind = 0
                for R in rvalues:
                    summ = np.double(0)
                    for i in range(self.total): # mean
                        x = points[i][0]
                        y = points[i][1]
                        MR = self.count(np.clip(x-R,0,Nx-1),np.clip(y-R,0,Ny-1),np.clip(x+R,0,Nx-1),np.clip(y+R,0,Ny-1),intImg)
                        summ+= down*((MR)**(q-1))
                    summ /= float(self.total) # mean
                    c[ind] = np.log(summ)
                    ind+=1

                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(sizes,c/float(q-1))
                plt.plot(sizes,c/float(q-1))

                print slope,"Residual for q = ", q," : " , "scipy:", r_value**2

            else: 
                # q = 1
                #< ln(M(R)/M0) >
                ind = 0
                for R in rvalues:
                    summ = np.double(0)
                    for i in range(self.total): # mean
                        x = points[i][0]
                        y = points[i][1]
                        MR = self.count(np.clip(x-R,0,Nx-1),np.clip(y-R,0,Ny-1),np.clip(x+R,0,Nx-1),np.clip(y+R,0,Ny-1),intImg)
                        summ += np.log(MR/m0)
                    summ /= float(self.total) # mean
                    c[ind] = summ
                    ind+=1

                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(sizes,c)
                plt.plot(sizes,c)
                print slope,"Residual for q = ", q," : " , "scipy:", r_value**2

            res[h] = slope
            h+=1

        plt.show()
        #print res.shape
        return res

