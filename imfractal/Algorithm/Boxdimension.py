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
#from pylab import plot, title, show , legend
import matplotlib
from matplotlib import pyplot as plt
import Image
import numpy as np
import sys
import os

class Boxdimension (Algorithm):
    
    """
    :Box dimension
    :version: 1.0
    :author: Rodrigo Baravalle
    """

    def __init__(self):
        pass

    def setDef(self,x,y):
        # x,y: for the white's binarizarion algorithm
        self.v = x
        self.b = y

    # returns the sum of (summed area) image pixels in the box between
    # (x1,y1) and (x2,y2)
    def mww(self,x1,y1,x2,y2,intImg):
        sum = intImg[x2][y2]
        if (x1>= 1 and y1 >= 1):
            sum = sum + intImg[x1-1][y1-1]
        if (x1 >= 1):
            sum = sum - intImg[x1-1][y2]
        if (y1 >= 1):
            sum = sum - intImg[x2][y1-1]
        return sum/((x2-x1+1)*(y2-y1+1))

    # constructs summed area table
    # img: original image
    # Nx: img size(x)
    # Ny: img size(y)
    # which: type of img: 
    #   'img': Image
    #   else:  array
    def sat(self,img,Nx,Ny,which):
        # summed area table, useful for speeding up the computing time by adding image pixels 
        intImg = np.zeros((Nx,Ny))
        img = np.array(img.getdata(),np.uint8).reshape(img.size[0], img.size[1])
        if(which == 'img'):
            intImg[0,0] = img[0,0]
            
            arrNx = range(1,Nx)
            arrNy = range(1,Ny)
            for h in arrNx:
                print img[h,0]
                intImg[h,0] = intImg[h-1,0] + img[h,0]
            
            for w in arrNy:
                intImg[0,w] = intImg[0,w-1] + img[0,w]
            
            for f in arrNx:
                for g in arrNy:
                   intImg[f,g] = img[f,g]+intImg[f-1,g]+intImg[f,g-1]-intImg[f-1,g-1]
        else:
            intImg[0,0] = img[0,0]
            
            arrNx = range(1,Nx)
            arrNy = range(1,Ny)
            for h in arrNx:
                intImg[h,0] = intImg[h-1,0] + img[h,0]
            
            for w in arrNy:
                intImg[0,w] = intImg[0,w-1] + img[0,w]
            
            for f in arrNx:
                for g in arrNy:
                   intImg[f,g] = img[f,g]+intImg[f-1,g]+intImg[f,g-1]-intImg[f-1,g-1]

        return intImg

    # white's algorithm
    # local threshold schema
    def white(self,img,Nx,Ny):
               
        im = np.zeros((Nx,Ny))
        
        intImg = self.sat(img,Nx,Ny,'img')
            
        vent = int(self.v)
        for i in range(Nx):
            for j in range(Ny):
                pix = img.getpixel((i,j))
                if(self.mww(max(0,i-vent),max(0,j-vent),min(Nx-1,i+vent),min(Ny-1,j+vent),intImg) >= pix*self.b ): 
                    im[i,j] = pix

        # do an opening operation to remove small elements
        return ndimage.binary_opening(im.T, structure=np.ones((1,1))).astype(np.int)

    def boxCount(self,e2,posx,posy,numBlocks,sx,sy,Nx,Ny):
        suma = 0
        for i in range(1,numBlocks+1):
            for j in range(1,numBlocks+1):

                xStart = posx+(i-1)*sx
                xEnd   = posx+i*sx - 1

                yStart = posy+(j-1)*sy
                yEnd   = posy+j*sy - 1

                dx = xEnd - (Nx-1) # sobrante en pixeles en x
                dy = yEnd - (Ny-1) # sobrante en pixeles en y
                if (dx > 0 and dy <= 0):
                    block1 = e2[xStart:Nx-1,yStart:yEnd]
                    block2 = e2[0:dx,yStart:yEnd]
                    suma+= np.sum(block1)+np.sum(block2) > 0

                if (dx <= 0 and dy > 0):
                    block1 = e2[xStart:xEnd,yStart:Ny-1]
                    block2 = e2[xStart:xEnd,0:dy]
                    suma+= np.sum(block1)+np.sum(block2) > 0

                if (dx > 0 and dy > 0):
                    block1 = e2[xStart:Nx-1,yStart:Ny-1]
                    block2 = e2[0:dx,yStart:Ny-1]
                    block3 = e2[xStart:Nx-1,0:dy]
                    block4 = e2[0:dx,0:dy]
                    suma+= np.sum(block1)+np.sum(block2)+np.sum(block3)+np.sum(block4) > 0

                if (dx <= 0 and dy <= 0): # todo el bloque esta en la grilla
                    block = e2[xStart:xEnd,yStart:yEnd]
                    #if(numBlocks > 
                    suma+= np.sum(block)>0 

        return suma

    def getFDs(self,filename):
        a = Image.open(filename)
        Nx, Ny = a.size
        L = Nx*Ny

        gray = a.convert('L') # rgb 2 gray

        IMG = True
        if(IMG):
            gray = self.white(gray,Nx,Ny) # local thresholding algorithm

            e2 = gray#np.array(gray.getdata(),np.uint8).reshape(gray.size[1], gray.size[0])
        else:
            e2 = np.array(gray.getdata(),np.uint8).reshape(gray.size[1], gray.size[0])

        plt.imshow(e2, cmap=matplotlib.cm.gray)
        plt.show()
        delta = []
        N = []
    
        for w in range(1,int(log(min(Nx,Ny))/log(2))):
            numBlocks = 2**w
            sx = np.floor(np.float32(Nx)/numBlocks)
            sy = np.floor(np.float32(Ny)/numBlocks)
            boxc = 0

            cant = 4#np.min((sx,8)).astype(int) # como maximo 16 (4x4) 
            suma = 0

            for c1 in range(cant): # promedio de casos
                posx = np.random.randint(1,sx)
                posy = np.random.randint(1,sy)
                print "POSS:", posx, posy
                temp = self.boxCount(e2,posx,posy,numBlocks,sx,sy,Nx,Ny)
                suma+=temp

                print "Proportion: ", temp, numBlocks*numBlocks
            boxc = np.floor(suma/cant)

            if(boxc > 0):
                delta.append(numBlocks) # numBlocks: (1/delta)
                N.append(boxc)

        
        x = np.log(delta)
        deltaA = np.vstack([x, np.ones(len(x))]).T
        print "delta:", x
        print "N:", np.log(N)
        m = np.linalg.lstsq(deltaA,np.log(N))

        fsize = 26
        plt.ylabel('$log(N_{\epsilon})$',fontsize=fsize)
        plt.xlabel('$log(1/\epsilon)$',fontsize=fsize)
        a = round(np.float32(str(m[0][0])),2)
        b = round(np.float32(str(m[1][0])),2)
        print np.array(x).shape
        print np.array(N).shape
        plt.plot(np.array(x),m[0][0]*np.array(x)+m[0][1],'r-',  label="Linear fit\nSlope (Box Dimension) = {0}\nR = {1}".format(a,np.float32(1.0)-b),linewidth=2.0)

        plt.plot(np.array(x),np.log(N),'bo',  label='Data')
        plt.legend(loc = 2)
        plt.show()

        return m
