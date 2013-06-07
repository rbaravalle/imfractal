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
        # summed area table, useful for speed up the computation by adding image pixels 
        intImg = [ [ 0 for i in range(Nx) ] for j in range(Ny) ]
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
                    im[j,i] = img.getpixel((i,j))

        # do an opening operation to remove small elements
        return ndimage.binary_opening(im, structure=np.ones((2,2))).astype(np.int)

    def getFDs(self,filename):
        a = Image.open(filename)
        Nx, Ny = a.size
        L = Nx*Ny

        gray = a.convert('L') # rgb 2 gray

        gray = self.white(gray,Nx,Ny) # local thresholding algorithm

        e2 = gray #np.array(gray.getdata(),np.uint8).reshape(gray.size[1], gray.size[0])
        #plt.imshow(gray, cmap=matplotlib.cm.gray)
        #plt.show()
        delta = []
        N = []
    
        for i in range(9):
            numBlocks = 2**i
            sizeBlocks_x = np.floor(np.float32(Nx)/numBlocks)
            sizeBlocks_y = np.floor(np.float32(Nx)/numBlocks)

            flag = np.zeros((numBlocks,numBlocks))

            posx = []
            posy = []
            boxCount = np.float32(0)
            if (sizeBlocks_x > 1 and numBlocks > 1):
                cant = np.min((sizeBlocks_x,4)).astype(int) # como maximo 16 (4x4) 

                for j in range(cant):
                    val = np.floor(np.random.random()*sizeBlocks_x)
                    while(val in posx):
                        val = np.floor(np.random.random()*sizeBlocks_x)

                    posx.append(val)

                    val = np.floor(np.random.random()*sizeBlocks_x)
                    while(val in posy):
                        val = np.floor(np.random.random()*sizeBlocks_x)

                    posy.append(val)

                

                for c1 in range(cant): # promedio de casos
                    for c2 in range(cant):     # promedio de casos
                        for i in range(1,numBlocks):
                            for j in range(1,numBlocks):
                
                                xStart = posx[c1]+(i-1)*sizeBlocks_x
                                xEnd   = posx[c1]+i*sizeBlocks_x - 1

                                yStart = posy[c2]+(j-1)*sizeBlocks_y
                                yEnd   = posy[c2]+j*sizeBlocks_y - 1

                                dx = xEnd - (Nx-1) # sobrante en pixeles en x
                                dy = yEnd - (Ny-1) # sobrante en pixeles en y
                                if (dx > 0 and dy <= 0):
                                    block1 = e2[xStart:Nx-1][yStart:yEnd]
                                    block2 = e2[0:dx][yStart:yEnd]
                                    flag[i][j] = np.sum(block1) > 0 or np.sum(block2) > 0

                                if (dx <= 0 and dy > 0):
                                    block1 = e2[xStart:xEnd][yStart:Ny-1]
                                    block2 = e2[xStart:xEnd][0:dy]
                                    flag[i][j] = np.sum(block1) > 0 or np.sum(block2) > 0

                                if (dx > 0 and dy > 0):
                                    block1 = e2[xStart:Nx-1][yStart:Ny-1]
                                    block2 = e2[0:dx][yStart:Ny-1]
                                    block3 = e2[xStart:Nx-1][0:dy]
                                    block4 = e2[0:dx][0:dy]

                                    flag[i][j] = np.sum(block1) > 0 or np.sum(block2) > 0 or np.sum(block3) > 0 or np.sum(block4) > 0

                                if (dx <= 0 and dy <= 0): # todo el bloque esta en la grilla
                                    block = e2[xStart:xEnd][yStart:yEnd]
                                    flag[i][j] = np.sum(block) # mark this if ANY part of block is true


                            boxCount = boxCount + flag.sum()
                    boxCount = np.floor(boxCount/(cant*cant))
                else:
                    for i in range(1,numBlocks):
                        for j in range(1,numBlocks):
                            xStart = (i-1)*sizeBlocks_x
                            xEnd   = i*sizeBlocks_x - 1

                            yStart = (j-1)*sizeBlocks_y
                            yEnd   = j*sizeBlocks_y - 1

                            block = e2[xStart:xEnd][yStart:yEnd]
                            flag[i][j] = np.sum(block) > 0 # mark this if ANY part of block is true

                    boxCount = flag.sum()

                if(boxCount != 0):
                    delta.append(numBlocks) # numBlocks; (1/delta)
                    N.append(boxCount)

        x = np.log(delta)
        deltaA = np.vstack([x, np.ones(len(x))]).T
        m = np.linalg.lstsq(deltaA,np.log(N))

        fsize = 26
        plt.ylabel('$log(N_{\epsilon})$',fontsize=fsize)
        plt.xlabel('$log(1/\epsilon)$',fontsize=fsize)
        a = round(np.float32(str(m[0][0])),2)
        b = round(np.float32(str(m[1][0])),2)
        plt.plot(x,m[0][0]*x+m[0][1],'r-',  label="Linear fit\nSlope (Box Dimension) = {0}\nR = {1}".format(a,np.float32(1.0)-b),linewidth=2.0)
        plt.plot(x,np.log(N),'bo',  label='Data')
        plt.legend(loc = 2)
        plt.show()

        return m
