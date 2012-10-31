from Algorithm import *
from random import randrange,randint
from math import log
from scipy import ndimage
import time
import Image
import numpy as np
import sys
import os

class Sandbox (Algorithm):

    """

    :sandbox multifractal spectrum
    :version: 1.0
    :author: Rodrigo Baravalle
    """

    def setDef(self,x,y):
        self.total = 10*10      # number of pixels for averaging
        self.P = 40             # window
        self.cant = 20+1           # number of fractal dimensions (-1 x2)
        self.v = x
        self.b = y

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
                
    # get multifractal dimentions
    def getFDs(self,filename):
        t = time.clock()
        tP = (self.P)   # tP : two raised to P
        x = tP+1
        y = tP+1
        cantSelected = 0

        a = Image.open(filename)
        Nx, Ny = a.size
        L = Nx*Ny

        points = []     # number of elements in the structure
        gray = a.convert('L') # rgb 2 gray

        gray = self.white(gray,Nx,Ny) # local thresholding algorithm
        #plt.imshow(gray, cmap=matplotlib.cm.gray)
        #plt.show()

        intImg = self.sat(gray,Nx,Ny,'array')

        m0 = intImg[Nx-1][Ny-1]#/float(L)

        while(gray[x][y] == 0):
            x = randint(tP,Nx-tP-1)
            y = randint(tP,Ny-tP-1)

        # list with selected points (the points should be in the "structure")
        # points shouldn't be close to the borders, in order for the windows to have the same size
        while cantSelected < self.total:
            while(([x,y] in points) or self.count(x-1,y-1,x+1,y+1,intImg) == 0):
                x = randint(tP,Nx-tP-1)
                y = randint(tP,Ny-tP-1)
            # new point, add to list
            points.append([x,y])
            cantSelected = cantSelected+1


        c = np.zeros((self.total+1,self.P), dtype=np.double ) # total+1 rows x P columns
        for i in range(self.total): # for each point randomly selected
            x = points[i][0]
            y = points[i][1]
            for h in range(1,self.P+1):
                # how many points in the box. M(R) in the literature
                c[i+1][h-1] = self.count(x-(h),y-(h),x+(h),y+(h),intImg)

        down = range(1,self.P+1)
        # Generalized Multifractal Dimentions 
        s = [0 for i in range(2*self.cant-2)]
        l = range(-self.cant+1,0)+  range(1,self.cant)
        j = 0
        for i in l:
            s[j] = self.Dq(c,i,L,m0,down)

            j = j+1

        #plot(s)
        #show()

        t =  time.clock()-t
        print "Time: ", t
        return s

    # D_{q} 
    def Dq(self,c,q,L,m0,down):
        # sum in each radius, all the points
        if q>0: # math representation issue
            aux1 = float(self.total)
            aux2 = 0
        else:
            aux1 = 1
            aux2 = log(float(self.total))

        for h in range(1,self.P+1):        
            for i in range(self.total):
                c[0][h-1] = c[0][h-1] + ((c[i+1][h-1]**q)/aux1) # mean of "total" points

        up = map(lambda i: log(i)-aux2-q*log(m0), c[0])
        up2 = map(lambda i: up[i]/(q*down[i]), range(len(up)))
        sizes = range(1,self.P+1)

        (ar,br)=np.polyfit(sizes,up2,1)
        return ar
        

