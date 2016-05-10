import numpy as np
cimport numpy as np
import scipy
import scipy.stats
import time
import Image


DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t
ctypedef np.float32_t DTYPE_tf
ctypedef np.int32_t DTYPE_ti
ctypedef np.double_t DTYPE_td

cdef extern from "math.h":
    float sqrt(float x)
    float round(float x)
    double pow(int x,double y)
    double log(double x)

# returns the sum of (summed area) image pixels in the box between
# (x1,y1) and (x2,y2)
def mww(int x1,int y1,int x2,int y2,np.ndarray[DTYPE_ti, ndim=2] intImg, int Nx, int Ny):
    cdef double sum
    if(x1 < 0): x1 = 0
    if(y1 < 0): y1 = 0
    if(x2 >= Nx): x2 = Nx-1
    if(y2 >= Ny): y2 = Ny-1
    sum = np.double(intImg[x2,y2])
    if (x1>= 1 and y1 >= 1):
        sum = sum + intImg[x1-1][y1-1]
    if (x1 >= 1):
        sum = sum - intImg[x1-1][y2]
    if (y1 >= 1):
        sum = sum - intImg[x2][y1-1]
    return sum/((x2-x1+1)*(y2-y1+1))

def white(np.ndarray[DTYPE_ti, ndim=2] img,int Nx,int Ny,int v,float b):               
    cdef int i, j
    cdef np.ndarray[DTYPE_ti, ndim=2] im = np.zeros((Nx,Ny),dtype=np.int32)
    cdef np.ndarray[DTYPE_ti, ndim=2] intImg = sat(img,Nx,Ny)
    for i from 0<= i < Nx:
        for j from 0<= j < Ny:
            if(mww(i-v,j-v,i+v,j+v,intImg,Nx,Ny) >= img[i,j]*b and img[i,j] > 0):
                    im[i,j] = 255
    return im
               

# constructs summed area table
def sat(np.ndarray[DTYPE_ti, ndim=2] img,int Nx,int Ny):
    cdef np.ndarray[DTYPE_ti, ndim=2] intImg = np.empty((Nx,Ny),dtype = np.int32)
    cdef int f,g
        
    intImg[0,0] = img[0,0]
    intImg[1:,0] = intImg[0:-1,0] + img[1:,0]    
    intImg[0,1:] = intImg[0,0:-1] + img[0,1:]
    
    for f from 1<=f<Nx:
        for g from 1<=g<Ny:
           intImg[f,g] = img[f,g]+intImg[f-1,g]+intImg[f,g-1]-intImg[f-1,g-1]

    return intImg

# sum of values in the region (x1,y1), (x2,y2) in intImg
# intImg: summed area table
def count(int x1,int y1,int x2,int y2,np.ndarray[DTYPE_ti, ndim=2] intImg, int Nx, int Ny):
    cdef int sum
    if(x1 < 0): x1 = 0
    if(y1 < 0): y1 = 0
    if(x2 >= Nx): x2 = Nx-1
    if(y2 >= Ny): y2 = Ny-1
    sum = intImg[x2,y2]
    if (x1>= 1 and y1 >= 1):
        sum += intImg[x1-1,y1-1]
    if (x1 >= 1):
        sum -= intImg[x1-1,y2]
    if (y1 >= 1):
        sum -= intImg[x2,y1-1]
    return sum

def aux(int P, int total, int Nx, int Ny,
        np.ndarray[DTYPE_ti, ndim=2] points,
        np.ndarray[DTYPE_ti, ndim=2] intImg,
        int m0, int cant):

    cdef double summ, down
    cdef int i, x, y, MR, ind, R, h, q, stepR, startR

    stepR = 1
    startR = 1
    cdef np.ndarray[DTYPE_ti, ndim=1] rvalues = np.array(range(startR,P+startR,stepR)).astype(np.int32)
    # ln (R/L)
    cdef np.ndarray[DTYPE_td, ndim=1] sizes = np.log(rvalues/float(Nx))

    cdef np.ndarray[DTYPE_td, ndim=1] c = np.zeros((len(rvalues)), dtype=np.double )
    cdef np.ndarray[DTYPE_td, ndim=1] res = np.zeros((cant*2+1), dtype=np.double )


    h = 0
    for q from -cant <= q < cant+1:
        down = 1.0/pow(m0,np.double(q-1))
        ind = 0
        if(q != 1):
            # ln< M(R)/M0 ** q-1 >
            for R in rvalues:
                summ = 0.0
                for i from 0<=i<total:
                    x = points[i,0]
                    y = points[i,1]
                    MR = count(x-R,y-R,x+R,y+R,intImg,Nx,Ny)
                    summ+= down*pow(MR,np.double(q-1))

                summ /= float(total) # mean
                if(summ > 0.0):
                    c[ind] = np.log(summ)/float(q-1)
                ind+=1
        else:
            #q = 1, < ln(M(R)/M0) >
            for R in rvalues:
                summ = 0.0
                for i from 0<=i<total:
                    x = points[i][0]
                    y = points[i][1]
                    MR = count(x-R,y-R,x+R,y+R,intImg,Nx,Ny)
                    if(MR > 0):
                        summ += np.log(MR/float(m0))
                summ /= float(total) # mean
                c[ind] = summ
                ind+=1

        slope, _, _, _, _ = scipy.stats.linregress(sizes,c)

        res[h] = slope
        h+=1

    return res


