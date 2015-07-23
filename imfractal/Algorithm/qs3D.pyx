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

    
# sum of values in the region (x1,y1,z1), (x2,y2,z2) in intImg
# intImg: summed area table
def count(int x1,int y1,int z1, int x2,int y2, int z2, np.ndarray[DTYPE_ti, ndim=3] intImg, int Nx, int Ny, int Nz):
    cdef int sum
    if(x1 < 0): x1 = 0
    if(y1 < 0): y1 = 0
    if(z1 < 0): z1 = 0
    if(x2 >= Nx): x2 = Nx-1
    if(y2 >= Ny): y2 = Ny-1
    if(z2 >= Nz): z2 = Nz-1
    sum = intImg[x2,y2,z2]
    
    sum -= intImg[x2,y2,z1]
    sum -= intImg[x2,y1,z2]
    sum -= intImg[x1,y2,z2]
    sum += intImg[x1,y1,z2]
    sum += intImg[x1,y2,z1]
    sum += intImg[x2,y1,z1]
    sum -= intImg[x1,y1,z1]
    #if (x1>= 1 and y1 >= 1):
    #    sum += intImg[x1-1,y1-1]
    #if (x1 >= 1):
    #    sum -= intImg[x1-1,y2]
    #if (y1 >= 1):
    #    sum -= intImg[x2,y1-1]
    return sum

def aux(int P,int total,int Nx,int Ny,int Nz, np.ndarray[DTYPE_ti, ndim=2] points,np.ndarray[DTYPE_ti, ndim=3] intImg, int m0, int cant):

    cdef double summ, down
    cdef int i,x,y,z,MR, ind, R,h,q, stepR, startR

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
                    z = points[i,2]
                    MR = count(x-R,y-R,z-R,x+R,y+R,z+R,intImg,Nx,Ny,Nz)
                    summ+= down*pow(MR,np.double(q-1))

                summ /= float(total) # mean
                c[ind] = np.log(summ)/float(q-1)
                ind+=1
        else:
            #q = 1, < ln(M(R)/M0) >
            for R in rvalues:
                summ = 0.0
                for i from 0<=i<total:
                    x = points[i][0]
                    y = points[i][1]
                    z = points[i][2]
                    MR = count(x-R,y-R,z-R,x+R,y+R,z+R,intImg,Nx,Ny,Nz)
                    summ += np.log(MR/float(m0))
                summ /= float(total) # mean
                c[ind] = summ
                ind+=1

        slope, _, _, _, _ = scipy.stats.linregress(sizes,c)

        res[h] = slope
        h+=1

    return res


