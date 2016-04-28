import numpy as np
cimport numpy as np
import scipy
import scipy.stats
import time
import Image
import matplotlib.pyplot as plt
from numpy import linalg

from libc.stdlib cimport rand

cdef extern from "limits.h":
    int INT_MAX

cdef extern from "math.h":
    float pow(int x ,float y)
    int floor(float x)


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

cdef saveField(field,filename):
    pp = 0
    N = field.shape[0]
    Nz = N
    I3 = Image.new('L',(N-2*pp,(N-2*pp)*(Nz)),0.0)

    for w in range(Nz):
        II = Image.frombuffer('L',(N-2*pp,N-2*pp), np.array(field[pp:N-pp,pp:N-pp,w]).astype(np.uint8),'raw','L',0,1)
        if(w == 216 and filename == "textures/system"): 
            II.save(filename+"/slice216.png")
        I3.paste(II,(0,(N-2*pp)*w))

    I3.save(filename+".png")
    print "Image "+filename+" saved"


def volume(np.ndarray[DTYPE_tf, ndim=1] params, int N, int Nz):
    cdef int param_a = int(params[0])
    cdef float param_b = params[1]
    cdef float param_c = params[2]
    cdef int param_d = int(params[3])
    cdef int param_e = int(params[4])

    print param_a,param_b,param_c,param_d,param_e
    cdef int r, v, i, j, k,x,y,z, maxrank
    cdef float cubr,rr
    cdef np.ndarray[DTYPE_t, ndim=3] field = np.zeros((N,N,Nz),dtype=DTYPE) + np.uint8(1)
    cubr = (param_b/float(20.0))*N*N*Nz
    for r from param_d <= r < param_e by param_a:
        maxrank = floor(cubr/(pow(r,param_c)))
        if(maxrank >=1.0):
            for v from 0<=v< maxrank:
                x = floor((rand() / float(INT_MAX))*N)
                y = floor((rand() / float(INT_MAX))*N)
                z = floor((rand() / float(INT_MAX))*Nz)
                rr = 4.0*np.random.random()
                for i from x-r<=i<=x+r:
                    for j from y-r<=j<=y+r:
                        for k from z-r<=k<=z+r:
                            if((x-i)*(x-i)+(y-j)*(y-j)+(z-k)*(z-k) <= r*r):
                                if(i < N and i >= 0 and j < N and j >= 0 and k < Nz and k >= 0 ):
                                    field[i,j,k] = 0

    saveField(255*field,'spec')
    return field
    
# sum of values in the region (x1,y1,z1), (x2,y2,z2) in intImg
# intImg: summed area table
cdef count(int x1,int y1,int z1, int x2,int y2, int z2, np.ndarray[DTYPE_ti, ndim=3] intImg, int Nx, int Ny, int Nz):
    cdef int sum
    # all up to 2,2,2
    sum = intImg[x2,y2,z2]
    # we take out three prisms
    sum -= intImg[x2,y2,z1]
    sum -= intImg[x2,y1,z2]
    sum -= intImg[x1,y2,z2]
    # we took out three extras, so we add them
    sum += intImg[x1,y1,z2]
    sum += intImg[x1,y2,z1]
    sum += intImg[x2,y1,z1]
    # we added one extra, so we take it out
    sum -= intImg[x1,y1,z1]
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
        down = 1.0/np.power(m0,np.double(q-1))
        ind = 0

        if(q != 1):
            for R in rvalues:
                # ln< M(R)/M0 ** q-1 >
                summ = 0.0
                for i from 0<=i<total:
                    x,y,z = points[i]
                    MR = count(x-R,y-R,z-R,x+R,y+R,z+R,intImg,Nx,Ny,Nz)
                    summ+= down*np.power(MR,np.double(q-1))

                summ /= np.float32(total) # mean
                c[ind] = np.log(summ)/np.float32(q-1)
                ind+=1

        slope, _, _, _, _ = scipy.stats.linregress(sizes,c)

        #if(q%5==0):
            #plt.plot(sizes,c,'o', label = str(q))

        print slope, " , q: ", q

        res[h] = slope
        h+=1

    #plt.legend(loc=2)
    #plt.show()


    res[cant+1] = (res[cant]+res[cant+2])/2.0

    return res


