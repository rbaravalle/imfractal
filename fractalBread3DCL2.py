import numpy as np
import random
import Image
import ImageDraw
import os
import time

#from lparams import *
from baking1D import calc
from mvc import mvc # mean value coordinates
import pyopencl as cl

N = 800
Nz = 100

def shape(x,y,z,r,field):
    r = np.round(r)
    rr = 0
    for i in range(x-r,x+r):
        for j in range(y-r,y+r):
            for k in range(z-r,z+r):
                if((x-i)*(x-i)+(y-j)*(y-j)+(z-k)*(z-k) < r*r):
                    if(i < N and i >= 0 and j < N and j >= 0 and k < Nz and k >= 0 ):
                        field[i,j,k] = rr

def f(c):
    return 20
    if(c==1): return 5
    if(c<=5): return 15
    if(c>5): return 7

def paint(i,N,field3,v,p):
    for x in range(int(i[0])-p,int(i[0])+p):
        for y in range(int(i[1])-p,int(i[1])+p):
            #for z in range(int(i[2])-p,int(i[2])+p):
            field3[np.clip(x,0,N-1),np.clip(y,0,N-1),24] = v
            field3[np.clip(x,0,N-1),np.clip(y,0,N-1),25] = v
            field3[np.clip(x,0,N-1),np.clip(y,0,N-1),26] = v
    return field3

def computeCages():
    p = int(60)
    # Arbitrary shape, mean value coordinates
    #cageOrig = np.array([[p,N-1-p],[N/2,N-1-p],[N/2+40,N-1-p],[N-1-p,N-1-p],[N-1-p,N/2+40],[N-1-p,N/2],[N-1-p,N/2-40],[N-1-p,p],[N/2+40,p],[N/2,p],[N/2-40,p],[p,p]]).astype(np.float32)
    cageOrig = np.array([[p,N-1-p],[p+50,N-1-p],[N/2,N-1-p],[N/2+50,N-1-p],[N-1-p,N-1-p],[N-1-p,N/2+50],[N-1-p,N/2],[N-1-p,N/2-50],[N-1-p,p],[N/2+50,p],[N/2,p],[N/2-50,p],[p,p],[p,N/2-50],[p,N/2],[p,N/2+50]]).astype(np.float32)

    cageReal = np.array(cageOrig)
    cageNew = np.array(cageOrig)

    # control points displacements
    # RIGHT, BOTTOM,LEFT, TOP < |> <_> <| > <->
    # X: BOTTOM - TOP
    # Y : LEFT - RIGHT
    #trs=[[0,-20],[0,0],[0,0],[0,0],[30,0],[45,0],[30,40],[-5,0],[30,10],[25,10],[45,20],[12,15]]
    #trs=[[0,-10],[0,0],[0,0],[0,0],[20,0],[35,0],[20,30],[-5,0],[20,10],[15,10],[25,10],[12,5]]
    defor = 20
    trs=[[0,10+defor],[0,25],[0,10+2*defor],[-10,2*defor+3],[-60,-10+defor],[12,2*defor],[13,10+2*defor],[14,-10+2*defor],[-60,defor],[0,-defor-5],[0,-defor],[0,-defor-5],[defor,defor],[-15,0],[-20,2*defor],[-15,defor]]

    for i in range(len(cageOrig)):
        cageReal[i] = cageOrig[i]+trs[i]
        cageNew[i] = cageOrig[i]-trs[i]

    print len(cageNew),len(cageOrig)

    #print cageOrig
    #print cageReal
    #print cageNew

    return cageReal, np.array(map (lambda i: np.array(i),cageNew)), cageOrig


def main(param_a,param_b,param_c,param_d,param_e):
    #print "Starting..."

    if not os.path.isdir('warp'): 
        os.mkdir ( 'warp' ) 

    if not os.path.isdir('warp/baked'): 
        os.mkdir ( 'warp/baked' ) 

    if not os.path.isdir('warp/warped'): 
        os.mkdir ( 'warp/warped' ) 

    arr = calc()
    gx, gy = np.gradient(arr)
    field = np.zeros((N,N,Nz)).astype(np.uint8) + np.uint8(255)

    #global v
    #maxx = 100000

    gx2 = np.zeros((N,N)) # vector field
    gy2 = np.zeros((N,N)) # vector field

    shx = gx.shape[0]-1
    shy = gy.shape[1]-1

    #print "Computing vector field from baking..."
    # the same foreach z
    for i in range(0,N):
        for j in range(0,N):
            dist = np.sqrt(((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2)))
            u = i*(shx/(np.float32(N)-1))
            v = j*(shy/(np.float32(N)-1))
            gx2[i,j] = gx[u,v]*dist
            gy2[i,j] = gy[u,v]*dist

    cageReal,cageNew,cageOrig = computeCages()

    I = Image.new('L',(N,Nz*N),0.0)
    #print "Proving..."
    for i in range(param_d,param_e,param_a):
        cubr = param_b/float(20.0)
        maxrank = cubr*N*N*Nz/(np.power(i,param_c))
        if(maxrank >=1):
            #print i, maxrank
            for j in range(0,int(np.floor(maxrank))):
                shape(np.random.randint(0,N),np.random.randint(0,N),np.random.randint(0,Nz),i,field)
        else: break

        if(i%5==0):
            #print "Saving proving..."
            for w in range(Nz):
                III = Image.frombuffer('L',(N,N), np.array(field[:,:,w]).astype(np.uint8),'raw','L',0,1)
                I.paste(III,(0,N*w))

            #print 'warp/proving'+str(i)+'.png'
            I.save('warp/proving'+str(i)+'.png')

    field2 = np.zeros((N,N,Nz)).astype(np.uint8) + np.uint8(255)
    k = float(15.0)


    #print "Baking..."
    for x in range(0,N):    
        for y in range(0,N):
            u = np.round(x+k*gx2[x,y])
            v = np.round(y+k*gy2[x,y])
            if(not(u < 0 or u >= N or v < 0 or v >= N)):
                field2[x,y] = field[u,v]


    Ires = 0
    I3 = Image.new('L',(N/2,(N/2)*(Nz/2)),0.0)
    for w in range(Nz):
        I2 = Image.frombuffer('L',(N/2,N/2), 255-np.array(field2[N/4:N*3/4,N/4:N*3/4,w]).astype(np.uint8),'raw','L',0,1)


        I2.save('warp/baked/bakedslice'+str(w)+'.png')
        I3.paste(I2,(0,N*w))

    #print 'warp/baked'+str(i)+'.png'
    I3.save('warp/baked'+str(i)+'.png')

    #return I2

    field3 = np.zeros((N,N,Nz)).astype(np.uint8) + np.uint8(255)
    print "Warping..."

    t2 = t1 = t3= 0
    p=int(N/4)
    eps = 10.0*np.nextafter(0,1)

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    mf = cl.mem_flags

    nSize = len(cageOrig)


    prg = cl.Program(ctx, """
    __kernel void main( __global int *cageOrig, __global int *cageNew, __global int *xn, const int nSize, const float eps, const int p, const int N) {
        int x = get_global_id(0)+p;
        int y = get_global_id(1)+p;

        float dest[17], dest2[17], dx, dy;
        int s[17*2];

        int i,h;

        for(i = 0; i < nSize; i++) {
            dest[i]   =   0;
            dest2[i]  =   0;
            s[2*i]    =   cageOrig[2*i]   - x;
            s[2*i+1]  =   cageOrig[2*i+1] - y;
        }

        int cut = 0; // FIX ME!!: we should skip one part or not 
        for(i = 0; i < nSize; i++) {
            dx  =   (float)s[2*i];
            dy  =   (float)s[2*i+1];
            int ip = (i+1)%nSize;
            float sC = s[2*ip];
            float sD = s[2*ip+1];

            float ri = sqrt( dx*dx + dy*dy );
            float Ai = 0.5*(dx*sD - sC*dy);
            float Di = sC*dx + sD*dy;


            if(ri <= eps) {dest[i] = 1.0; cut = 1; break;}
            else {
                if( fabs((float)Ai) == 0.0 && Di < 0.0){
                    dx = (float)(cageOrig[2*ip] - cageOrig[2*i]);
                    dy = (float)(cageOrig[2*ip+1] - cageOrig[2*i+1]);
                    float dl = sqrt( dx*dx + dy*dy);
                    dx = (float)(x - cageOrig[2*i]);
                    dy = (float)(y - cageOrig[2*i+1]);
                    float mu = sqrt(dx*dx + dy*dy)/dl;
                    dest[i]  = 1.0-mu;
                    dest[ip] = mu;
                    cut = 1;
                    break;
                }
            }

            float rp = sqrt( sC*sC + sD*sD );
            if(Ai == 0.0) dest2[i] = 0.0;
            else dest2[i] = (ri*rp - Di)/(2.0*Ai);
        }

        if(cut == 0) {
            float wsum = 0.0;
            for(i = 0; i < nSize; i++) {
                dx  =   (float)(cageOrig[2*i]   - x);
                dy  =   (float)(cageOrig[2*i+1] - y);
                float ri = sqrt( dx*dx + dy*dy );
                int im = (nSize-1+i)%nSize;
                float wi = 2.0*( dest2[i] + dest2[im] )/ri;
                dest[i] = wi;
                wsum += wi;
            }


            if( fabs(wsum) > 0.0)
                for(i = 0; i < nSize; i++) 
                    dest[i] /= wsum;
        }

        // dest : computed barycoords
        float msumx = 0.0;
        float msumy = 0.0;
        for(i = 0; i < nSize; i++) {
            msumx += dest[i]*cageNew[2*i];
            msumy += dest[i]*cageNew[2*i+1];
        }

        int isumx = (int)msumx;
        int isumy = (int)msumy;

        if(isumx > N-1) isumx = N-1;
        if(isumy > N-1) isumy = N-1;
        if(isumx < 0) isumx = 0;
        if(isumy < 0) isumy = 0;

        xn[x+y*(N-2*p)] = isumx+isumy*(N-2*p);

    }
    """).build()

    a = np.array(cageOrig).astype(np.int32)
    cageOrig_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a)

    a = np.array(cageNew).astype(np.int32)
    cageNew_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a)

    s = np.zeros((nSize,2)).astype(np.float32)
    s_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a)

    p = 0
    t = time.clock()

    xn = np.zeros(((N-2*p)*(N-2*p))).astype(np.int32)
    d2 = np.zeros((nSize)).astype(np.float32)
    d = np.zeros((nSize)).astype(np.float32)
    dest_xn_buf = cl.Buffer(ctx, mf.WRITE_ONLY, xn.nbytes)

    prg.main(queue, (N-2*p,N-2*p), None, cageOrig_buf, cageNew_buf, dest_xn_buf, np.int32(nSize), np.float32(eps), np.int32(p), np.int32(N))
    cl.enqueue_read_buffer(queue, dest_xn_buf, xn).wait()

    print "TIEMPO MVC: ", time.clock()-t

    t = time.clock()
    for x in range(p,N-p):
       for y in range(p,N-p):
           xxn = xn[x+y*(N-2*p)]
           a = xxn%(N-2*p)
           b = xxn/(N-2*p)
           field3[x,y] = field2[a+p,b+p]

    print "TIEMPO MVC2: ", time.clock()-t

    for h in range(len(cageReal)):
        field3 = paint(cageReal[h],N,field3,0,4)
        field3 = paint(cageOrig[h],N,field3,120,6)
    #I3 = Image.new('L',(N,Nz*N),0.0)
    print "Saving Image..."
    Ires = 0
    for w in range(Nz):
        II = Image.frombuffer('L',(N,N), 255-np.array(field3[:,:,w]).astype(np.uint8),'raw','L',0,1)
        if(w >= 20 and w <= 29):
            II.save('images/test4/bread/warpedslice'+str(w)+'.png')
            Ires = II
        II.save('warp/warped/warpedslice'+str(w)+'.png')
        I3.paste(II,(0,N*w))
        #II = Image.frombuffer('L',(N,N), np.array(field3[:,:,w]).astype(np.uint8),'raw','L',0,1)
        #I3.paste(II,(0,N*w))

    I3.save('warp/warped'+str(i)+'.png')
    print "Image "+ str(i) +" saved"
    return Ires


main(1,0.06,2.7,2,20)
