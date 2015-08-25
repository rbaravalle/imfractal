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
from imfractal import *

import Image
import time
import matplotlib.pyplot as plt
from pylab import *

import os
import sys
sys.path.append('/home/rodrigo/imfractal/imfractal/Algorithm/')

import qs3D

def createFolders(dirs):

    for f in dirs:
        if not os.path.isdir(f): 
            os.mkdir (f) 


def do_test():



    arr = [ "imfractal/Algorithm/qs3D", "imfractal/Algorithm/qs"]

    print "WTF"

    for i in range(len(arr)):

        command1 = "cython "+arr[i]+".pyx "
        command2 = "gcc -c -fPIC -I/usr/include/python2.7/ "+arr[i]+".c"+" -o "+arr[i]+".o"
        command3 = "gcc -shared "+arr[i]+".o -o "+arr[i]+".so"

        print command1
        os.system(command1)
        print command2
        os.system(command2)
        print command3
        os.system(command3)


    # load array object file
    #res = np.load("mfss.npy")
    res = np.load("mfss.npy")
    #res = np.load("mfssM90-100-P20-DIFF.npy")
    #res = np.load("mfssM90-100-P31-DIFF.npy")

    patients = ["5c", "6b", "8b", "8c", "V12"]

    # scans except xct
    scans = ["01", "02", "03","M1", "M2"]
    # all scans
    #scans = ["01", "02", "03","M1", "M2", "xct"]


    dims = 10
    vois = 27

    # compute R2
    # res = [patient, scan, voi] -> MFS
    # compute spj for j in voi's and i scans of patient p and for each mfdimension
    for p in range(len(patients)):
        sp = np.zeros((vois,2*dims+1)).astype(np.float32)
        spg = np.zeros(2*dims+1).astype(np.float32)
        for j in range(vois):      
           temp = np.zeros(2*dims+1).astype(np.float32)
           for k in range(len(scans)):
                temp += res[p][k][j]
                spg += res[p][k][j]

           sp[j] = temp/len(scans)
        spg = spg/(len(scans)*vois)

        debug = False
        if(debug):
            print sp

            plt.plot(spg)
            plt.show()
            for i in range(27):
                plt.plot(sp[i])
            plt.show()


        r2 = np.zeros(2*dims+1)
        # R2 for each dimension
        for r in range(len(r2)):
            
            sumup = 0
            sumdown = 0
        
            for j in range(vois):
                for s in range(len(scans)):
                    # sumup
                    temp = res[p][s][j][r] - sp[j][r]
                    sumup = sumup + temp*temp

                    # sumdown
                    temp = res[p][s][j][r] - spg[r]
                    sumdown = sumdown + temp*temp

            r2[r] = 1.0 - sumup / sumdown

        print r2
    

    exit()
    

    pp = 0
    for p in patients:
        ss = 0
        for s in scans:
            for i in range(1,26,2):
                clf()
                plt.figure(pp*len(patients)+ss*len(scans)+i+1)
                plt.subplot(221)
                plt.ylim(ymax = 3.6, ymin = 2.0)
                #plt.title("XCT 5c_XtremeCTSlices")
                plt.title(p + "-" + s + " - HRCT - VOI "+str(i))
                plt.plot(res[pp][ss][i])

                plt.subplot(222)
                plt.ylim(ymax = 3.6, ymin = 2.0)
                plt.title(p + " - XCT - VOI "+str(i))
                plt.plot(res[pp][5][i])

                plt.subplot(223)
                plt.ylim(ymax = 3.6, ymin = 2.0)
                
                plt.title(p + "-" + s + " - HRCT - VOI "+ str(i+1))
                plt.plot(res[pp][ss][i+1])

                plt.subplot(224)
                plt.ylim(ymax = 3.6, ymin = 2.0)
                plt.title(p + " - XCT - VOI "+ str(i+1))
                plt.plot(res[pp][5][i+1])

                
                savefig("exps/figs/"+p+s+"VOI"+str(i)+"-"+str(i+1)+".png")

            ss = ss+1
        pp = pp+1

    
