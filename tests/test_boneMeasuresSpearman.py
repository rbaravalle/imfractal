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
import scipy

def createFolders(dirs):

    for f in dirs:
        if not os.path.isdir(f): 
            os.mkdir (f) 

def do_test():

    #path = "/home/rodrigo/rodrigo/rodrigo/boneMeasures/SiCopyRodrigo/Win32Release/BVTV/"

    measures = np.load('measures.npy')
    mfss = np.load('mfss.npy')

    dims = 10  

    patients = ["5c", "6b", "8b", "8c", "V12"]
    scans = ["01", "02", "03", "M1", "M2", "XCT"]
    
    # amount of volumes of interest (voi)
    vois = 27
    fsize = 18
    e = 0.1
    s = 10.0
    x1 = 1
    x2 = 20
    x = np.arange(2*dims+1)+1



    for m in range(measures.shape[3]):


        correls = np.zeros((len(patients),vois,2*dims+1))
        correlsAll = np.zeros((len(patients),2*dims+1))
        for p in range(len(patients)):
            for j in range(vois):      
               for i in range(2*dims+1):
                   correls[p,j,i] = scipy.stats.stats.spearmanr(mfss[p,:,j,i],measures[p,:,j,m])[0]

        # for all vois
        for p in range(len(patients)):
            for i in range(2*dims+1):
                allmeasures = measures[p,:,0,m]
                allvois = mfss[p,:,0,i]
                for j in range(1,vois):
                   allvois = np.hstack((allvois,mfss[p,:,j,i]))
                   allmeasures = np.hstack((allmeasures,measures[p,:,j,m]))

                print allvois.shape, allmeasures.shape
                correlsAll[p,i] = scipy.stats.stats.spearmanr(allvois,allmeasures)[0]

        print correlsAll

        plt.ylim((-1, 1))
        plt.xlim(x1,x2)
        plt.title('ALL - Measure ' + str(m))
        plt.ylabel(r'$\rho$',fontsize=fsize)
        plt.xlabel(r'$q$',fontsize=fsize)   
        plt.plot(x, correlsAll[0], 'kD--', label='5c', linewidth=2.0,markeredgewidth=e, markersize=s)
        plt.plot(x, correlsAll[1], 'rs--', label='6b', linewidth=2.0,markeredgewidth=e, markersize=s)
        plt.plot(x, correlsAll[2], 'b^--', label='8b', linewidth=2.0,markeredgewidth=e, markersize=s)
        plt.plot(x, correlsAll[3], 'r*--', label='8c', linewidth=2.0,markeredgewidth=e, markersize=s)
        plt.plot(x, correlsAll[4], 'gv--', label='V12', linewidth=2.0,markeredgewidth=e, markersize=s)
        plt.legend(loc = 2) # loc 4: bottom, right
        plt.show()
        
        print correls


        figure(2)

        if(False):
            xticks(x,range(-10,10)) # translate
            for j in range(vois):
                plt.ylim((-1, 1))
                plt.xlim(x1,x2)
                plt.title('VOI ' + str(j+1) + ' Measure ' + str(m))
                plt.ylabel(r'$\rho$',fontsize=fsize)
                plt.xlabel(r'$q$',fontsize=fsize)   
                #for p in range(len(patients)): 
                    #plt.plot(x, correls[p,j], 'kD--', linewidth=1.5, markeredgewidth=e, markersize=s)
                plt.plot(x, correls[0,j], 'kD--', label='5c', linewidth=2.0,markeredgewidth=e, markersize=s)
                plt.plot(x, correls[1,j], 'rs--', label='6b', linewidth=2.0,markeredgewidth=e, markersize=s)
                plt.plot(x, correls[2,j], 'b^--', label='8b', linewidth=2.0,markeredgewidth=e, markersize=s)
                plt.plot(x, correls[3,j], 'r*--', label='8c', linewidth=2.0,markeredgewidth=e, markersize=s)
                plt.plot(x, correls[4,j], 'gv--', label='V12', linewidth=2.0,markeredgewidth=e, markersize=s)
                plt.legend(loc = 3) # loc 4: bottom, right
                plt.show()




    np.save("correlsMeasures",correls)

    
