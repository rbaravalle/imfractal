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

import sys
sys.path.append('/home/rodrigo/imfractal/imfractal/Algorithm/')

import qs3D



def do_test():

    # load array object file
    res = np.load("mfss.npy")

    patients = ["5c", "6b", "8b", "8c", "V12"]

    # scans except xct
    scans = ["M1", "M2", "01", "02", "03"]
    

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

    
