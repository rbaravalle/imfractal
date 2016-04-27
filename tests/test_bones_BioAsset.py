"""
Copyright (c) 2016 Rodrigo Baravalle
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

#import Image
#import time
#import matplotlib.pyplot as plt
from pylab import *

import sys
sys.path.append('/almacen/imfractal/imfractal/Algorithm/')

import qs3D



def do_test():

    path = "/almacen/members.imaglabs.org/felix.thomsen/Rodrigo/BioAsset/mats/all/"

    patients = ["01"] #"5c", "6b", "8b", "8c", "V12"]
    scans = ["1"]#, "02", "03", "M1", "M2", "XCT"]
    # amount of volumes of interest
    vois = 1#27
    dims = 10

    # BioAsset bone's multifractal spectra database
    mfss = np.zeros([len(patients),len(scans),vois,2*dims+1])

    aux = CSandbox3D(dims)

    slices = "slice"
    masks = "mask"

    ii = 0
    for i in patients:
        jj = 0
        for j in scans:
            for k in range(1,vois+1):
                fmask =path + "M_" + masks + "BA" + i + "_120_" + j + ".mat"
                # set voi number and mask filename
                xct = "no"
                if(j == 5): xct = "xct"
                params = [1, 0.75, 3.7, 1, 15, k, fmask, xct, 'M', 'M']
                aux.setDef(40,1.02,True,params)

                filename = path + "M_" + slices + "BA" + i + "_120_" + j + ".mat"

                print i,j,"voi: ",k
                print fmask
                print filename

                mfss[ii,jj,k-1] = aux.getFDs(filename)

            jj = jj+1

        ii = ii+1


    np.save("mfs_BioAsset",mfss)


    
