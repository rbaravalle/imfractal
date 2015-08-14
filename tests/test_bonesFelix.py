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

    # load array object file
    #res = np.load("mfss.npy")
    res = np.load("mfss.npy")


    patients = ["5c", "6b", "8b", "8c", "V12"]

    # scans except xct
    #scans = ["01", "02", "03","M1", "M2"]
    # all scans
    scans = ["01", "02", "03","M1", "M2", "xct"]


    dims = 10
    vois = 27


    print "test"
    with open("results.txt", "w") as text_file:
        for p in range(len(patients)):
            for k in range(len(scans)):
                for j in range(vois):      
                    s = str(p+1)+ " " + str(k+1)+ " " + str(j+1) + " "+  ' '.join(map(str, res[p][k][j])) + "\n"
                    text_file.write(s)
