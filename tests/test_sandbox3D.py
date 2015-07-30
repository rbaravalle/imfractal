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
    filename = 'images/train2/bread/b1.png'
    #filename = 'images/warpbin.png'

    dims = 10
    
    params = np.array([1,0.75,3.7,1,15]).astype(np.float32)
    #arr = qs3D.volume(params,256,256)

    i = CSandbox3D(dims)
    i.setDef(40,1.02,True,params)
    fds2 = i.getFDs(filename)
    arr2 = i.readDicom("/home/rodrigo/dicom")

    for h in range(2):
        print "Computing 3D Cython Sandbox Multifractal Spectrum..."
        t =  time.clock()
        fds2 = np.vstack((fds2,i.getFDs(filename)))
        t =  time.clock()-t
    plt.figure(1)
    plt.subplot(121)
    plt.ylim(ymax = 3.8, ymin = 2.4)
    plt.title(str(params))
    plt.boxplot(fds2, sym='')#, 'b*', label='synthetic',linewidth=2.0)
    plt.subplot(122)

    plt.imshow(arr2[:,:,220],plt.get_cmap('gray'))
    plt.show()

    
